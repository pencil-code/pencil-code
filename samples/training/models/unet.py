import torch
import torch.nn as nn
import torch.nn.functional as F

class Normalizer(torch.nn.Module):
    def __init__(self, size, name, max_accumulation=10**8, std_epsilon=1e-8):
        super(Normalizer, self).__init__()
        self.name = name
        self.max_accumulation = max_accumulation
        self.std_epsilon = std_epsilon

        self.register_buffer('acc_count', torch.zeros(1, dtype=torch.float32))
        self.register_buffer('num_acc', torch.zeros(1, dtype=torch.float32))
        self.register_buffer('acc_sum', torch.zeros(size, dtype=torch.float32))
        self.register_buffer('acc_sum_squared', torch.zeros(size, dtype=torch.float32))

        
    def forward(self, batched_data: torch.Tensor, accumulate: bool) -> torch.Tensor:
        """Normalizes input data and accumulates statistics."""
        if accumulate and self.num_acc < self.max_accumulation:
            self._accumulate(batched_data)
        batched_data = (batched_data - self._mean()) / self._std_with_eps()
        return batched_data
    
    def inverse(self, normalized_batch_data):
        """Inverse transformation of the normalizer."""
        return normalized_batch_data * self._std_with_eps() + self._mean()
        
    def _accumulate(self, batched_data):
        data_sum = torch.sum(batched_data, dim=[2,3,4], keepdim=True)
        data_sum_squared = torch.sum(batched_data**2, dim=[2,3,4], keepdim=True)

        self.acc_sum += data_sum
        self.acc_sum_squared += data_sum_squared
        self.acc_count += torch.tensor(batched_data[0,0,:,:,:].numel()).to(self.num_acc.device)
        self.num_acc += torch.tensor(1.0).to(self.num_acc.device)
        
    def _mean(self):
        safe_count = torch.maximum(self.acc_count, torch.tensor(1.0).to(self.acc_count.device))
        return self.acc_sum / safe_count
    
    def _std_with_eps(self):
        safe_count = torch.maximum(self.acc_count, torch.tensor(1.0).to(self.acc_count.device))
        std = torch.sqrt(self.acc_sum_squared / safe_count - self._mean()**2)
        return torch.maximum(std, torch.tensor(self.std_epsilon).to(self.acc_count.device))
    


class UNet3D(nn.Module):
    def __init__(self, in_channels=3, out_channels=6, init_features=16):
        super(UNet3D, self).__init__()
        features = init_features

        self.normalizer = Normalizer(size=[1, in_channels, 1, 1, 1], name="input")

        self.encoder1 = UNet3D._block(in_channels, features)
        self.pool1 = nn.MaxPool3d(kernel_size=2, stride=2)

        self.encoder2 = UNet3D._block(features, features * 2)
        self.pool2 = nn.MaxPool3d(kernel_size=2, stride=2)

        self.bottleneck = UNet3D._block(features * 2, features * 4)

        self.upconv2 = nn.ConvTranspose3d(features * 4, features * 2, kernel_size=2, stride=2)
        self.decoder2 = UNet3D._block(features * 4, features * 2)

        self.upconv1 = nn.ConvTranspose3d(features * 2, features, kernel_size=2, stride=2)
        self.decoder1 = UNet3D._block(features * 2, features)

        self.conv = nn.Conv3d(features, out_channels, kernel_size=1)

    def forward(self, x):
        x = self.normalizer(x, accumulate=True).float()
        enc1 = self.encoder1(x)
        enc2 = self.encoder2(self.pool1(enc1))

        bottleneck = self.bottleneck(self.pool2(enc2))

        dec2 = self.upconv2(bottleneck)
        
        dec2 = F.interpolate(dec2, size=enc2.shape[2:], mode='trilinear', align_corners=True)
        dec2 = torch.cat((dec2, enc2), dim=1)
        dec2 = self.decoder2(dec2)

        dec1 = self.upconv1(dec2)
        
        dec1 = F.interpolate(dec1, size=enc1.shape[2:], mode='trilinear', align_corners=True)
        dec1 = torch.cat((dec1, enc1), dim=1)
        dec1 = self.decoder1(dec1)

        return self.conv(dec1).double()


    
    def _block(in_channels, out_channels):
        return nn.Sequential(
            nn.Conv3d(in_channels, out_channels, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
            nn.Conv3d(out_channels, out_channels, kernel_size=3, padding=1),
            nn.ReLU(inplace=True)
        )
    

class CNN3D(nn.Module):
    def __init__(self, in_channels=3, out_channels=6, init_features=16):
        super(CNN3D, self).__init__()
        
        self.normalizer = Normalizer(size=[1, in_channels, 1, 1, 1], name="input")
        self.layer = nn.Sequential(
            nn.Conv3d(in_channels, 12, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
            nn.Conv3d(12, 24, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
            nn.Conv3d(24, out_channels, kernel_size=3, padding=1)
        )

    def forward(self, x):
        return self.layer(self.normalizer(x, accumulate=True).float()).double()
