import torch
import torch.nn as nn
import torch.nn.functional as F
import os

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


class SpectralConv3d(nn.Module):
    def __init__(self, in_channels, out_channels, modes1, modes2, modes3):
        super(SpectralConv3d, self).__init__()

        """
        3D Fourier layer. It does FFT, linear transform, and Inverse FFT.    
        """

        self.in_channels = in_channels
        self.out_channels = out_channels
        self.modes1 = modes1 #Number of Fourier modes to multiply, at most floor(N/2) + 1
        self.modes2 = modes2
        self.modes3 = modes3

        self.scale = (1 / (in_channels * out_channels))
        self.weights1 = nn.Parameter(self.scale * torch.rand(in_channels, out_channels, self.modes1, self.modes2, self.modes3, dtype=torch.cfloat))
        self.weights2 = nn.Parameter(self.scale * torch.rand(in_channels, out_channels, self.modes1, self.modes2, self.modes3, dtype=torch.cfloat))
        self.weights3 = nn.Parameter(self.scale * torch.rand(in_channels, out_channels, self.modes1, self.modes2, self.modes3, dtype=torch.cfloat))
        self.weights4 = nn.Parameter(self.scale * torch.rand(in_channels, out_channels, self.modes1, self.modes2, self.modes3, dtype=torch.cfloat))

    # Complex multiplication
    def compl_mul3d(self, input, weights):
        # (batch, in_channel, x,y,t ), (in_channel, out_channel, x,y,t) -> (batch, out_channel, x,y,t)
        return torch.einsum("bixyz,ioxyz->boxyz", input, weights)

    def forward(self, x):
        batchsize = x.shape[0]
        #Compute Fourier coeffcients up to factor of e^(- something constant)
        x_ft = torch.fft.rfftn(x, dim=[-3,-2,-1])

        # Multiply relevant Fourier modes
        out_ft = torch.zeros(batchsize, self.out_channels, x.size(-3), x.size(-2), x.size(-1)//2 + 1, dtype=torch.cfloat, device=x.device)
        out_ft[:, :, :self.modes1, :self.modes2, :self.modes3] = \
            self.compl_mul3d(x_ft[:, :, :self.modes1, :self.modes2, :self.modes3], self.weights1)
        out_ft[:, :, -self.modes1:, :self.modes2, :self.modes3] = \
            self.compl_mul3d(x_ft[:, :, -self.modes1:, :self.modes2, :self.modes3], self.weights2)
        out_ft[:, :, :self.modes1, -self.modes2:, :self.modes3] = \
            self.compl_mul3d(x_ft[:, :, :self.modes1, -self.modes2:, :self.modes3], self.weights3)
        out_ft[:, :, -self.modes1:, -self.modes2:, :self.modes3] = \
            self.compl_mul3d(x_ft[:, :, -self.modes1:, -self.modes2:, :self.modes3], self.weights4)

        #Return to physical space
        x = torch.fft.irfftn(out_ft, s=(x.size(-3), x.size(-2), x.size(-1)))
        return x
    
class MLP(nn.Module):
    def __init__(self, in_channels, out_channels, mid_channels):
        super(MLP, self).__init__()
        self.mlp1 = nn.Conv3d(in_channels, mid_channels, 1)
        self.mlp2 = nn.Conv3d(mid_channels, out_channels, 1)

    def forward(self, x):
        x = self.mlp1(x)
        x = F.gelu(x)
        x = self.mlp2(x)
        return x
    

class FNO(nn.Module):
    def __init__(self, modes1 = 19, modes2 = 19, modes3 = 19, width=128, padding=6, in_channels=3, out_channels=6):
        super(FNO, self).__init__()

        """
        The overall network. It contains 4 layers of the Fourier layer.
        1. Lift the input to the desire channel dimension by self.fc0 .
        2. 4 layers of the integral operators u' = (W + K)(u).
            W defined by self.w; K defined by self.conv .
        3. Project from the channel space to the output space by self.fc1 and self.fc2 .
        
        input shape: (Batch=5, Channel=3, uux=38, uuy=38, uuz=38)
        output shape: (batchsize, Channel=6, Txx=38, Tyy=38, Tzz=38, Txy=38, Tyz=38, Txz=38)
        """


        self.modes1 = modes1
        self.modes2 = modes2
        self.modes3 = modes3
        self.width = width
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.padding = 6 # pad the domain if input is non-periodic

        self.normalizer = Normalizer(size=[1, in_channels, 1, 1, 1], name="input")

        self.p = nn.Linear(self.in_channels, self.width)

        self.conv0 = SpectralConv3d(self.width, self.width, self.modes1, self.modes2, self.modes3)
        self.conv1 = SpectralConv3d(self.width, self.width, self.modes1, self.modes2, self.modes3)
        self.conv2 = SpectralConv3d(self.width, self.width, self.modes1, self.modes2, self.modes3)
        self.conv3 = SpectralConv3d(self.width, self.width, self.modes1, self.modes2, self.modes3)

        self.mlp0 = MLP(self.width, self.width, self.width)
        self.mlp1 = MLP(self.width, self.width, self.width)
        self.mlp2 = MLP(self.width, self.width, self.width)
        self.mlp3 = MLP(self.width, self.width, self.width)

        self.w0 = nn.Conv3d(self.width, self.width, 1)
        self.w1 = nn.Conv3d(self.width, self.width, 1)
        self.w2 = nn.Conv3d(self.width, self.width, 1)
        self.w3 = nn.Conv3d(self.width, self.width, 1)

        self.q = MLP(self.width, self.out_channels, self.width * 4)

    
    def forward(self, x):
        x = self.normalizer(x, accumulate=True).float()

        x = x.permute(0, 2, 3, 4, 1)

        x = self.p(x)
 
        x = x.permute(0, 4, 1, 2, 3)

        x = F.pad(x, [0,self.padding]) # pad the domain if input is non-periodic

        x1 = self.conv0(x)
        x1 = self.mlp0(x1)
        x2 = self.w0(x)
        x = x1 + x2
        x = F.gelu(x)

        x1 = self.conv1(x)
        x1 = self.mlp1(x1)
        x2 = self.w1(x)
        x = x1 + x2
        x = F.gelu(x)

        x1 = self.conv2(x)
        x1 = self.mlp2(x1)
        x2 = self.w2(x)
        x = x1 + x2
        x = F.gelu(x)

        x1 = self.conv3(x)
        x1 = self.mlp3(x1)
        x2 = self.w3(x)
        x = x1 + x2

        x = x[..., :-self.padding]
        x = self.q(x)
        x = x.permute(0, 1, 2, 3, 4) # pad the domain if input is non-periodic
        return x.double()
