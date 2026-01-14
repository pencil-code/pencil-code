import torch
from torch import nn
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
        device = batched_data.device 
        data_sum = torch.sum(batched_data, dim=[2,3,4], keepdim=True)
        data_sum_squared = torch.sum(batched_data**2, dim=[2,3,4], keepdim=True)
        
        self.acc_sum = self.acc_sum.to(device)
        self.acc_sum_squared = self.acc_sum_squared.to(device)
        self.acc_count = self.acc_count.to(device)
        self.num_acc = self.num_acc.to(device)

        self.acc_sum += data_sum
        self.acc_sum_squared += data_sum_squared
        self.acc_count += torch.tensor(batched_data[0,0,:,:,:].numel(), device=device)
        self.num_acc += torch.tensor(1.0, device=device)
        
    def _mean(self):
        safe_count = torch.maximum(self.acc_count, torch.tensor(1.0).to(self.acc_count.device))
        return self.acc_sum / safe_count
    
    def _std_with_eps(self):
        safe_count = torch.maximum(self.acc_count, torch.tensor(1.0).to(self.acc_count.device))
        std = torch.sqrt(self.acc_sum_squared / safe_count - self._mean()**2)
        return torch.maximum(std, torch.tensor(self.std_epsilon).to(self.acc_count.device))



class CustomLoss(torch.nn.Module):
  def __init__(self):
    super(CustomLoss, self).__init__()
    self.normalizer = Normalizer([1,6,1,1,1], "output")
    self.counter = 0 # Needed to save the normalisation stats for output
    self.mse = nn.MSELoss()

  def forward(self, prediction, label):
    label = label.double()
    rank = int(str(label.device).split(":")[1]) 
    label = self.normalizer(label, True)
    if self.counter % 100 == 0:
        torch.save(
            {"acc_count": self.normalizer.acc_count, "num_acc": self.normalizer.num_acc, 
             "acc_sum": self.normalizer.acc_sum, "acc_sum_squared": self.normalizer.acc_sum_squared}, 
            f"stats_rank{rank}.pt"
        )
    self.counter += 1

    return self.mse(prediction, label) 
