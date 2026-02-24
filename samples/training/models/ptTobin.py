import numpy as np
import struct
import torch

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

stats = torch.load("stats_rank0.pt", weights_only=False, map_location=device)


arrays = {
    "acc_count": stats["acc_count"].cpu().numpy(),
    "num_acc": stats["num_acc"].cpu().numpy(),
    "acc_sum": stats["acc_sum"].cpu().numpy(),
    "acc_sum_squared": stats["acc_sum_squared"].cpu().numpy()
}

with open("normalizer.bin", "wb") as f:
    for name, arr in arrays.items():
        # write name length and name (for debugging or identification)
        name_bytes = name.encode("utf-8")
        f.write(struct.pack("I", len(name_bytes)))
        f.write(name_bytes)
        
        # write shape info
        f.write(struct.pack("I", arr.ndim))
        f.write(struct.pack("I" * arr.ndim, *arr.shape))
        
        # write the raw float data
        f.write(arr.astype(np.float32).tobytes())
