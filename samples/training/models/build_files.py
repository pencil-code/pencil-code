import torch
import sys
import os
import yaml
import random
import re
import glob
import struct
import numpy as np

from afno import AFNO
from fno import FNO
from unet import UNet3D
from unet import CNN3D
from normalizer_loss import *

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")



def build_model(model_save_dir, stats_dir, model_name):

    try:
        with open(f"{os.environ['PENCIL_HOME']}/samples/training/models/model_config.yaml", "r") as f:
            config = yaml.safe_load(f)
    except Exception as e:
        print(f"Failed to load model configuration file : {e}")
        sys.exit(1)
        
    
    models = {}
    # Initialize model
    if(model_name == 'afno' ):
    
        afno_params = config['afno']
    
        model_afno = AFNO(**afno_params).to(device)
        
        models['afno'] = model_afno
    
    if(model_name == 'fno' ):
        
        fno_params = config['fno']
        model_fno = FNO(**fno_params).to(device)
    
        models['fno'] = model_fno
    
    
    
    if(model_name == 'unet' ):
    
        model_unet = UNet3D().to(device)
    
        models['unet'] = model_unet
        
    
    
    if(model_name == 'cnn' ):
    
        model_cnn = CNN3D().to(device)
    
        models['cnn'] = model_cnn
    
    
    for i in models:
        # Example input (Batch=1, Channels=3, Depth=70, Height=70, Width=70)
        x = torch.randn(1, 3, 70, 70, 70).to(device)
    
        # Forward pass
        out = models[i](x)
        print(f"Model {i} output shape: {out.shape}")
        
        # Saving as torchscript
        try:
            if torch.os.path.exists(os.path.join(stats_dir, "stats_current_input.pt")):
                models[i].normalizer.load_stats(os.path.join(stats_dir, "stats_current_input.pt"))
                print("creating model with current norm values")
            else:
                models[i].normalizer.load_stats(os.path.join(stats_dir, "stats_restart_input.pt")) 
                print("creating model with new norm values")
            model_scripted = torch.jit.script(models[i])
            save_path = os.path.join(model_save_dir, f"{i}.pt")
            model_scripted.save(save_path)
            print(f"Successfully saved {i}")
        except Exception as e:
            print(f"Failed to script model {i}: {e}")



def build_loss(loss_save_dir, stats_dir):
    loss = CustomLoss()
    if torch.os.path.exists(os.path.join(stats_dir, "stats_current_output.pt")):
        loss.normalizer.load_stats(os.path.join(stats_dir, "stats_current_output.pt"))
        print("creating loss with current norm values")
    else:
        loss.normalizer.load_stats(os.path.join(stats_dir, "stats_restart_output.pt"))
        print("creating loss with new norm values")

    try:
        # Move model to GPU, JIT, and save
        loss.to(device)
        print(f"Loss module moved to {device}")
    except Exception as e:
        print(f"Error moving to device: {e}")
        sys.exit(1)

    try:
        
        loss_jit = torch.jit.script(loss)
        save_path = os.path.join(loss_save_dir, "loss_torchscript_new.pt")
        loss_jit.save(save_path)
        print(f"Successfully saved loss TorchScript to {save_path}")
    except Exception as e:
        print(f"Failed to script loss to dir {save_path}: {e}")
        sys.exit(1)


def build_torchfort_config(save_dir, model_name=None):

    grad_acc = random.randint(5, 10)
    config_data = {
        "model": {
            "type": "torchscript",
            "parameters": {
                "filename":  f"data/training/{model_name.lower()}.pt"
            }
        },
        "loss": {
            "type": "torchscript",
            "parameters": {
                "filename":  "data/training/loss_torchscript_new.pt"
            }
        },
        "optimizer": {
            "type": "adam",
            "parameters": {
                "learning_rate": 1e-3,
                "beta1": 0.9,
                "beta2": 0.999,
                "weight_decay": 0,
                "eps": 1e-8,
                "amsgrad": 0
            },
            "general": {
                "grad_accumulation_steps": grad_acc
            }
        },
        "lr_scheduler": {
            "type": "cosine_annealing",
            "parameters": {
                "T_max": 100000
            }
        }
    }

    config_path = os.path.join(save_dir, "config_torchfort.yaml")
    try:
        with open(config_path, "w") as f:
            yaml.dump(config_data, f, default_flow_style=False, sort_keys=False)
        print(f"Successfully created TorchFort config at: {config_path}")
    except Exception as e:
        print(f"Failed to create TorchFort config: {e}")

def rand_dt_train(run_dir, min, max):
    new_dt = f"{random.uniform(min, max):.10e}"

    lines = []
    with open(run_dir, 'r') as f:
        for line in f:
            if 'dt_train' in line and '=' in line:
                print(line)
                line = re.sub(r'=\s*[\d.eE+-]+', f'= {new_dt}', line)
            lines.append(line)
    
    with open(run_dir, 'w') as f:
        f.writelines(lines)



def avg_stats(stats_dir):
    #sync all norm stats for input
    fpattern = os.path.join(stats_dir, "stats_current_input_rank_*.pt")
        
    files = glob.glob(fpattern)    

    synchronized_data = torch.load(files[0], weights_only=False, map_location=device)

    for f in files[1:]:
        data = torch.load(f, weights_only=False, map_location=device)
        synchronized_data["acc_count"] += data["acc_count"]
        synchronized_data["num_acc"] += data["num_acc"]
        synchronized_data["acc_sum"] += data["acc_sum"]
        synchronized_data["acc_sum_squared"] += data["acc_sum_squared"]

    output = os.path.join(stats_dir,"stats_current_input.pt" )

    torch.save(synchronized_data, output)

 

    #sync all de-norm stats for output
    fpattern = os.path.join(stats_dir, "stats_current_output_rank_*.pt")
        
    files = glob.glob(fpattern)    

    synchronized_data = torch.load(files[0], weights_only=False, map_location=device)

    for f in files[1:]:
        data = torch.load(f, weights_only=False, map_location=device)
        synchronized_data["acc_count"] += data["acc_count"]
        synchronized_data["num_acc"] += data["num_acc"]
        synchronized_data["acc_sum"] += data["acc_sum"]
        synchronized_data["acc_sum_squared"] += data["acc_sum_squared"]

    output = os.path.join(stats_dir,"stats_current_output.pt" )

    torch.save(synchronized_data, output)

def ptTObin(stats_dir):
    stats_file_open = os.path.join(stats_dir, "stats_current_output.pt")
    stats = torch.load(stats_file_open, weights_only=False, map_location=device)

    arrays = {
        "acc_count": stats["acc_count"].cpu().numpy(),
        "num_acc": stats["num_acc"].cpu().numpy(),
        "acc_sum": stats["acc_sum"].cpu().numpy(),
        "acc_sum_squared": stats["acc_sum_squared"].cpu().numpy()
    }

    acc_count = arrays["acc_count"][0]  # 652,044,350.0
    acc_sum = arrays["acc_sum"]        # Shape (6, 1, 1, 1, 1)
    acc_sum_sq = arrays["acc_sum_squared"] 
    epsilon = 1e-8

    # 2. Calculate Mean (acc_sum / acc_count)
    # This matches your _mean() method
    means = acc_sum / acc_count

    # 3. Calculate Std (sqrt(acc_sum_sq / acc_count - mean^2))
    # This matches your _std_with_eps() method
    variance = (acc_sum_sq / acc_count) - (means**2)
    stds = np.sqrt(np.maximum(variance, 0)) # Clip to 0 for numerical stability
    stds = np.maximum(stds, epsilon)        # Apply epsilon

    # 4. Display Results (Flattened for readability)
    component_means = means.flatten()
    component_stds = stds.flatten()

    for i in range(len(component_means)):
        print(f"Component {i}: Mean = {component_means[i]:.8f}, Std = {component_stds[i]:.8f}")



    stats_file_save = os.path.join(stats_dir, "normalizer.bin")
    with open(stats_file_save, "wb") as f:
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
