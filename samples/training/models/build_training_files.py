import torch
import sys
import os
import yaml

from afno import AFNO
from fno import FNO
from unet import UNet3D
from unet import CNN3D
from normalizer_loss import *

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def build_model(model_save_dir):

    do_all=False
    try:
        with open(f"{os.environ['PENCIL_HOME']}/samples/training/models/model_config.yaml", "r") as f:
            config = yaml.safe_load(f)
    except Exception as e:
        print(f"Failed to load model configuration file : {e}")
        sys.exit(1)
        
    
    
    if len(sys.argv) > 2:
        model_name = sys.argv[2]
    else:
        do_all = True 
    
    models = {}
    # Initialize model
    if(do_all==True or  model_name == 'AFNO' ):
    
        afno_params = config['afno']
    
        model_afno = AFNO(**afno_params).to(device)
        
        models['afno'] = model_afno
    
    if(do_all==True or model_name == 'FNO' ):
        
        fno_params = config['fno']
        model_fno = FNO(**fno_params).to(device)
    
        models['fno'] = model_fno
    
    
    
    if(do_all==True or model_name == 'UNET' ):
    
        model_unet = UNet3D().to(device)
    
        models['unet'] = model_unet
        
    
    
    if(do_all==True or model_name == 'CNN' ):
    
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
            model_scripted = torch.jit.script(models[i])
            save_path = os.path.join(model_save_dir, f"{i}.pt")
            model_scripted.save(save_path)
            print(f"Successfully saved {i}")
        except Exception as e:
            print(f"Failed to script model {i}: {e}")



def build_loss(loss_save_dir):
    loss = CustomLoss()

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
    if model_name is None:
        print("Multiple models created. A specific model name must be provided to the generated a torchfort config file.")
        
    model_name="<model_name>"
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
                "grad_accumulation_steps": 5
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


if len(sys.argv) < 2:
    print("Usage: create_model.py <save_dir> [model_name]")
    sys.exit(1)
    
save_dir=sys.argv[1]

build_model(save_dir)
build_loss(save_dir)
selected_model = sys.argv[2] if len(sys.argv) > 2 else None
build_torchfort_config(save_dir, selected_model)

