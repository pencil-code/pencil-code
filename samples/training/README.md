# Machine Learning Training & Inference Samples

This directory contains samples and scripts to compile/run and create deep learning models.

---

## Samples

### 1. **TG (Taylor-Green Vortex)**

### 2. **Conv-Slab**

### 3. **Helical-MHDturb**

---

## Machine Learning Models

### 1. **UNet3D/CNN3D** `unet.py`

### 2. **Fourier Neural Operator (FNO)** `fno.py`

### 3. **Adaptive Fourier Neural Operator (AFNO)** `afno.py`

---

## Workflow: Training & Inference

### Initial Setup
Navigate to your chosen sample and create necessary symbolic links and add correct modules:

**Prerequisites:** Both the `scripts directory` and the `container image` be symbolicly linked in the current working sample subdirectory.

```bash
module purge
module load gcc/11.3.0
module load /appl/opt/nvhpc-hpcx-cuda11/24.11
module load cmake

cd <sample name>/<training/inference>
ln -s ../../scripts .
```


### 1. Compile the sample

**Configuration:**
- `pc_build -f NVCC_MPI+torch`
- `sbatch scripts/run_native.sh`


### 2. Generate Training/Inference Files
This step creates the three essential files needed for TorchFort:

1. **ML model file** (`.pt`) - TorchScript compiled neural network
2. **Custom loss function** (`.pt`) - TorchScript compiled loss function
3. **Configuration file** (`config_torchfort.yaml`) - Training hyperparameters

**Configuration:**
Edit `scripts/create_training_files.sh` and verify:
- `sample_name`: Matches your current sample
- `data_src`: Directory where training files will be generated (must exist)

**Execution:**
```bash
sbatch scripts/create_training_files.sh
```

**Advanced Usage:**

The script internally calls `build_training_files.py` with the following syntax:

```bash
python build_training_files.py <output_directory> [model_type]
```

**Arguments:**
- `output_directory` (required): Path where files will be generated
- `model_type` (optional): One of `AFNO`, `FNO`, `UNET`, or `CNN`

**Examples:**
```bash
# Generate all four model types
python build_training_files.py ./training_files

# Generate only FNO model
python build_training_files.py ./training_files FNO

# Generate only AFNO model
python build_training_files.py ./training_files AFNO
```

**Note:** If `model_type` is omitted, all four models are created. You must then specify which model to use in `config_torchfort.yaml`.


---

### 3. Run Training or Inference

Once runtime compilation succeeds, use the script `run_native.sh` that calls `pc_run` for subsequent runs.

**Execution:**
```bash
sbatch scripts/run_native.sh
```

---

## Model Configuration

### Configuring Hyperparameters

FNO and AFNO models require hyperparameters to be specified before model creation. Edit `models/model_config.yaml` to customize these settings.

### Fourier Neural Operator (FNO) Parameters

Location: `models/model_config.yaml`

```yaml
FNO:
  modes1: 19          # Fourier modes along x-axis
  modes2: 19          # Fourier modes along y-axis
  modes3: 19          # Fourier modes along z-axis
  width: 128          # Embedding dimension
  padding: 6          # Boundary padding
  in_channels: 3      # Input field channels
  out_channels: 6     # Output field channels
```

---

### Adaptive Fourier Neural Operator (AFNO) Parameters

Location: `models/model_config.yaml`

```yaml
AFNO:
  inp_shape: [38, 38, 38]      # 3D grid dimensions [D, H, W]
  in_channels: 3               # Input field channels
  out_channels: 6              # Output field channels
  patch_size: [1, 1, 1]        # Patch dimensions
  embed_dim: 256               # Embedding dimension
  depth: 4                     # Number of AFNO layers
  mlp_ratio: 4.0               # MLP expansion ratio
  drop_rate: 0.0               # Dropout probability
  num_blocks: 16               # Spectral filter blocks
  sparsity_threshold: 0.01     # Soft-shrink threshold
  hard_thresholding_fraction: 1.0  # Mode retention fraction
```
---

## File Structure

```
samples/training/
├── scripts/
│   ├── compile_submit.sh          # Compilation script
│   ├── create_training_files.sh   # Training file generation
│   ├── build_training_files.py    # Python model builder
│   ├── run_submit.sh              # Single run script
│   ├── run_submit_epoch.sh        # Multi-epoch training script
│   └── run_native.sh              # For running torchfort without the container
├── models/
│   ├── model_config.yaml          # Model hyperparameters
│   ├── unet.py                    # UNet/CNN implementation
│   ├── fno.py                     # FNO implementation
│   └── afno.py                    # AFNO implementation
├── TG/
│   ├── training/
│   └── inference/
├── conv-slab/
│   ├── training/
│   └── inference/
└── helical-MHDturb/
    ├── training/
    └── inference/
```
---


