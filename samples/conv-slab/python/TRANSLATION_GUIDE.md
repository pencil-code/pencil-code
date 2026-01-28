# IDL to Python Translation Guide for Pencil Code Analysis

## Overview
This guide explains how to translate IDL analysis scripts (like `thermo.pro` and `pc_fluxz.pro`) to Python using the **Pencil Code Python library** (`pencil-code/python/pencil/`).

---

## Key Python Functions from Pencil Code

### 1. **Reading Parameters**
Two approaches are available:

**Option A: Using `pc.read.param()` (direct, more flexible)**
```python
import pencil as pc

par = pc.read.param(datadir='data', quiet=False)
cp = par.cp           # Access parameter as object attribute
gamma = par.gamma
cs20 = getattr(par, 'cs20', 1.0)  # with default for missing params
```

**Option B: Using `sim.param` (recommended, faster, cached)**
```python
sim = pc.sim.simulation('path/to/sim')
par = sim.param       # Returns a dictionary (cached)
cp = par.get('cp', 1.0)        # Use .get() with default value
gamma = par.get('gamma', 1.6666666)
cs20 = par.get('cs0', 1.0)**2  # Handle parameter name variations
```

**Comparison Table:**

| Aspect | `pc.read.param()` | `sim.param` |
|--------|-------------------|-----------|
| **Type** | Object with attributes | Dictionary |
| **Access** | `par.cp` or `getattr(par, 'cp', default)` | `par.get('cp', default)` |
| **Speed** | Slower (file I/O every call) | Faster (cached in simulation object) |
| **Default handling** | Need `getattr()` | Built-in `.get()` method |
| **Best use** | When you need fresh data | Most analysis scripts |

**Equivalent IDL:**
```idl
pc_read_param, obj=par
cp = par.cp
```

### 2. **Reading Grid Dimensions**

**Option A: Using `pc.read.dim()` (direct)**
```python
dim = pc.read.dim(datadir='data')
print(dim.nx, dim.ny, dim.nz)
z = dim.z  # z-coordinate array
```

**Option B: Using `sim.dim` (recommended, cached)**
```python
sim = pc.sim.simulation('path/to/sim')
dim = sim.dim          # Returns an object (cached)
nx = getattr(dim, 'nx', None)
z = dim.z              # z-coordinate array
nz = len(dim.z)
```

Both return the same object type, but `sim.dim` is cached by the simulation object and faster for repeated access.

### 3. **Reading Time Series** (`pencil.read.ts`)
```python
ts = pc.read.ts(sim=sim)
print(ts.t)      # time array
print(ts.urms)   # RMS velocity
print(ts.umax)   # max velocity
```

**Equivalent IDL:**
```idl
pc_read_ts, obj=ts
```

### 4. **Reading Snapshot Data** (`pencil.read.var`)
```python
var = pc.read.var(sim=sim, trimall=True)
lnrho = var.lnrho     # log density (3D array: nx x ny x nz)
ss = var.ss           # specific entropy (3D array)
uu = var.uu           # velocity (4D array: nx x ny x nz x 3)
```

**Equivalent IDL:**
```idl
pc_read_var, obj=var, /trimall
```

### 5. **Reading Xy-Averaged Data** (`pencil.read.aver`)
```python
aver = pc.read.aver(sim=sim, quiet=False)
t = aver.t           # time array
z = aver.z           # z-coordinate array
fradz = aver.fradz   # radiative flux (nz x nt)
fkinz = aver.fkinz   # kinetic energy flux (nz x nt)
fconvz = aver.fconvz # convective flux (nz x nt)
fcoolz = aver.dcoolz # cooling rate (nz x nt)
```

This reads from `xyaverages.dat` which contains horizontally-averaged (xy-averaged) quantities at each (z, t) position. Shape is typically `(nz, nt)` for 1D profiles over time.

**Equivalent IDL:**
```idl
pc_read_xyaver, obj=aver
; or access via rxyaver.pro which loads into common block variables
```

### 6. **Reading Slices** (`pencil.read.slices`)
```python
slices = pc.read.slices(sim=sim)
uuz_xy = slices.xy.uu3        # vertical velocity at top (xy plane, shape: nt x nx x ny)
rho_xz = np.exp(slices.xz.lnrho)  # density in vertical plane (nt x nx x nz)
t_slice = slices.xy.t         # time array for slices
```

**Equivalent IDL:**
```idl
pc_read_slices, obj=slices
```

---

## Common IDL → Python Translation Patterns

### **Parameter Access**
```idl
; IDL - with default value
cp = par.cp
if n_elements(cp) eq 0 then cp = 1.0  ; default value
```

```python
# Python - Method 1: Using object attribute
cp = getattr(par, 'cp', 1.0)  # returns 1.0 if cp doesn't exist

# Python - Method 2: Using dictionary (sim.param)
cp = par.get('cp', 1.0)  # cleaner with dictionary
```

### **Array Operations**
```idl
; IDL - conditional indexing
top = where(z ge z2)
Tinit[top] = expression
```

```python
# Python - boolean masking
top_mask = z >= z2
Tinit[top_mask] = expression
```

### **Array Initialization**
```idl
; IDL
ssinit = (lnrhoinit = (Tinit = 0.*z))
```

```python
# Python
Tinit = np.zeros_like(z)
lnrhoinit = np.zeros_like(z)
ssinit = np.zeros_like(z)
```

### **Logarithms and Exponentials**
```idl
; IDL
alog(x)    ; natural log
exp(x)     ; exponential
```

```python
# Python
np.log(x)   # natural logarithm
np.exp(x)   # exponential
```

### **Maximum/Minimum Operations**
```idl
; IDL
value > threshold   ; minimum (clamp low)
value < threshold   ; maximum (clamp high)
```

```python
# Python
np.maximum(value, threshold)  # max
np.minimum(value, threshold)  # min
```

### **Vertical Integration**
```idl
; IDL - using numerical integration
result = integr(data, x=z)
```

```python
# Python - using NumPy trapezoid rule
result = np.trapz(data, z)  # scalar result
# For 2D (nz, nt) integrating over z axis:
result = np.trapz(data_2d, z, axis=0)  # returns (nt,) array
```

### **Plotting Function Structure**
Both scripts organize plots inside functions:

```python
# Python - put plotting inside function
def plot_temperature_profile(z, T, title="Temperature"):
    """Plot temperature vs height."""
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(z, T, 'b-', linewidth=2)
    ax.set_xlabel('Height (z)')
    ax.set_ylabel('Temperature')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    return fig, ax

# Call the function
fig, ax = plot_temperature_profile(z, T_init)
fig.savefig('temperature.png', dpi=150)
plt.show()
```

---

## Translation of `thermo.pro` to Python

### **Key Equations**
From the ideal gas equation of state:

```
pp  = cs20 * rho0 / gamma * exp(gamma*(ss/cp + lnrho - lnrho0))
cs2 = cs20 * exp(gamma*ss/cp + gamma_m1*(lnrho - lnrho0))
TT  = cs2 / gamma_m1 / cp
```

### **Python Implementation Steps**

1. **Read parameters from simulation**
   ```python
   import pencil as pc
   sim = pc.sim.simulation('path/to/sim')
   par = sim.param  # Recommended: use sim.param (dictionary)
   dim = sim.dim
   ```

2. **Extract thermodynamic constants**
   ```python
   cp = par.get('cp', 1.0)
   rho0 = par.get('rho0', 1.0)
   cs20 = par.get('cs0', 1.0)**2
   gamma = par.get('gamma', 1.6666666)
   gamma_m1 = gamma - 1.0
   ```

3. **Read snapshot data** (optional)
   ```python
   var = pc.read.var(sim=sim)
   lnrho = var.lnrho
   ss = var.ss
   ```

4. **Compute thermodynamic variables**
   ```python
   lnrho0 = np.log(rho0)
   pp = cs20 * rho0 / gamma * np.exp(gamma * (ss/cp + lnrho - lnrho0))
   cs2 = cs20 * np.exp(gamma*ss/cp + gamma_m1*(lnrho - lnrho0))
   TT = cs2 / gamma_m1 / cp
   ```

---

## Translation of `pc_fluxz.pro` to Python

### **Key Steps**

1. **Read xy-averaged data containing fluxes**
   ```python
   aver = pc.read.aver(sim=sim)
   t = aver.t
   z = aver.z
   
   # Extract available flux components
   fradz = aver.fradz if hasattr(aver, 'fradz') else None
   fkinz = aver.fkinz if hasattr(aver, 'fkinz') else None
   fconvz = aver.fconvz if hasattr(aver, 'fconvz') else None
   dcoolz = aver.dcoolz if hasattr(aver, 'dcoolz') else None
   ```

2. **Compute cooling flux by integration**
   ```python
   # cooling flux: fcool = -∫ dcool dz
   fcool = np.zeros_like(dcoolz)
   for it in range(1, len(t)):
       fcool[:, it] = -np.trapz(dcoolz[:, it], z)
   ```

3. **Compute total flux**
   ```python
   ftot = fradz + fkinz + fconvz + fcool
   ```

4. **Create visualization functions** (see pc_flux.py for examples)

---

## Useful Pencil Python Modules

### **Module: `pencil.read`**
- `param()` - Read parameters (start.in, run.in)
- `ts()` - Read time series (time_series.dat)
- `var()` - Read 3D snapshots (var files)
- `slices()` - Read 2D slices (video files)
- `dim()` - Read grid dimensions
- `grid()` - Read grid coordinates (x, y, z)
- `aver()` - Read xy-averaged data (xyaverages.dat)

### **Module: `pencil.calc`**
- `curl()` - Compute curl of a vector field
- `Reynolds()` - Compute Reynolds number components
- `gravity.g_eff()` - Effective gravity in stratified atmosphere

### **Module: `pencil.visu`**
- `animate_interactive()` - Interactive 2D animation of time-evolving data
- `contour_plot()` - Enhanced contour plotting

### **Module: `pencil.sim`**
- `simulation()` - Create simulation object (manages paths and metadata, caches dim and param)

---

## Debugging Tips

### **Check available parameters**
```python
# Method 1: List all dictionary keys
par = sim.param
print(par.keys())

# Method 2: Using object attributes
par = pc.read.param(datadir=sim.datadir, quiet=False)
print([attr for attr in dir(par) if not attr.startswith('_')])
```

### **Check available averaged quantities**
```python
aver = pc.read.aver(sim=sim, quiet=False)
print([attr for attr in dir(aver) if not attr.startswith('_')])
```

### **Verify data shapes**
```python
print(f"lnrho shape: {var.lnrho.shape}")  # should be (nx, ny, nz)
print(f"uu shape: {var.uu.shape}")        # should be (nx, ny, nz, 3)
print(f"aver frad shape: {aver.fradz.shape}")  # should be (nz, nt)
```

### **Handle missing attributes gracefully**
```python
# With sim.param (dictionary)
cs20 = par.get('cs0', 1.0)**2

# With pc.read.param (object)
cs20 = getattr(par, 'cs0', 1.0)**2
```

---

## Common Errors and Solutions

| Error | Cause | Solution |
|-------|-------|----------|
| `FileNotFoundError` | Wrong `datadir` path | Use `sim.datadir` from simulation object |
| `KeyError: 'cp'` | Parameter not in sim.param dict | Use `par.get('cp', 1.0)` with default |
| `AttributeError: 'Aver' has no attribute 'fradz'` | Flux not in simulation | Check with `hasattr(aver, 'fradz')` first |
| `ValueError: operands could not be broadcast` | Array shape mismatch | Check dimensions with `.shape` |
| `ModuleNotFoundError: pencil` | Pencil not installed | Run: `pip install pencil-code` |

---

## Running Your Scripts

```bash
# From command line
python thermo.py
python pc_flux.py

# From IPython (interactive)
ipython
%run thermo.py
%run pc_flux.py

# Check script help
python -c "import pc_flux; help(pc_flux.plot_fluxes_final_time)"
```

---

## File Organization Best Practices

```
conv-slab/
├── python/
│   ├── thermo.py                  # Thermodynamic analysis
│   ├── pc_flux.py                 # Energy flux visualization
│   ├── TRANSLATION_GUIDE.md       # This file
│   └── README.md                  # Local Python analysis guide
├── idl/
│   ├── thermo.pro                 # Original IDL script
│   └── pc_fluxz.pro               # Original IDL script
└── data/
    ├── xyaverages.dat             # Horizontal averages
    ├── time_series.dat            # Time series
    └── var files...
```

---

## References

- **Pencil Code Python**: `pencil-code/python/README`
- **Example Scripts**: `pencil-code/samples/conv-slab/python/`
- **Pencil Code Manual**: http://pencil-code.googlecode.com/
- **NumPy Documentation**: https://numpy.org/doc/
- **Matplotlib Guide**: https://matplotlib.org/stable/index.html
