# exec(open('/home/simon/pencil-code/simon/hopf_kink/plot_streamlines_blender.py').read())

import numpy as np
import blendaviz as blt
import sys
import os

# Add Pencil Code directory to the path.
sys.path.append('/home/simon/pencil-code/python')
import pencil as pc

# Define the root directory ewhere we store the data.
#root_dir = '/home/iomsn/pencil-code/simon/hopf_kink'
# root_dir = '/media/iomsn/2f9ec593-2396-4a79-9e65-d245c63f67c2/pencil-code/simon/hopf_kink/dag'
root_dir = '/media/simon/2f9ec593-2396-4a79-9e65-d245c63f67c2/pencil-code/simon/hopf_kink/dag'

# Define the run directory with the data to be plotted.
#run_dir = 'test_b_ext'
run_dir = 'n_256_eta_0_Lxyz_4_mehta'

# Define the time index.
t_idx = 50

# Define the streamline parameters.
n_seeds = 300
integration_time = 100
seed_radius = 1.1


def generate_seeds(n_seeds, seed_radius, random_seed=None):
    """
    Generate n_seeds points uniformly distributed inside a sphere of radius seed_radius.

    Parameters:
        n_seeds (int): Number of points to generate
        seed_radius (float): Radius of the sphere
        random_seed (int, optional): Seed for reproducibility

    Returns:
        points (np.ndarray): Array of shape (n_seeds, 3) with 3D coordinates
    """
    if random_seed is not None:
        np.random.seed(random_seed)  # Set the seed for reproducibility

    # Random angles
    phi = np.random.uniform(0, 2 * np.pi, n_seeds)
    costheta = np.random.uniform(-1, 1, n_seeds)
    u = np.random.uniform(0, 2, n_seeds)  # For radius scaling

    theta = np.arccos(costheta)
    r = seed_radius * np.cbrt(u)  # Cube root ensures uniform density

    # Convert spherical to Cartesian coordinates
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    points = np.vstack((x, y, z)).T

    return points


seeds = generate_seeds(n_seeds, seed_radius, 0)

# Read the data.
print('loading {0}'.format(t_idx))
var = pc.read.var(datadir=os.path.join(root_dir, run_dir, 'data'), magic='bb', var_file='VAR{0}'.format(t_idx))

# Reshape the data.
bb = var.bb.swapaxes(1, 3)
b2 = bb[0]**2 + bb[1]**2 + bb[2]**2

# Plot the data.
streamlines = blt.streamlines_array(var.x, var.y, var.z, bb[0], bb[1], bb[2],
                                    integration_time=integration_time,
                                    seeds=seeds,
                                    atol=1e-3, rtol=1e-3, integration_steps=400, radius=0.01,
                                    color_scalar='magnitude')
