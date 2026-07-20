import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ——————————————————————
# 1) Paths to your data directories
new_dir = '/gpfs/scratch/arpereira/TABmomentum/model_2/new_no_reverse'   # new files (Δt = 0.05)
old_dir = '/gpfs/scratch/arpereira/TABmomentum/model_2/old_no_reverse'  # old files (Δt = 0.5)

# ——————————————————————
# 2) Gather file lists
new_files = sorted(glob.glob(os.path.join(new_dir, 'pos*.dat')))
old_files = sorted(glob.glob(os.path.join(old_dir, 'pos*.dat')))

# ——————————————————————
# 3) Load all data into numpy arrays
#    Each array has shape (n_steps, 3) for columns [t, x1, x2]
new_data = [np.loadtxt(f, skiprows=1) for f in new_files]
old_data = [np.loadtxt(f, skiprows=1) for f in old_files]

# Extract the time axes
times_new = new_data[0][:, 0]   # 0.00, 0.05, 0.10, …
times_old = old_data[0][:, 0]   # 0.00, 0.50, 1.00, …

# ——————————————————————
# 4) Build index map: for each old time, find its index in new_times
#    Assumes exact matching entries exist.
idx_map = [np.where(times_new == t)[0][0] for t in times_old]

# Stack positions into arrays of shape (n_files, n_steps, 2)
new_coords = np.stack([d[:, 1:] for d in new_data], axis=0)  # (n_new_files, n_new_steps, 2)
old_coords = np.stack([d[:, 1:] for d in old_data], axis=0)  # (n_old_files, n_old_steps, 2)

# Number of animation frames equals number of old time‐steps
n_frames = len(times_old)

# ——————————————————————
# 5) Set up the matplotlib figure
fig, ax = plt.subplots(figsize=(6,6))
ax.set_xlabel('x1')
ax.set_ylabel('x2')
ax.set_title('New (○) vs Old (×) Trajectories without Reversing')

# Create two scatter artists, one for each data set
scatter_new = ax.scatter([], [], marker='o', label='new with no reverse')
scatter_old = ax.scatter([], [], marker='x', label='old with no reverse')

# Text handle for displaying the current time
time_text = ax.text(
    0.95, 0.95,
    '',
    transform=ax.transAxes,
    ha='right', va='top'
)

# Auto‐scale axes to fit all points
all_x = np.concatenate([new_coords[...,0].ravel(), old_coords[...,0].ravel()])
all_y = np.concatenate([new_coords[...,1].ravel(), old_coords[...,1].ravel()])
ax.set_xlim(all_x.min(), all_x.max())
ax.set_ylim(all_y.min(), all_y.max())
ax.legend(loc='lower right')

# ——————————————————————
# 6) Animation functions
def init():
    scatter_new.set_offsets(np.empty((0,2)))
    scatter_old.set_offsets(np.empty((0,2)))
    time_text.set_text('')
    return scatter_new, scatter_old, time_text

def update(frame):
    # Current time from old_times
    t_cur = times_old[frame]

    # Old: direct lookup
    pts_old = old_coords[:, frame, :]        # shape (n_old_files, 2)

    # New: subsample at idx_map
    j = idx_map[frame]
    pts_new = new_coords[:, j, :]            # shape (n_new_files, 2)

    # Update scatters and time text
    scatter_old.set_offsets(pts_old)
    scatter_new.set_offsets(pts_new)
    time_text.set_text(f't = {t_cur:.2f}')
    return scatter_new, scatter_old, time_text

# ——————————————————————
# 7) Build and save the animation
ani = FuncAnimation(
    fig,
    update,
    frames=range(n_frames),
    init_func=init,
    interval=100,   # milliseconds between frames
    blit=True
)

# Save as MP4 (requires ffmpeg installed)
ani.save('trajectories_overlay_aligned.mp4', writer='ffmpeg', dpi=200)
print("Saved aligned overlay as 'trajectories_model2_no_reverse.mp4'")
