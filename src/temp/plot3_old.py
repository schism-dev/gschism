from pathlib import Path
import re, glob
import netCDF4 as nc
import numpy as np
import matplotlib
# matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from os.path import join
import os
from scipy.spatial import cKDTree as KDTree


def sorted_out2d_files(folder, start_num=1, end_num=None):
    """Get sorted out2d files with optional range limits."""
    files = glob.glob(str(Path(folder) / "out2d_*.nc"))
    if not files:
        raise FileNotFoundError(f"No out2d_*.nc files in {folder}")
    
    # Extract numbers and filter by range
    file_numbers = []
    for f in files:
        match = re.search(r"out2d_(\d+)\.nc$", f)
        if match:
            num = int(match.group(1))
            if start_num <= num and (end_num is None or num <= end_num):
                file_numbers.append((num, f))
    
    # Sort by number and return file paths
    return [f for num, f in sorted(file_numbers, key=lambda x: x[0])]


def load_velocity_files(folder):
    """Load velocity data from separate files."""
    # Find all velocity files
    vel_x_files = sorted(glob.glob(os.path.join(folder, "horizontalVelX_*.nc")))
    vel_y_files = sorted(glob.glob(os.path.join(folder, "horizontalVelY_*.nc")))
    vert_vel_files = sorted(glob.glob(os.path.join(folder, "verticalVelocity_*.nc")))
    z_files = sorted(glob.glob(os.path.join(folder, "zCoordinates_*.nc")))
    
    print(f"Found {len(vel_x_files)} horizontalVelX files")
    if vel_y_files:
        print(f"Found {len(vel_y_files)} horizontalVelY files")
    if vert_vel_files:
        print(f"Found {len(vert_vel_files)} verticalVelocity files")
    if z_files:
        print(f"Found {len(z_files)} zCoordinates files")
    
    # Load velocity data
    all_velocities_x = []
    all_velocities_y = []
    all_vert_velocities = []
    all_z_coordinates = []
    
    # Load horizontal X velocities
    for file_path in vel_x_files:
        print(f"Reading {os.path.basename(file_path)}")
        with nc.Dataset(file_path, 'r') as dataset:
            velocities_x = dataset.variables['horizontalVelX'][:]  # Shape: (time, nodes, layers)
            velocities_x = np.where(velocities_x >= invalid_threshold, np.nan, velocities_x)
            all_velocities_x.append(velocities_x)
    
    # Load horizontal Y velocities
    if vel_y_files:
        for file_path in vel_y_files:
            print(f"Reading {os.path.basename(file_path)}")
            with nc.Dataset(file_path, 'r') as dataset:
                velocities_y = dataset.variables['horizontalVelY'][:]  # Shape: (time, nodes, layers)
                velocities_y = np.where(velocities_y >= invalid_threshold, np.nan, velocities_y)
                all_velocities_y.append(velocities_y)
    else:
        print(f"No horizontalVelY files found in {folder}")
        if all_velocities_x:
            all_velocities_y = [np.zeros_like(all_velocities_x[0])]  # Create dummy Y velocities
    
    # Load vertical velocities
    if vert_vel_files:
        for file_path in vert_vel_files:
            print(f"Reading {os.path.basename(file_path)}")
            with nc.Dataset(file_path, 'r') as dataset:
                vert_velocities = dataset.variables['verticalVelocity'][:]  # Shape: (time, nodes, layers)
                vert_velocities = np.where(vert_velocities >= invalid_threshold, np.nan, vert_velocities)
                all_vert_velocities.append(vert_velocities)
    else:
        print(f"No vertical velocity files found in {folder}")
        if all_velocities_x:
            all_vert_velocities = [np.zeros_like(all_velocities_x[0])]  # Create dummy vertical velocities
    
    # Load z coordinates
    if z_files:
        for file_path in z_files:
            print(f"Reading {os.path.basename(file_path)}")
            with nc.Dataset(file_path, 'r') as dataset:
                z_coords = dataset.variables['zCoordinates'][:]  # Shape: (time, nodes, layers)
                z_coords = np.where(z_coords >= invalid_threshold, np.nan, z_coords)
                all_z_coordinates.append(z_coords)
    else:
        print(f"No zCoordinates files found in {folder}")
        if all_velocities_x:
            all_z_coordinates = [np.zeros_like(all_velocities_x[0])]  # Create dummy z coordinates
    
    # Concatenate all velocity data
    velocities_x = np.concatenate(all_velocities_x, axis=0) if all_velocities_x else None
    velocities_y = np.concatenate(all_velocities_y, axis=0) if all_velocities_y else None
    vert_velocities = np.concatenate(all_vert_velocities, axis=0) if all_vert_velocities else None
    z_coordinates = np.concatenate(all_z_coordinates, axis=0) if all_z_coordinates else None
    
    return velocities_x, velocities_y, vert_velocities, z_coordinates


def load_run_concat(folder, max_files=None):
    """Memory-efficient loading with explicit cleanup."""
    files = sorted_out2d_files(folder, start_num=1, end_num=max_files)
    
    # Pre-allocate arrays if we know the total size
    all_times = []
    all_elevs = []
    all_dry_flags = []
    x = y = None
    depth = None
    
    for path in files:
        try:
            with nc.Dataset(path) as ds:
                # Read data
                tvar = ds.variables["time"]
                t_py = nc.num2date(tvar[:], units=tvar.units)
                t_np = np.array(t_py, dtype="datetime64[ns]")
                
                elev        = np.array(ds.variables["elevation"][:], copy=True)
                dryFlagNode = np.array(ds.variables["dryFlagNode"][:], copy=True)
                elev = np.where((elev >= invalid_threshold) | (dryFlagNode == 1), np.nan, elev)
                
                if x is None:
                    x = np.array(ds.variables["SCHISM_hgrid_node_x"][:], copy=True)
                    y = np.array(ds.variables["SCHISM_hgrid_node_y"][:], copy=True)
                    depth = np.array(ds.variables["depth"][:], copy=True)
                    depth = np.where(depth >= invalid_threshold, np.nan, depth)
                
                # Append to lists
                all_times.append(t_np)
                all_elevs.append(elev)
                all_dry_flags.append(dryFlagNode)
                
        except Exception as e:
            print(f"Error reading {path}: {e}")
            continue
    
    # Concatenate and immediately delete intermediate data
    time_all = np.concatenate(all_times, axis=0)
    elev_all = np.concatenate(all_elevs, axis=0)
    dry_flags_all = np.concatenate(all_dry_flags, axis=0)
    
    # Explicit cleanup
    del all_times, all_elevs, all_dry_flags
    import gc
    gc.collect()  # Force garbage collection
    
    # Load velocity data from separate files
    velocities_x, velocities_y, vert_velocities, z_coordinates = load_velocity_files(folder)
    
    tree = KDTree(np.column_stack((x, y)).astype(np.float64))
    return {
        "x": x, 
        "y": y, 
        "depth": depth, 
        "time": time_all, 
        "elev": elev_all, 
        "dryFlagNode": dry_flags_all,
        "velocities_x": velocities_x,
        "velocities_y": velocities_y,
        "vert_velocities": vert_velocities,
        "z_coordinates": z_coordinates,
        "tree": tree
    }


# Statistics matrices: rows = cases, columns = stations (nodes)
# Allocated on first case once we know num_stations
stats_std_diff = None   # shape (len(case_names), num_stations)
stats_max_amp1 = None   # max |first run| per (case, station)
stats_max_amp2 = None   # max |second run| per (case, station)

invalid_threshold = 1e10
dpi = 450
case_names = ['RUN02a', 'RUN02a_1']
for case_idx, case_name in enumerate(case_names):
    outdir = r'C:\Users\jmgonzalez\src\temp\%s' % case_name

    # ------------------------------- CONFIG -------------------------------
    # Runs to compare (name + fol
der with out2d_*.nc)
    path = r"C:\Users\jmgonzalez\src\schism_test_suit\%s" % case_name
    runs = [
        # {"name": "GPU", "dir": r"/sciclone/home/jgonzalezondin/scr10/RUN01b"},
        # {"name": "CPU", "dir": r"/sciclone/home/yinglong/DISKS/schism10/External_Exchanges/EX98/RUN01b/outputs"},
        # {"name": "CPU High Res", "dir": join(path, 'outputs_small_nudge')},
        {"name": "CPU", "dir": join(path, 'outputs')},
        {"name": "GPU", "dir": join(path, 'gpuoutput')}, 
        # {"name": "GPU High Res", "dir": join(path, 'gpuoutput_hires')},    
        # add more as needed...
    ]




    run_data = {}
    for r in runs:
        print(f"Loading {r['name']} from {r['dir']} ...")
        rd = load_run_concat(r["dir"])
        rd["name"] = r["name"]
        rd["dir"]  = r["dir"]
        run_data[r["name"]] = rd

    # Plot comparisons of the elevation for each run (cpu and cpu) at the stations.

    # Load station data
    stations_file = Path(__file__).parent / "RUN2Astations.dat"
    stations_data = np.loadtxt(stations_file, skiprows=1)  # Skip header
    station_indices = stations_data[:, 0].astype(int)
    station_x = stations_data[:, 1]
    station_y = stations_data[:, 2]
    station_z = stations_data[:, 3]

    # Group stations by (x, y) to handle multiple z-levels per location
    from collections import defaultdict
    station_groups = defaultdict(list)
    for i, (idx, x, y, z) in enumerate(stations_data):
        # Use rounded coordinates as key to group nearby stations
        key = (round(x, 2), round(y, 2))
        station_groups[key].append((idx, x, y, z, i))

    # Find nearest mesh nodes for each unique station location
    station_node_indices = {}
    for (x_key, y_key), group in station_groups.items():
        # Use first station in group for x, y (they should be the same)
        x_station, y_station = group[0][1], group[0][2]
        
        # Find nearest node using KDTree (use first run's tree)
        first_run = list(run_data.values())[0]
        dist, idx = first_run["tree"].query([x_station, y_station], k=1)
        station_node_indices[(x_key, y_key)] = idx
        print(f"Station at ({x_station:.2f}, {y_station:.2f}) -> Node {idx} (distance: {dist:.2f})")

    # Extract elevation time series for each station from each run
    # Create multiple pages with 32 plots per page

    plots_per_page = 4
    rows_per_page = 2
    cols_per_page = 2  # 4x8 = 32 plots per page

    # Get sorted station groups - sort by minimum station index in each group
    sorted_stations = sorted(station_groups.items(), key=lambda item: min(st[0] for st in item[1]))
    num_stations = len(sorted_stations)
    num_pages = (num_stations + plots_per_page - 1) // plots_per_page  # Ceiling division

    # Allocate stats matrices on first case (rows = cases, columns = stations)
    if stats_std_diff is None:
        stats_std_diff = np.full((len(case_names), num_stations), np.nan)
        stats_max_amp1 = np.full((len(case_names), num_stations), np.nan)
        stats_max_amp2 = np.full((len(case_names), num_stations), np.nan)

    print(f"Creating {num_pages} pages with up to {plots_per_page} plots each")

    for page_num in range(num_pages):
        start_idx = page_num * plots_per_page
        end_idx = min(start_idx + plots_per_page, num_stations)
        stations_on_page = sorted_stations[start_idx:end_idx]
        
        # Create figure for this page
        fig, axes = plt.subplots(rows_per_page, cols_per_page, 
                                figsize=(cols_per_page * 3.5, rows_per_page * 2.5))
        axes = axes.flatten()  # Flatten to 1D array for easier indexing
        
        # # Plot each station on this page
        # for plot_idx, ((x_key, y_key), group) in enumerate(stations_on_page):
        #     ax = axes[plot_idx]
        #     node_idx = station_node_indices[(x_key, y_key)]
        #     x_station, y_station = group[0][1], group[0][2]
            
        #     # Plot zCoordinates at level 0 from each run
        #     lineStyles = ['-', '-']
        #     lineColors = ['blue', 'red']
            
        #     # Store data for difference calculation
        #     first_data = None
        #     second_data = None
            
        #     for run_idx, (run_name, run) in enumerate(run_data.items()):
        #         time_data = run["time"]
        #         # Use zCoordinates at vertical level 0 instead of elevation
        #         z_coords_data = run["z_coordinates"][:, node_idx, 1]  # Level 0 at this node over time
                
        #         linestyle = lineStyles[run_idx]
        #         color = lineColors[run_idx]  # None uses default color cycle
        #         ax.plot(time_data, z_coords_data, label=f"{run_name}", 
        #             linestyle=linestyle, color=color, alpha=0.7, linewidth=0.5)
                
        #         # Store data for difference (first and second runs)
        #         if run_idx == 0:
        #             first_data = z_coords_data
        #         elif run_idx == 1:
        #             second_data = z_coords_data
            
        #     # Plot difference scaled by 10 (convert to dm)
        #     if first_data is not None and second_data is not None:
        #         diff_dm = 10.0 * (first_data - second_data)  # Difference in dm
        #         ax.plot(time_data, diff_dm, label="Difference (10x, dm)", 
        #             linestyle='--', color='green', alpha=0.6, linewidth=0.5)
            
        #     # Compute max absolute value from only first_data and set ylim
        #     if first_data is not None:
        #         line_data = first_data.flatten()
        #         try:
        #             max_abs = np.nanmax(np.abs(line_data))
        #         except Exception:
        #             max_abs = None
        #         if max_abs is not None and not np.isnan(max_abs) and max_abs > 0:
        #             ylim = np.ceil(max_abs * 10) / 10.0
        #             ax.set_ylim(-ylim, ylim)
        #     ax.set_xlabel("Time", fontsize=4)
        #     ax.set_ylabel("Elevation (m)", fontsize=4)
        #     ax.set_title(f"{case_name} - Stn {group[0][0]:.0f}\n({x_station:.0f}, {y_station:.0f})", fontsize=4.5)
        #     ax.legend(fontsize=3.5)
        #     ax.grid(True, alpha=0.3)
        #     ax.tick_params(axis='both', labelsize=3.5)

        #     # Std. dev. of difference and max amplitude of first and second data (on the side of each plot)
        #     if first_data is not None and second_data is not None:
        #         run_names = list(run_data.keys())
        #         name1 = run_names[0] if len(run_names) > 0 else "Run 1"
        #         name2 = run_names[1] if len(run_names) > 1 else "Run 2"
        #         std_diff = np.nanstd(first_data - second_data)
        #         max_amp1 = np.nanmax(np.abs(first_data))
        #         max_amp2 = np.nanmax(np.abs(second_data))
        #         # Store in matrices: rows = cases, columns = stations
        #         station_global_idx = page_num * plots_per_page + plot_idx
        #         stats_std_diff[case_idx, station_global_idx] = std_diff
        #         stats_max_amp1[case_idx, station_global_idx] = max_amp1
        #         stats_max_amp2[case_idx, station_global_idx] = max_amp2

        #         stats_text = (
        #             f"σ(diff) = {std_diff:.4f} m\n"
        #             f"max|{name1}| = {max_amp1:.4f} m\n"
        #             f"max|{name2}| = {max_amp2:.4f} m"
        #         )
        #         ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=3,
        #                 verticalalignment='top', horizontalalignment='right',
        #                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        # # Hide unused subplots
        # for plot_idx in range(len(stations_on_page), plots_per_page):
        #     axes[plot_idx].axis('off')


        
        plt.suptitle(f"Station Elevation Comparison - Page {page_num + 1} of {num_pages}", 
                    fontsize=7, y=0.995)
        plt.tight_layout(rect=[0, 0, 1, 0.99])  # Leave room for suptitle
        
        # Save this page
        run_names = [r['name'].lower().replace(' ', '_') for r in runs]
        filename_suffix = '_vs_'.join(run_names)
        filename = f"station_elevation_comparison_page_{page_num + 1:02d}_{filename_suffix}.png"
        # filename = f"station_elevation_comparison_page_{page_num + 1:02d}_gpu_vs_cpu_high_res.png"
        filename = join(outdir, filename)
        # INSERT_YOUR_CODE
        # Make sure output directory exists before saving
        os.makedirs(outdir, exist_ok=True)

        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        print(f"Saved {filename}")
        plt.close()  # Close figure to free memory

# Save stats matrices (rows = cases, columns = stations)
if stats_std_diff is not None:
    outdir_root = Path(r'C:\Users\jmgonzalez\src\temp')
    outdir_root.mkdir(parents=True, exist_ok=True)
    np.save(outdir_root / 'stats_std_diff.npy', stats_std_diff)
    np.save(outdir_root / 'stats_max_amp1.npy', stats_max_amp1)
    np.save(outdir_root / 'stats_max_amp2.npy', stats_max_amp2)
    np.savetxt(outdir_root / 'stats_std_diff.csv', stats_std_diff, delimiter=',',
               header=','.join(case_names), comments='')
    np.savetxt(outdir_root / 'stats_max_amp1.csv', stats_max_amp1, delimiter=',',
               header=','.join(case_names), comments='')
    np.savetxt(outdir_root / 'stats_max_amp2.csv', stats_max_amp2, delimiter=',',
               header=','.join(case_names), comments='')
    print(f"Saved stats matrices: shape (cases, stations) = {stats_std_diff.shape}")

    # Plot stats: one set of pages, 4 subplots per page; each subplot has 3 lines (σ(diff), max|run1|, max|run2|) vs runs
    num_cases, num_stations = stats_std_diff.shape
    from collections import defaultdict
    stations_data_plot = np.loadtxt(Path(__file__).parent / "RUN2Astations.dat", skiprows=1)
    station_groups_plot = defaultdict(list)
    for i, (idx, x, y, z) in enumerate(stations_data_plot):
        key = (round(x, 2), round(y, 2))
        station_groups_plot[key].append((idx, x, y, z, i))
    sorted_stations_plot = sorted(station_groups_plot.items(), key=lambda item: min(st[0] for st in item[1]))
    plots_per_page_s = 4
    rows_per_page_s = 2
    cols_per_page_s = 2
    num_pages_s = (num_stations + plots_per_page_s - 1) // plots_per_page_s

    for page_num in range(num_pages_s):
        start_idx = page_num * plots_per_page_s
        fig, axes = plt.subplots(rows_per_page_s, cols_per_page_s, figsize=(cols_per_page_s * 3.5, rows_per_page_s * 2.5))
        axes = axes.flatten()
        for plot_idx in range(plots_per_page_s):
            ax = axes[plot_idx]
            station_global_idx = start_idx + plot_idx
            if station_global_idx >= num_stations:
                ax.axis('off')
                continue
            (x_key, y_key), group = sorted_stations_plot[station_global_idx]
            x_station, y_station = group[0][1], group[0][2]
            stn_id = group[0][0]
            x = np.arange(len(case_names))
            # ax.plot(x, stats_std_diff[:, station_global_idx], 'o-', label='σ(diff)', color='#1f77b4', markersize=4)
            ax.plot(x, stats_max_amp1[:, station_global_idx], 's-', label='CPU', color='#ff7f0e', markersize=4)
            ax.plot(x, stats_max_amp2[:, station_global_idx], '^-', label='GPU', color='#2ca02c', markersize=4)
            ax.set_xticks(x)
            ax.set_xticklabels(case_names, fontsize=6, rotation=15)
            ax.set_ylabel('(m)', fontsize=6)
            ax.set_title(f"Stn {stn_id:.0f} ({x_station:.0f}, {y_station:.0f})", fontsize=6)
            ax.legend(fontsize=5)
            ax.grid(True, alpha=0.3)
        plt.suptitle(f"Stats by run - Page {page_num + 1} of {num_pages_s}", fontsize=8, y=0.995)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        fname = outdir_root / f"stats_by_run_page_{page_num + 1:02d}.png"
        plt.savefig(fname, dpi=dpi, bbox_inches='tight')
        plt.close()
        print(f"Saved {fname}")


print(f"All {num_pages} pages saved successfully!")

