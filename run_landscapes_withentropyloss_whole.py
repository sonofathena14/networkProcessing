import gudhi
import numpy as np
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt
from gudhi.representations import Landscape
import gudhi.representations
from itertools import combinations
import seaborn as sns
from gudhi.hera import wasserstein_distance
from gudhi.hera import bottleneck_distance
from collections import deque, defaultdict
from scipy.interpolate import interp1d
from scipy.stats import ks_2samp
from scipy.spatial.distance import pdist
from landscape_heatmap_functions import generate_all_landscape_heatmaps
def nodesExtractorC(name): #extracts nodes and their corresponding information
    file_path = 'Networks/Network_Vessels_' + name +'.mat'
    matlab_data = scipy.io.loadmat(file_path)
    # Extract the 'connectivity' field from the 'Data' structured array
    data_structure = matlab_data['nodesC2']
    # Reshape or ensure it's a proper 2D array (if required)
    nodes_data = data_structure.squeeze()
    # Create a DataFrame from the connectivity data
    nodes_df = pd.DataFrame(nodes_data, columns=['NodeID', 'X', 'Y', 'Z', 'Degree'])
    # Save the DataFrame to inspect it
    return nodes_df

def nodesExtractorH(name): #extracts nodes and their corresponding information
    file_path = 'Pruned/Pruned_Network_' + name +'.mat'
    matlab_data = scipy.io.loadmat(file_path)
    # Extract the 'connectivity' field from the 'Data' structured array
    data_structure = matlab_data['nodesC3']
    # Reshape or ensure it's a proper 2D array (if required)
    nodes_data = data_structure.squeeze()
    # Create a DataFrame from the connectivity data
    nodes_df = pd.DataFrame(nodes_data, columns=['NodeID', 'X', 'Y', 'Z', 'Degree'])
    # Save the DataFrame to inspect it
    return nodes_df

def nodesToArrayC(name):
    nodes_df = nodesExtractorC(name)
    nodes_loc = nodes_df.loc[:,['X','Y','Z']]/1000
    loc_array = nodes_loc.to_numpy()
    return loc_array

def nodesToArrayH(name):
    nodes_df = nodesExtractorH(name)
    nodes_loc = nodes_df.loc[:,['X','Y','Z']]/1000
    loc_array = nodes_loc.to_numpy()
    return loc_array

def compute_pointcloud_diameter(points):
    # Returns the maximum Euclidean distance (diameter)
    dists = pdist(points)
    return np.max(dists)

def compute_persistence(points):
    """
    Compute persistence diagram and track the last death in each dimension.
    Excludes infinite values in H1 and H2.
    """
    # --- Build alpha complex ---
    alpha_complex = gudhi.AlphaComplex(points=points)
    simplex_tree = alpha_complex.create_simplex_tree()
    simplex_tree.compute_persistence()
    
    # --- Collect persistence pairs ---
    persistence_pairs = simplex_tree.persistence()
    
    # DIAGNOSTIC: Check for infinite bars
    for dim, (birth, death) in persistence_pairs:
        if np.isinf(death):
            print(f"WARNING: Infinite bar in H{dim}: birth={birth}, death=∞")
    
    # --- Filter infinite bars for H1/H2 ---
    diag = []
    for dim, (birth, death) in persistence_pairs:
        diag.append((dim, (birth, death)))
    
    # --- Track last death in each dimension ---
    last_deaths = {}
    for dim, (birth, death) in diag:
        if np.isfinite(death):
            if dim not in last_deaths:
                last_deaths[dim] = death
            else:
                last_deaths[dim] = max(last_deaths[dim], death)
    
    return diag, last_deaths


def average_curves(curve_dicts, resolution=100000):
    """
    Generic averaging function for any type of curve dicts (Betti, lifespan, etc.)
    curve_dicts: list of dicts {dim: (grid, curve)}
    Returns: dict {dim: (common_grid, avg_curve)}
    """
    average = {}

    dims = set().union(*(d.keys() for d in curve_dicts))

    for dim in dims:
        all_grids, all_curves = [], []
        for cdict in curve_dicts:
            if dim not in cdict:
                continue
            grid, curve = cdict[dim]
            if len(grid) == 0:
                continue
            all_grids.append(grid)
            all_curves.append(curve)

        if not all_grids:
            average[dim] = (np.array([]), np.zeros(resolution))
            continue

        global_min = min(g[0] for g in all_grids)
        global_max = max(g[-1] for g in all_grids)
        common_grid = np.linspace(global_min, global_max, resolution)

        interpolated = []
        for grid, curve in zip(all_grids, all_curves):
            f = interp1d(grid, curve, bounds_error=False, fill_value=0)
            interpolated.append(f(common_grid))

        avg_curve = np.mean(interpolated, axis=0)
        average[dim] = (common_grid, avg_curve)

    return average

def plot_curves(curves_dicts, labels, dimension, title='Betti Curves Comparison', limit = [0,100]):
    """
    curves_dicts: list of betti_curve dicts {dim: (grid, curve)}
    labels: list of labels for each curve
    dimension: the Betti number (int) to plot
    """
    plt.figure(figsize=(10, 6))

    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'cyan', 'magenta']
    linestyles = ['-', '--', '-.', ':']

    for idx, (curve_dict, label) in enumerate(zip(curves_dicts, labels)):
        if dimension in curve_dict:
            grid, curve = curve_dict[dimension]
            if len(grid) > 0:
                color = colors[idx % len(colors)]
                linestyle = linestyles[idx % len(linestyles)]
                plt.plot(grid, curve, color=color, linestyle=linestyle, label=f'{label} (Dimension-{dimension})')
        else:
            print(f"Warning: {label} does not contain Betti-{dimension}")

    plt.title(f'{title} (Betti-{dimension})', fontsize=28)
    plt.xlabel('Filtration Value', fontsize=20)
    plt.ylabel(f'Betti-{dimension}', fontsize=20)
    plt.xlim(limit)
    plt.legend(fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(False)
    plt.tight_layout()
    plt.show()

def save_curves(curves_dicts, labels, dimension, filename, title='Betti Curves Comparison', limit=[0, 100], dpi=300):
    """
    curves_dicts: list of betti_curve dicts {dim: (grid, curve)}
    labels: list of labels for each curve
    dimension: the Betti number (int) to plot
    filename: path to save the plot (e.g., 'output.png' or 'output.jpg')
    dpi: resolution of saved image (default 300)
    """
    plt.figure(figsize=(10, 6))
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'cyan', 'magenta']
    linestyles = ['-', '--', '-.', ':']
    
    for idx, (curve_dict, label) in enumerate(zip(curves_dicts, labels)):
        if dimension in curve_dict:
            grid, curve = curve_dict[dimension]
            if len(grid) > 0:
                color = colors[idx % len(colors)]
                linestyle = linestyles[idx % len(linestyles)]
                plt.plot(grid, curve, color=color, linestyle=linestyle, label=f'{label} (Dimension-{dimension})')
        else:
            print(f"Warning: {label} does not contain Betti-{dimension}")
    
    plt.title(f'{title} (Betti-{dimension})', fontsize=28)
    plt.xlabel('Filtration Value', fontsize=20)
    plt.ylabel(f'Betti-{dimension}', fontsize=20)
    plt.xlim(limit)
    plt.legend(fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    plt.close()  # Close the figure to free memory
    
    print(f"Saved plot to {filename}")

def find_global_k_per_dimension(
    persistence_diagrams,
    t_min,
    t_max,
    max_dim=2,
    resolution=500,
    energy_threshold=0.95
):
    """
    Find optimal k for each dimension across ALL datasets using cumulative energy threshold.
    Returns dict {dim: k_value}
    
    Parameters:
    -----------
    persistence_diagrams : dict
        Dictionary of persistence diagrams
    t_min, t_max : float
        Range for landscape computation
    max_dim : int
        Maximum homology dimension to consider
    resolution : int
        Resolution for landscape functions
    energy_threshold : float
        Cumulative energy threshold (default 0.95 for 95%)
    """
    def compute_optimal_k_single(diagram, t_min, t_max, resolution, energy_threshold):
        """Helper to find k for a single diagram using cumulative energy"""
        if diagram.shape[0] == 0:
            return 1
        
        max_k = min(21,diagram.shape[0])
        
        landscape = Landscape(
            num_landscapes=max_k,
            resolution=resolution,
            sample_range=(t_min, t_max)
        )
        landscapes = landscape.fit_transform([diagram])[0]
        landscapes = landscapes.reshape(max_k, resolution)
        
        # Compute L2 norms
        norms = np.linalg.norm(landscapes, axis=1)
        
        if norms.sum() == 0:
            return 1
        
        # Compute cumulative energy
        total_energy = np.sum(norms)
        cumulative_energy = np.cumsum(norms) / total_energy
        
        # Find first k where cumulative energy exceeds threshold
        k_optimal = np.argmax(cumulative_energy >= energy_threshold) + 1
        
        if k_optimal == 1 and cumulative_energy[0] < energy_threshold:
            # Threshold not reached, use all landscapes
            k_optimal = max_k
        
        return max(1, k_optimal)
    
    # Collect k values for each dimension across all datasets
    k_per_dim = {dim: [] for dim in range(max_dim + 1)}
    
    for name, pers_diag in persistence_diagrams.items():
        for dim in range(max_dim + 1):
            diagram = np.array([
                (birth, death)
                for d, (birth, death) in pers_diag
                if d == dim and np.isfinite(death)
            ])
            
            k = compute_optimal_k_single(diagram, t_min, t_max, resolution, energy_threshold)
            k_per_dim[dim].append(k)
    
    # Take maximum k for each dimension across all datasets
    global_k = {dim: max(k_list) if k_list else 1 for dim, k_list in k_per_dim.items()}
    
    return global_k

def compute_landscape_curves(
    persistence_diagram,
    t_min,
    t_max,
    max_dim=2,
    resolution=500,
    k_values=None
):
    """
    Compute persistence landscapes for all dimensions.
    
    Args:
        persistence_diagram: List of (dim, (birth, death)) tuples
        t_min, t_max: Range for landscape functions
        max_dim: Maximum dimension to compute
        resolution: Number of sample points
        k_values: Dict {dim: k} specifying k for each dimension
    
    Returns:
        dict of {dim: (t_vals, landscapes)} where landscapes is shape (k, resolution)
    """
    landscape_curves = {}
    
    for dim in range(max_dim + 1):
        diagram = np.array([
            (birth, death)
            for d, (birth, death) in persistence_diagram
            if d == dim and np.isfinite(death)
        ])
        
        # Determine k for this dimension
        if k_values is not None:
            k = k_values[dim]
        else:
            k = 5  # default fallback
        
        if diagram.shape[0] == 0:
            # No features in this dimension
            t_vals = np.linspace(t_min, t_max, resolution)
            landscapes = np.zeros((k, resolution))
            landscape_curves[dim] = (t_vals, landscapes)
            continue
        
        # Adjust k if diagram has fewer intervals
        k = min(k, diagram.shape[0])
        
        # Compute landscapes with specified k
        landscape = Landscape(
            num_landscapes=k,
            resolution=resolution,
            sample_range=(t_min, t_max)
        )
        landscapes = landscape.fit_transform([diagram])[0]
        landscapes = landscapes.reshape(k, resolution)
        
        t_vals = np.linspace(t_min, t_max, resolution)
        landscape_curves[dim] = (t_vals, landscapes)
    
    return landscape_curves


def save_landscape(t_vals, landscapes, dim, filename, title=None, limit=[0, 100], dpi=300):
    """
    Save persistence landscape plot to file.
    landscapes: shape (k, resolution) array of k landscape functions
    """
    plt.figure(figsize=(10, 6))
    
    for i, layer in enumerate(landscapes):
        plt.plot(t_vals, layer, label=f"$\\lambda_{{{i+1}}}$")
    
    plt.xlabel("Filtration Value", fontsize=20)
    plt.ylabel("Landscape Value", fontsize=20)
    plt.title(title or f"Persistence Landscape (H_{dim})", fontsize=28)
    plt.xlim(limit)
    plt.legend(fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(False)
    plt.tight_layout()
    
    plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"Saved landscape plot to {filename}")

def plot_landscape(t_vals, landscapes, dim, title=None):
    plt.figure(figsize=(8, 5))
    for i, layer in enumerate(landscapes):
        plt.plot(t_vals, layer, label=f"$\\lambda_{{{i+1}}}$")

    plt.xlabel("t")
    plt.ylabel("Landscape value")
    plt.title(title or f"Persistence Landscape (H_{dim})")
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_persistence_diagram(persistence_diagram, title="Persistence Diagram"):
    """
    Plot H0 and H1 persistence diagrams side by side.
    
    Parameters:
    -----------
    persistence_diagram : list or dict
        If list: expects format [(dim, (birth, death)), ...]
        If dict: expects format {0: [(birth, death), ...], 1: [(birth, death), ...]}
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Parse the persistence diagram format
    if isinstance(persistence_diagram, dict):
        h0_points = persistence_diagram.get(0, [])
        h1_points = persistence_diagram.get(1, [])
    else:
        # Assume format [(dim, (birth, death)), ...]
        h0_points = [(b, d) for dim, (b, d) in persistence_diagram if dim == 0]
        h1_points = [(b, d) for dim, (b, d) in persistence_diagram if dim == 1]
    
    # Convert to numpy arrays
    h0_points = np.array(h0_points) if len(h0_points) > 0 else np.array([]).reshape(0, 2)
    h1_points = np.array(h1_points) if len(h1_points) > 0 else np.array([]).reshape(0, 2)
    
    # Plot H0
    ax = axes[0]
    if len(h0_points) > 0:
        # Filter out infinite death times
        finite_h0 = h0_points[np.isfinite(h0_points[:, 1])]
        infinite_h0 = h0_points[~np.isfinite(h0_points[:, 1])]
        
        if len(finite_h0) > 0:
            ax.scatter(finite_h0[:, 0], finite_h0[:, 1], alpha=0.6, s=50, label='Finite')
            max_val = max(finite_h0[:, 0].max(), finite_h0[:, 1].max())
        else:
            max_val = 1
            
        if len(infinite_h0) > 0:
            ax.scatter(infinite_h0[:, 0], [max_val * 1.1] * len(infinite_h0), 
                      marker='^', s=100, c='red', label='Infinite', alpha=0.6)
            max_val = max_val * 1.2
        
        # Plot diagonal
        ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='Diagonal')
        ax.set_xlim(0, max_val)
        ax.set_ylim(0, max_val)
    
    ax.set_xlabel('Birth', fontsize=12)
    ax.set_ylabel('Death', fontsize=12)
    ax.set_title('H0 (Connected Components)', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    # Plot H1
    ax = axes[1]
    if len(h1_points) > 0:
        # Filter out infinite death times
        finite_h1 = h1_points[np.isfinite(h1_points[:, 1])]
        infinite_h1 = h1_points[~np.isfinite(h1_points[:, 1])]
        
        if len(finite_h1) > 0:
            ax.scatter(finite_h1[:, 0], finite_h1[:, 1], alpha=0.6, s=50, 
                      c='orange', label='Finite')
            max_val = max(finite_h1[:, 0].max(), finite_h1[:, 1].max())
        else:
            max_val = 1
            
        if len(infinite_h1) > 0:
            ax.scatter(infinite_h1[:, 0], [max_val * 1.1] * len(infinite_h1), 
                      marker='^', s=100, c='red', label='Infinite', alpha=0.6)
            max_val = max_val * 1.2
        
        # Plot diagonal
        ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='Diagonal')
        ax.set_xlim(0, max_val)
        ax.set_ylim(0, max_val)
    
    ax.set_xlabel('Birth', fontsize=12)
    ax.set_ylabel('Death', fontsize=12)
    ax.set_title('H1 (Loops/Cycles)', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    plt.suptitle(title, fontsize=16, y=1.02)
    plt.tight_layout()
    plt.show()


def compute_landscape_metrics(t_vals, landscape):
    """
    Compute peak_x, peak_y, and auc for a single landscape function.
    
    Parameters:
    -----------
    t_vals : array
        The t values (filtration values)
    landscape : array
        The landscape values
    
    Returns:
    --------
    dict : {'peak_x': float, 'peak_y': float, 'auc': float}
    """
    # Find peak
    peak_idx = np.argmax(landscape)
    peak_x = t_vals[peak_idx]
    peak_y = landscape[peak_idx]
    
    # Compute AUC using trapezoidal rule
    auc = np.trapz(landscape, t_vals)
    
    return {'peak_x': peak_x, 'peak_y': peak_y, 'auc': auc}


def create_landscape_summary_csv(landscape_curves_dict, global_k_by_pressure, output_filename='landscape_metrics_summary.csv'):
    """
    Create a CSV file with landscape metrics for all datasets.
    
    Parameters:
    -----------
    landscape_curves_dict : dict
        Dictionary of {dataset_name: {dim: (t_vals, landscapes)}}
    global_k_by_pressure : dict
        Dictionary of {pressure: {dim: k_value}} specifying number of landscapes per dimension per pressure
    output_filename : str
        Name of the output CSV file
    """
    # Find the maximum k for each dimension across all pressures (for column creation)
    max_k_per_dim = {}
    for pressure, k_dict in global_k_by_pressure.items():
        for dim, k in k_dict.items():
            if dim not in max_k_per_dim:
                max_k_per_dim[dim] = k
            else:
                max_k_per_dim[dim] = max(max_k_per_dim[dim], k)
    
    # Create column names based on maximum k across all pressures
    columns = ['dataset', 'pressure']
    
    for dim in sorted(max_k_per_dim.keys()):
        k = max_k_per_dim[dim]
        for landscape_num in range(1, k + 1):
            columns.extend([
                f'dim{dim}_landscape{landscape_num}_peak_x',
                f'dim{dim}_landscape{landscape_num}_peak_y',
                f'dim{dim}_landscape{landscape_num}_auc'
            ])
    
    # Collect data rows
    rows = []
    
    for dataset_name, landscape_data in landscape_curves_dict.items():
        # Parse dataset name to extract base name and pressure
        # Format is like 'm1053007_p1'
        if '_p' in dataset_name:
            base_name, pressure = dataset_name.rsplit('_p', 1)
        else:
            base_name = dataset_name
            pressure = ''
        
        row = {'dataset': base_name, 'pressure': pressure}
        
        # Get the k values for this specific pressure
        k_for_pressure = global_k_by_pressure.get(pressure, {})
        
        # Process each dimension
        for dim in sorted(max_k_per_dim.keys()):
            max_k = max_k_per_dim[dim]
            k_this_pressure = k_for_pressure.get(dim, 0)
            
            if dim in landscape_data:
                t_vals, landscapes = landscape_data[dim]
                
                # Process each landscape up to max_k
                for landscape_num in range(max_k):
                    # Only extract metrics if this landscape exists AND is within k for this pressure
                    if landscape_num < k_this_pressure and landscape_num < landscapes.shape[0]:
                        metrics = compute_landscape_metrics(t_vals, landscapes[landscape_num])
                        row[f'dim{dim}_landscape{landscape_num+1}_peak_x'] = metrics['peak_x']
                        row[f'dim{dim}_landscape{landscape_num+1}_peak_y'] = metrics['peak_y']
                        row[f'dim{dim}_landscape{landscape_num+1}_auc'] = metrics['auc']
                    else:
                        # If this landscape is beyond k for this pressure, fill with NaN
                        row[f'dim{dim}_landscape{landscape_num+1}_peak_x'] = np.nan
                        row[f'dim{dim}_landscape{landscape_num+1}_peak_y'] = np.nan
                        row[f'dim{dim}_landscape{landscape_num+1}_auc'] = np.nan
            else:
                # If dimension doesn't exist, fill all landscapes with NaN
                for landscape_num in range(max_k):
                    row[f'dim{dim}_landscape{landscape_num+1}_peak_x'] = np.nan
                    row[f'dim{dim}_landscape{landscape_num+1}_peak_y'] = np.nan
                    row[f'dim{dim}_landscape{landscape_num+1}_auc'] = np.nan
        
        rows.append(row)
    
    # Create DataFrame and save to CSV
    df = pd.DataFrame(rows, columns=columns)
    df.to_csv(output_filename, index=False)
    print(f"\nSaved landscape metrics summary to {output_filename}")
    
    return df


# Dictionary to store all landscape curves across all pressures
all_landscape_curves = {}
# Store global_k for each pressure
global_k_by_pressure = {}

for pressure in ['1','2','3','4']:
    datasets = {
        'm1053007': nodesToArrayH('m1p'+pressure+'_053007'),
        'm2053007': nodesToArrayH('m2p'+pressure+'_053007'),
        'm1053107': nodesToArrayH('m1p'+pressure+'_053107'),
        'm2053107': nodesToArrayH('m2p'+pressure+'_053107'),
        'm1060107': nodesToArrayH('m1p'+pressure+'_060107'),
        'm1060407': nodesToArrayC('m1p'+pressure+'_060407'),
        'm2060407': nodesToArrayC('m2p'+pressure+'_060407'),
        'm3060407': nodesToArrayC('m3p'+pressure+'_060407'),
        'm1060507': nodesToArrayC('m1p'+pressure+'_060507'),
        'm2060507': nodesToArrayC('m2p'+pressure+'_060507'),
        'm3060507': nodesToArrayC('m3p'+pressure+'_060507'),
        'm2060607': nodesToArrayC('m2p'+pressure+'_060607'),
        'm3060607': nodesToArrayC('m3p'+pressure+'_060607'),
    }

    # Example usage across all datasets
    persistence_results = {
        name: compute_persistence(points)
        for name, points in datasets.items()
    }

    # Separate results
    persistence_diagrams = {name: res[0] for name, res in persistence_results.items()}
    last_deaths = {name: res[1] for name, res in persistence_results.items()}
    # Compute global max death across all datasets for each dimension
    global_last_deaths = {}
    for name, deaths_by_dim in last_deaths.items():
        for dim, death in deaths_by_dim.items():
            if dim not in global_last_deaths:
                global_last_deaths[dim] = death
            else:
                global_last_deaths[dim] = max(global_last_deaths[dim], death)

    """for name, diagram in persistence_diagrams.items():
        plot_persistence_diagram(diagram, title=f"Persistence Diagram - {name}")"""
    # Step 1: Find global k using elbow method
    global_k = find_global_k_per_dimension(
        persistence_diagrams,
        t_min=0,
        t_max=max(global_last_deaths.values()),
        max_dim=2,
        resolution=1000, energy_threshold=.95
    )

    # Step 3: Compute landscapes with elbow-selected k
    landscape_curves = {}
    for name, diag in persistence_diagrams.items():
        landscape_curves[name] = compute_landscape_curves(
            diag,
            t_min=0,
            t_max=max(global_last_deaths.values()),
            max_dim=2,
            resolution=500,
            k_values=global_k
        )
    
    # Store global_k for this pressure
    global_k_by_pressure[pressure] = global_k
    
    # Store landscape curves with pressure info for CSV generation later
    for name, curves in landscape_curves.items():
        full_name = name + '_p' + pressure
        all_landscape_curves[full_name] = curves
        
    keys = list(landscape_curves.keys())
    graph = 'Persistence_Landscape'
    for i in range(13):
        mouse_name = keys[i]
        for dim in [0, 1, 2]:
            if i< 5:
                filename = 'fullLungGraphs/Pressure'+str(pressure)+'/'+graph+'_'+str(dim)+'_Hyper_'+mouse_name+'_P'+str(pressure)+'_fullLung'
            else:
                    filename = 'fullLungGraphs/Pressure'+str(pressure)+'/'+graph+'_'+str(dim)+'_Control_'+mouse_name+'_P'+str(pressure)+'_fullLung'
            t_vals, landscapes = landscape_curves[mouse_name][dim]
            save_landscape(
                t_vals,
                landscapes,
                dim=dim,
                filename=filename,
                title=f'Persistence Landscape {mouse_name} Dimension {dim} Full Lung',
                limit=[0, global_last_deaths[dim]],
                dpi=150
            )
    
    print("Completed Pressure "+pressure)

# After your loop finishes and all_landscape_curves is populated:
generate_all_landscape_heatmaps(
    all_landscape_curves, 
    global_k_by_pressure,
    output_base_dir='landscape_heatmaps'
)

# After all pressures are processed, create the summary CSV
print("\n=== Creating Landscape Metrics Summary CSV ===")
summary_df = create_landscape_summary_csv(
    all_landscape_curves, 
    global_k_by_pressure, 
    output_filename='landscape_metrics_summary_full.csv'
)

print("\n=== Summary Statistics ===")
print(f"Total datasets processed: {len(all_landscape_curves)}")
print("Landscapes per dimension by pressure:")
for pressure, k_dict in global_k_by_pressure.items():
    print(f"  Pressure {pressure}: {k_dict}")
print("\nFirst few rows of the summary:")
print(summary_df.head())