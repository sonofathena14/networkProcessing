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
import warnings
warnings.filterwarnings('ignore', message='.*tight_layout.*')

#Global functions
def nodesExtractor(name): #extracts nodes and their corresponding information
    file_path = 'Complement/' + name +'.csv'
    # Extract the 'connectivity' field from the 'Data' structured array
    points = pd.read_csv(file_path)
    points = points.iloc[:, 1:]
    # Save the DataFrame to inspect it
    return points

def nodesToArray(name):
    nodes_df = nodesExtractor(name)
    nodes_loc = nodes_df.loc[:,['x','y','z']]
    loc_array = nodes_loc.to_numpy()
    return loc_array

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

def compute_pairwise_l2_distances(curve_dicts, names, dims=[0, 1, 2], normalize=True, 
                                 group1_indices=None, group2_indices=None):
    """
    Compute pairwise L2 distances between all curves for each dimension.
    NOW INCLUDES distances to group averages and overall average.
    
    Args:
        curve_dicts: list of curve dictionaries {dim: (grid, curve)}
        names: list of names for each curve (e.g., mouse IDs)
        dims: list of dimensions to compute (default [0, 1, 2])
        normalize: if True, normalize by sqrt(n) for per-point distance
        group1_indices: indices for group 1 (for computing group 1 average)
        group2_indices: indices for group 2 (for computing group 2 average)
    
    Returns:
        dict of {dim: distance_matrix} where distance_matrix includes average columns
        extended_names: names with added average labels
    """
    n_curves = len(curve_dicts)
    
    # Compute group averages if indices provided
    avg_group1 = None
    avg_group2 = None
    avg_all = None
    
    if group1_indices is not None:
        group1_curves = [curve_dicts[i] for i in group1_indices]
        avg_group1 = average_curves(group1_curves)
    
    if group2_indices is not None:
        group2_curves = [curve_dicts[i] for i in group2_indices]
        avg_group2 = average_curves(group2_curves)
    
    # Always compute overall average
    avg_all = average_curves(curve_dicts)
    
    # Determine number of extra columns
    n_extra = 0
    extra_names = []
    if avg_group1 is not None:
        n_extra += 1
        extra_names.append('Hyper_Avg')
    if avg_group2 is not None:
        n_extra += 1
        extra_names.append('Control_Avg')
    n_extra += 1  # Always add overall average
    extra_names.append('Avg_All')
    
    extended_names = names + extra_names
    
    distance_matrices = {}
    
    for dim in dims:
        # Initialize distance matrix with extra columns
        dist_matrix = np.zeros((n_curves, n_curves + n_extra))
        
        # Compute pairwise distances (original logic)
        for i in range(n_curves):
            for j in range(n_curves):
                if i == j:
                    dist_matrix[i, j] = 0
                elif dim in curve_dicts[i] and dim in curve_dicts[j]:
                    _, curve_i = curve_dicts[i][dim]
                    _, curve_j = curve_dicts[j][dim]
                    
                    # Compute L2 distance
                    l2_dist = np.linalg.norm(curve_i - curve_j)
                    
                    if normalize:
                        l2_dist = l2_dist / np.sqrt(len(curve_i))
                    
                    dist_matrix[i, j] = l2_dist
                else:
                    # If dimension doesn't exist in one of the curves
                    dist_matrix[i, j] = np.nan
        
        # Add columns for distances to averages
        col_idx = n_curves
        
        if avg_group1 is not None:
            for i in range(n_curves):
                if dim in curve_dicts[i] and dim in avg_group1:
                    _, curve_i = curve_dicts[i][dim]
                    _, avg_curve = avg_group1[dim]
                    
                    l2_dist = np.linalg.norm(curve_i - avg_curve)
                    if normalize:
                        l2_dist = l2_dist / np.sqrt(len(curve_i))
                    
                    dist_matrix[i, col_idx] = l2_dist
                else:
                    dist_matrix[i, col_idx] = np.nan
            col_idx += 1
        
        if avg_group2 is not None:
            for i in range(n_curves):
                if dim in curve_dicts[i] and dim in avg_group2:
                    _, curve_i = curve_dicts[i][dim]
                    _, avg_curve = avg_group2[dim]
                    
                    l2_dist = np.linalg.norm(curve_i - avg_curve)
                    if normalize:
                        l2_dist = l2_dist / np.sqrt(len(curve_i))
                    
                    dist_matrix[i, col_idx] = l2_dist
                else:
                    dist_matrix[i, col_idx] = np.nan
            col_idx += 1
        
        # Always add overall average column
        for i in range(n_curves):
            if dim in curve_dicts[i] and dim in avg_all:
                _, curve_i = curve_dicts[i][dim]
                _, avg_curve = avg_all[dim]
                
                l2_dist = np.linalg.norm(curve_i - avg_curve)
                if normalize:
                    l2_dist = l2_dist / np.sqrt(len(curve_i))
                
                dist_matrix[i, col_idx] = l2_dist
            else:
                dist_matrix[i, col_idx] = np.nan
        
        distance_matrices[dim] = dist_matrix
    
    return distance_matrices, extended_names


def plot_distance_heatmap(distance_matrix, names, dim, title=None, figsize=(10, 8), 
                          cmap='viridis', annot=True, fmt='.3f',
                          group_stats=None, group1_indices=None, group2_indices=None):
    """
    Plot a heatmap of pairwise distances with optional group statistics.
    NOW HANDLES extended distance matrix with average columns.
    """
    if group_stats is not None:
        # Create figure with two subplots: heatmap and stats
        fig = plt.figure(figsize=(figsize[0] + 4, figsize[1]))
        gs = fig.add_gridspec(1, 2, width_ratios=[3, 1], wspace=0.3)
        ax_heat = fig.add_subplot(gs[0])
        ax_stats = fig.add_subplot(gs[1])
    else:
        fig, ax_heat = plt.subplots(figsize=figsize)
    
    # Plot heatmap
    sns.heatmap(distance_matrix, 
                xticklabels=names, 
                yticklabels=names[:distance_matrix.shape[0]],  # Only row names
                cmap=cmap, 
                annot=annot, 
                fmt=fmt,
                cbar_kws={'label': 'L2 Distance'},
                ax=ax_heat)
    
    ax_heat.set_title(title or f'Pairwise L2 Distances (Dimension {dim})', fontsize=16)
    ax_heat.set_xlabel('Curves', fontsize=12)
    ax_heat.set_ylabel('Curves', fontsize=12)
    ax_heat.tick_params(axis='x', rotation=45)
    plt.setp(ax_heat.get_xticklabels(), rotation=45, ha='right')
    
    # Add visual separation lines between groups if indices provided
    if group1_indices is not None and group2_indices is not None:
        split_point = len(group1_indices)
        ax_heat.axhline(y=split_point, color='red', linewidth=2, linestyle='--')
        # Vertical line only at split between actual samples (not before averages)
        n_samples = distance_matrix.shape[0]
        ax_heat.axvline(x=split_point, color='red', linewidth=2, linestyle='--')
        # Add vertical line before average columns
        ax_heat.axvline(x=n_samples, color='blue', linewidth=2, linestyle='--')
    
    # Add statistics panel if provided
    if group_stats is not None:
        ax_stats.axis('off')
        
        # Calculate statistics
        within1_mean = group_stats['within_group1'].mean()
        within1_std = group_stats['within_group1'].std()
        within2_mean = group_stats['within_group2'].mean()
        within2_std = group_stats['within_group2'].std()
        
        stats_text = f"Group Statistics\n{'='*25}\n\n"
        stats_text += f"Within Group 1 (Hyper):\n"
        stats_text += f"  Mean: {within1_mean:.4f}\n"
        stats_text += f"  Std:  {within1_std:.4f}\n"
        stats_text += f"  N:    {len(group_stats['within_group1'])}\n\n"
        
        stats_text += f"Within Group 2 (Control):\n"
        stats_text += f"  Mean: {within2_mean:.4f}\n"
        stats_text += f"  Std:  {within2_std:.4f}\n"
        stats_text += f"  N:    {len(group_stats['within_group2'])}\n\n"
        
        ax_stats.text(0.05, 0.95, stats_text, 
                     transform=ax_stats.transAxes,
                     fontsize=10,
                     verticalalignment='top',
                     fontfamily='monospace',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.tight_layout()
    plt.show()


def save_distance_heatmap(distance_matrix, names, dim, filename, title=None, 
                         figsize=(10, 8), cmap='viridis', annot=True, fmt='.3f', dpi=300,
                         group_stats=None, group1_indices=None, group2_indices=None):
    """
    Save a heatmap of pairwise distances to file with optional group statistics.
    NOW HANDLES extended distance matrix with average columns.
    """
    if group_stats is not None:
        # Create figure with two subplots: heatmap and stats
        fig = plt.figure(figsize=(figsize[0] + 4, figsize[1]))
        gs = fig.add_gridspec(1, 2, width_ratios=[3, 1], wspace=0.3)
        ax_heat = fig.add_subplot(gs[0])
        ax_stats = fig.add_subplot(gs[1])
    else:
        fig, ax_heat = plt.subplots(figsize=figsize)
    
    # Plot heatmap
    sns.heatmap(distance_matrix, 
                xticklabels=names, 
                yticklabels=names[:distance_matrix.shape[0]],  # Only row names
                cmap=cmap, 
                annot=annot, 
                fmt=fmt,
                cbar_kws={'label': 'L2 Distance'},
                ax=ax_heat)
    
    ax_heat.set_title(title or f'Pairwise L2 Distances (Dimension {dim})', fontsize=16)
    ax_heat.set_xlabel('Curves', fontsize=12)
    ax_heat.set_ylabel('Curves', fontsize=12)
    ax_heat.tick_params(axis='x', rotation=45)
    plt.setp(ax_heat.get_xticklabels(), rotation=45, ha='right')
    
    # Add visual separation lines between groups if indices provided
    if group1_indices is not None and group2_indices is not None:
        split_point = len(group1_indices)
        ax_heat.axhline(y=split_point, color='red', linewidth=2, linestyle='--')
        # Vertical line only at split between actual samples (not before averages)
        n_samples = distance_matrix.shape[0]
        ax_heat.axvline(x=split_point, color='red', linewidth=2, linestyle='--')
        # Add vertical line before average columns
        ax_heat.axvline(x=n_samples, color='red', linewidth=2, linestyle='-')
    
    # Add statistics panel if provided
    if group_stats is not None:
        ax_stats.axis('off')
        
        # Calculate statistics
        within1_mean = group_stats['within_group1'].mean()
        within1_std = group_stats['within_group1'].std()
        within2_mean = group_stats['within_group2'].mean()
        within2_std = group_stats['within_group2'].std()
        
        stats_text = f"Group Statistics\n{'='*25}\n\n"
        stats_text += f"Within Group 1 (Hyper):\n"
        stats_text += f"  Mean: {within1_mean:.4f}\n"
        stats_text += f"  Std:  {within1_std:.4f}\n"
        stats_text += f"  N:    {len(group_stats['within_group1'])}\n\n"
        
        stats_text += f"Within Group 2 (Control):\n"
        stats_text += f"  Mean: {within2_mean:.4f}\n"
        stats_text += f"  Std:  {within2_std:.4f}\n"
        stats_text += f"  N:    {len(group_stats['within_group2'])}\n\n"
        
        ax_stats.text(0.05, 0.95, stats_text, 
                     transform=ax_stats.transAxes,
                     fontsize=10,
                     verticalalignment='top',
                     fontfamily='monospace',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.tight_layout()
    plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved heatmap to {filename}")


def compute_group_comparison(curve_dicts, names, group1_indices, group2_indices, 
                            dims=[0, 1, 2], normalize=True):
    """
    Compute L2 distances specifically comparing two groups.
    
    Args:
        curve_dicts: list of curve dictionaries
        names: list of names
        group1_indices: indices of group 1 (e.g., [0,1,2,3,4] for hyper)
        group2_indices: indices of group 2 (e.g., [5,6,7,8,9,10,11,12] for control)
        dims: dimensions to compute
        normalize: whether to normalize distances
    
    Returns:
        dict of {dim: (within_group1, within_group2, between_groups)}
    """
    results = {}
    
    for dim in dims:
        # Within group 1
        within1 = []
        for i in group1_indices:
            for j in group1_indices:
                if i < j and dim in curve_dicts[i] and dim in curve_dicts[j]:
                    _, curve_i = curve_dicts[i][dim]
                    _, curve_j = curve_dicts[j][dim]
                    dist = np.linalg.norm(curve_i - curve_j)
                    if normalize:
                        dist /= np.sqrt(len(curve_i))
                    within1.append(dist)
        
        # Within group 2
        within2 = []
        for i in group2_indices:
            for j in group2_indices:
                if i < j and dim in curve_dicts[i] and dim in curve_dicts[j]:
                    _, curve_i = curve_dicts[i][dim]
                    _, curve_j = curve_dicts[j][dim]
                    dist = np.linalg.norm(curve_i - curve_j)
                    if normalize:
                        dist /= np.sqrt(len(curve_i))
                    within2.append(dist)
        
        # Between groups
        between = []
        for i in group1_indices:
            for j in group2_indices:
                if dim in curve_dicts[i] and dim in curve_dicts[j]:
                    _, curve_i = curve_dicts[i][dim]
                    _, curve_j = curve_dicts[j][dim]
                    dist = np.linalg.norm(curve_i - curve_j)
                    if normalize:
                        dist /= np.sqrt(len(curve_i))
                    between.append(dist)
        
        results[dim] = {
            'within_group1': np.array(within1),
            'within_group2': np.array(within2)
            #'between_groups': np.array(between)
        }
    
    return results

def extract_intervals(persistence_diagram, dim):
    return np.array([
        (b, d)
        for d_dim, (b, d) in persistence_diagram
        if d_dim == dim and np.isfinite(d)
    ])

def compute_lifespan_curves(
    persistence_diagram,
    t_min,
    t_max,
    max_dim=2,
    resolution=500
):
    lifespan_curves = {}
    
    for dim in range(max_dim + 1):
        intervals = extract_intervals(persistence_diagram, dim)
        
        if intervals.size == 0:
            lifespan_curves[dim] = (np.zeros(resolution), np.zeros(resolution))
            continue
        
        births = intervals[:, 0]
        deaths = intervals[:, 1]
        lifespans = deaths - births
        
        # Use GLOBAL range, not per-dimension range
        t_vals = np.linspace(t_min, t_max, resolution)
        LC = np.zeros_like(t_vals)
        
        for i, t in enumerate(t_vals):
            alive = (births <= t) & (t < deaths)
            LC[i] = lifespans[alive].sum()
        
        lifespan_curves[dim] = (t_vals, LC)
    
    return lifespan_curves

def compute_betti_curve(diag, max_dim=2, resolution=100000, cutoff=None):
    betti_curves = {}
    
    # Determine global filtration range across ALL dimensions
    all_values = []
    for dim, (birth, death) in diag:
        all_values.append(birth)
        if np.isfinite(death):
            all_values.append(death)
    
    global_min = np.min(all_values) if all_values else 0
    global_max = np.max(all_values) if all_values else 1
    
    if cutoff is not None:
        global_max = cutoff
    
    for dim in range(max_dim + 1):
        diag_dim = np.array([pt[1] for pt in diag if pt[0] == dim])
        
        if len(diag_dim) == 0:
            betti_curves[dim] = (np.zeros(resolution), np.zeros(resolution))
            continue
        
        # Use GLOBAL range, not per-dimension range
        grid = np.linspace(global_min, global_max, resolution)
        curve = np.zeros_like(grid)
        
        for birth, death in diag_dim:
            if np.isinf(death):
                # Infinite bar: alive from birth onward
                curve += (grid >= birth)
            else:
                # Finite bar: alive in [birth, death]
                curve += (grid >= birth) & (grid <= death)
        
        betti_curves[dim] = (grid, curve)
    
    return betti_curves

def compute_norm_lifespan_curves(
    persistence_diagram,
    t_min,
    t_max,
    max_dim=2,
    resolution=500
):
    lifespan_curves = {}
    
    for dim in range(max_dim + 1):
        intervals = extract_intervals(persistence_diagram, dim)
        
        if intervals.size == 0:
            lifespan_curves[dim] = (np.zeros(resolution), np.zeros(resolution))
            continue
        
        births = intervals[:, 0]
        deaths = intervals[:, 1]
        lifespans = deaths - births
        
        # Compute total lifespan for normalization
        total_lifespan = lifespans.sum()
        
        # Use GLOBAL range, not per-dimension range
        t_vals = np.linspace(t_min, t_max, resolution)
        LC = np.zeros_like(t_vals)
        
        for i, t in enumerate(t_vals):
            alive = (births <= t) & (t < deaths)
            # Normalized: sum of lifespans alive at t / total lifespan
            if total_lifespan > 0:
                LC[i] = lifespans[alive].sum() / total_lifespan
            else:
                LC[i] = 0
        
        lifespan_curves[dim] = (t_vals, LC)
    
    return lifespan_curves

def run_pairwise(curves, graphtype, lobe, pressure, all_metrics):
    # Compute pairwise distances for all mice
    keys = list(curves.keys())
    curve_list = [curves[k] for k in keys]
    
    # Define group indices
    group1_indices = list(range(5))  # First 5 are hyper
    group2_indices = list(range(5, 13))  # Last 8 are control
    
    # MODIFIED: Pass group indices to compute distances to averages
    distance_matrices, names = compute_pairwise_l2_distances(
        curve_list, 
        keys, 
        dims=[0, 1, 2],
        normalize=True,
        group1_indices=group1_indices,
        group2_indices=group2_indices
    )
    
    # Compare groups specifically
    group_stats = compute_group_comparison(
        curve_list,
        keys,
        group1_indices=group1_indices,
        group2_indices=group2_indices,
        dims=[0, 1, 2],
        normalize=True
    )
    
    for dim in [0,1,2]:
        # Save heatmap visualization
        save_distance_heatmap(distance_matrices[dim], names, dim, 
            'PairwiseDistances/Lobe/'+graphtype+'_'+lobe+'_Pressure'+pressure+'_Dim'+str(dim), 
            title=f"{graphtype} Curve Pairwise $L^2$ Distances (Dimension {dim}, Lobe {lobe})", 
            figsize=(12, 8), cmap='viridis', annot=True, fmt='.3f', dpi=300,
            group_stats=group_stats[dim], group1_indices=group1_indices, group2_indices=group2_indices)
        
        # CHANGE 5: Save distance matrix as CSV
        csv_filename = 'PairwiseDistances/Lobe/'+graphtype+'_'+lobe+'_Pressure'+pressure+'_Dim'+str(dim)+'.csv'
        df_distances = pd.DataFrame(
            distance_matrices[dim],
            index=keys,  # Row names are the original mouse IDs
            columns=names  # Column names include mice + averages
        )
        df_distances.to_csv(csv_filename)
        print(f"Saved distance matrix to {csv_filename}")
    
    # Compute metrics for each dataset and add to all_metrics
    for name, curve_dict in curves.items():
        metrics = compute_curve_metrics(curve_dict, dims=[0, 1, 2])
        
        # Create key for this dataset-pressure-lobe combination
        key = (name, pressure, lobe)
        if key not in all_metrics:
            all_metrics[key] = {'dataset': name, 'pressure': pressure, 'lobe': lobe}
        
        # Add metrics with curve type suffix
        for metric_name, metric_value in metrics.items():
            all_metrics[key][f'{metric_name}_{graphtype}'] = metric_value

def compute_curve_metrics(curve_dict, dims=[0, 1, 2]):
    """
    Compute peak x, peak y, and area under curve for each dimension.
    
    Returns:
        dict with keys like 'dim0_peak_x', 'dim0_peak_y', 'dim0_auc', etc.
    """
    metrics = {}
    
    for dim in dims:
        if dim in curve_dict:
            grid, curve = curve_dict[dim]
            
            if len(grid) > 0 and len(curve) > 0:
                # Peak y value and its x location
                peak_idx = np.argmax(curve)
                peak_x = grid[peak_idx]
                peak_y = curve[peak_idx]
                
                # Area under curve using trapezoidal rule
                auc = np.trapezoid(curve, grid)
                
                metrics[f'dim{dim}_peak_x'] = peak_x
                metrics[f'dim{dim}_peak_y'] = peak_y
                metrics[f'dim{dim}_auc'] = auc
            else:
                metrics[f'dim{dim}_peak_x'] = np.nan
                metrics[f'dim{dim}_peak_y'] = np.nan
                metrics[f'dim{dim}_auc'] = np.nan
        else:
            metrics[f'dim{dim}_peak_x'] = np.nan
            metrics[f'dim{dim}_peak_y'] = np.nan
            metrics[f'dim{dim}_auc'] = np.nan
    
    return metrics

# Initialize dictionary to collect all metrics
all_metrics = {}
for pressure in ['1','2','3','4']:
    for lobe in ['left','middle','superior','inferior','postcaval']:
        print('Starting: Pressure '+pressure+' Lobe '+lobe)
        datasets = {
            'm1053007': nodesToArray('m1p'+pressure+'_053007_'+lobe),
            'm2053007': nodesToArray('m2p'+pressure+'_053007_'+lobe),
            'm1053107': nodesToArray('m1p'+pressure+'_053107_'+lobe),
            'm2053107': nodesToArray('m2p'+pressure+'_053107_'+lobe),
            'm1060107': nodesToArray('m1p'+pressure+'_060107_'+lobe),
            'm1060407': nodesToArray('m1p'+pressure+'_060407_'+lobe),
            'm2060407': nodesToArray('m2p'+pressure+'_060407_'+lobe),
            'm3060407': nodesToArray('m3p'+pressure+'_060407_'+lobe),
            'm1060507': nodesToArray('m1p'+pressure+'_060507_'+lobe),
            'm2060507': nodesToArray('m2p'+pressure+'_060507_'+lobe),
            'm3060507': nodesToArray('m3p'+pressure+'_060507_'+lobe),
            'm2060607': nodesToArray('m2p'+pressure+'_060607_'+lobe),
            'm3060607': nodesToArray('m3p'+pressure+'_060607_'+lobe)
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

        
        lifespan_curves = {}
        for name, diag in persistence_diagrams.items():
            points = datasets[name]
            lifespan_curves[name] = compute_lifespan_curves(diag, t_min = 0, t_max=global_last_deaths[dim], resolution = 100000)

        betti_curves = {}
        for name, diag in persistence_diagrams.items():
            points = datasets[name]
            betti_curves[name] = compute_betti_curve(diag, cutoff=global_last_deaths[dim])

        norm_lifespan_curves = {}
        for name, diag in persistence_diagrams.items():
            points = datasets[name]
            norm_lifespan_curves[name] = compute_norm_lifespan_curves(diag, t_min = 0, t_max=global_last_deaths[dim], resolution = 100000)

        run_pairwise(lifespan_curves, 'Lifespan', lobe, pressure, all_metrics)
        run_pairwise(betti_curves, 'Betti', lobe, pressure, all_metrics)
        run_pairwise(norm_lifespan_curves, 'Norm_Lifespan', lobe, pressure, all_metrics)
        print('Completed: Pressure '+pressure+' '+lobe)

# Save all metrics to CSV
metrics_list = list(all_metrics.values())
metrics_df = pd.DataFrame(metrics_list)

# Sort columns: first dataset, pressure, lobe (if exists), then all metrics alphabetically
base_cols = ['dataset', 'pressure']
if 'lobe' in metrics_df.columns:
    base_cols.append('lobe')
metric_cols = sorted([col for col in metrics_df.columns if col not in base_cols])
column_order = base_cols + metric_cols

metrics_df = metrics_df[column_order]
metrics_df.to_csv('PairwiseDistances/Lobe_Complement/curve_metrics_summary.csv', index=False)
print(f"\nSaved curve metrics summary to PairwiseDistances/Lobe_Complement/curve_metrics_summary.csv")
print(f"Total rows: {len(metrics_df)}")
print(f"Columns: {list(metrics_df.columns)}")