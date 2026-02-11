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

# Data extractor functions
def connectivityExtractor(name, pruned):
    if pruned == 0:
        file_path = 'Networks/Network_Vessels_' + name +'.mat'
    elif pruned == 1: 
        file_path = 'Pruned/Pruned_Network_' + name +'.mat'
    matlab_data = scipy.io.loadmat(file_path)
    # Extract the 'connectivity' field from the 'Data' structured array
    data_structure = matlab_data['Data']
    connectivity_raw = data_structure['connectivity'][0, 0]  # Access the data (adjust indexing if needed)
    # Reshape or ensure it's a proper 2D array (if required)
    connectivity_data = connectivity_raw.squeeze()
    # Create a DataFrame from the connectivity data
    connectivity_df = pd.DataFrame(connectivity_data, columns=['Parent', 'Daughter1', 'Daughter2', 'Daughter3'])
    connectivity_df.replace(0, np.nan, inplace=True) #ensure all nonexistent vessels have NaN
    connectivity_df.at[0,'Parent']=0 #make sure first vessel is 0 (purposefully removed in last step for ease)
    # Save the DataFrame to inspect it
    return connectivity_df

def nodesExtractor(name, pruned): #extracts nodes and their corresponding information
    if pruned == 0:
        file_path = 'Networks/Network_Vessels_' + name +'.mat'
    elif pruned == 1: 
        file_path = 'Pruned/Pruned_Network_' + name +'.mat'
    matlab_data = scipy.io.loadmat(file_path)
    # Extract the 'connectivity' field from the 'Data' structured array
    if pruned == 0:
        data_structure = matlab_data['nodesC2']
    elif pruned == 1:
        data_structure = matlab_data['nodesC3']
    # Reshape or ensure it's a proper 2D array (if required)
    nodes_data = data_structure.squeeze()
    # Create a DataFrame from the connectivity data
    nodes_df = pd.DataFrame(nodes_data, columns=['NodeID', 'X', 'Y', 'Z', 'Degree'])
    # Save the DataFrame to inspect it
    return nodes_df

def edgesExtractor(name, pruned): #extracts segments to create a dataframe of from and to nodes
    if pruned == 0:
        file_path = 'Networks/Network_Vessels_' + name +'.mat'
    elif pruned == 1: 
        file_path = 'Pruned/Pruned_Network_' + name +'.mat'
    matlab_data = scipy.io.loadmat(file_path)
    # Extract the 'segments' field
    data_structure = matlab_data['segments']
    # Reshape or ensure it's a proper 2D array (if required)
    edges_data = data_structure.squeeze()
    # Create a DataFrame from the connectivity data
    edge_df = pd.DataFrame(edges_data, columns=['ID', 'From', 'To'])
    # Save the DataFrame to inspect it
    return edge_df
    
def findInputVessel(segments,fromnode,to):
    vessel = segments[((segments['From'] == fromnode)&(segments['To']==to))|((segments['From'] == to)&(segments['To']==fromnode))]
    return int(vessel['ID'].iloc[0])

def mapIDExtractor(name, pruned):
    if pruned == 0:
        file_path = 'Networks/Network_Vessels_' + name +'.mat'
    elif pruned == 1: 
        file_path = 'Pruned/Pruned_Network_' + name +'.mat'
    matlab_data = scipy.io.loadmat(file_path)
    # Extract the 'mapID' field from the 'Data' structured array
    data_structure = matlab_data['Data']
    map_raw = data_structure['mapIDs'][0, 0]  # Access the data (adjust indexing if needed)
    # Reshape or ensure it's a proper 2D array (if required)
    map_data = map_raw.squeeze()
    # Create a DataFrame from the connectivity data
    map_df = pd.DataFrame(map_data, columns=['New', 'Old'])
    # Save the DataFrame to inspect it
    return map_df

def lobeExtractor(name, vesID,pruned):
    data = connectivityExtractor(name,pruned)
    
    tree = defaultdict(list)
    for _,row in data.iterrows():
        parent = row['Parent']
        for daughter_col in ['Daughter1','Daughter2','Daughter3']:
            daughter = row[daughter_col]
            if pd.notna(daughter):
                tree[parent].append(daughter)

    visited = set()
    queue = deque([vesID])

    while queue:
        current = queue.popleft()
        if current not in visited:
            visited.add(current)
            queue.extend(tree.get(current,[]))
    
    visited.discard(vesID)  # Remove vesID from visited
    downstream_df = data[data['Parent'].isin(visited)]
    return downstream_df

def node_loc(name,lobe_nodes,pruned):
    nodes = nodesExtractor(name,pruned)
    lobe = nodes[nodes['NodeID'].isin(lobe_nodes)]
    return lobe[['X','Y','Z']]

def lobeTermLoc(name,fromnode,tonode,pruned):
    segments = edgesExtractor(name, pruned)
    maps = mapIDExtractor(name, pruned)
    vesID = findInputVessel(segments,fromnode,tonode)
    newID = int(maps[maps['Old']==vesID]['New'].iloc[0])
    lobe_ves = lobeExtractor(name,newID, pruned)
    new_lobe_ves_ID = lobe_ves['Parent'].to_numpy()
    oldID = maps[maps['New'].isin(new_lobe_ves_ID)]['Old'].to_numpy()
    fromnodes = segments[segments['ID'].isin(oldID)]['From'].to_numpy()
    tonodes = segments[segments['ID'].isin(oldID)]['To'].to_numpy()
    lobe_nodes = np.unique(np.concatenate((fromnodes,tonodes))).astype(int)
    lobe_node_loc = node_loc(name,lobe_nodes,pruned)/1000
    return lobe_node_loc.to_numpy()

def twoInputVessels(name,from1,to1,from2,to2,pruned):
    points1 = lobeTermLoc(name,from1,to1,pruned)
    points2 = lobeTermLoc(name,from2,to2,pruned)
    return np.unique(np.concatenate((points1,points2),axis=0),axis=0)

    #TDA functions
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

        if lobe == "middle":
            if pressure == '1':
                m1053007 = [1841, 1945, 1841, 1964]
                m2053007 = [2463, 2497, 2502, 2501]
                m1053107 = [4179, 4304, 4018, 4379]
                m2053107 = [2920, 2921, 3150, 3149]
                m1060107 = [2766, 2744, 2992, 1332]
                m1060407 = [2456, 2459, 2455, 1061]
                m2060407 = [2963, 2568, 2964, 2995]
                m3060407 = [2509, 2613, 2507, 408]
                m1060507 = [1992, 1929, 2030, 2065]
                m2060507 = [2390, 2227, 2284, 2283]
                m3060507 = [2624, 2377, 2568, 1365]
                m2060607 = [1916, 2052, 2057, 1991]
                m3060607 = [2604, 2466, 2677, 2607]
            if pressure == '2':
                m1053007 = [1841, 1957, 1841, 2025]
                m2053007 = [1831, 1854, 1908, 1919]
                m1053107 = [2738, 2739, 2959, 1567]
                m2053107 = [1645, 1644, 1642, 1641]
                m1060107 = [2164, 2098, 2153, 927]
                m1060407 = [1289, 1269, 1253, 1252]
                m2060407 = [1464, 1466, 1465, 592]
                m3060407 = [1907, 1962, 1939, 337]
                m1060507 = [1322, 1400, 1320, 1425]
                m2060507 = [1565, 1566, 1578, 1577]
                m3060507 = [1732, 1589, 1749, 1109]
                m2060607 = [1441, 1221, 1460, 221]
                m3060607 = [1925, 1946, 2002, 2094]
            if pressure == '3':
                m1053007 = [1721, 1720, 1721, 95]
                m2053007 = [1731, 1737, 1682, 1733]
                m1053107 = [3103, 2859, 2892, 844]
                m2053107 = [1621, 1622, 1714, 1709]
                m1060107 = [1651, 1727, 1611, 1610]
                m1060407 = [1110, 1109, 1140, 426]
                m2060407 = [854, 283, 891, 50]
                m3060407 = [968, 938, 967, 976]
                m1060507 = [1470, 1475, 1469, 1370]
                m2060507 = [1041, 1121, 1042, 1163]
                m3060507 = [1000, 451, 1062, 1016]
                m2060607 = [607, 595, 616, 640]
                m3060607 = [1261, 1202, 1257, 1268]
            if pressure == '4':
                m1053007 = [1954, 1803, 1954, 1788]
                m2053007 = [1223, 161, 1273, 1342]
                m1053107 = [2059, 2186, 2058, 1054]
                m2053107 = [1520, 1532, 1555, 1565]
                m1060107 = [1652, 1619, 1653, 1529]
                m1060407 = [1055, 1054, 1084, 516]
                m2060407 = [692, 127, 690, 313]
                m3060407 = [865, 61, 878, 880]
                m1060507 = [1243, 1244, 1212, 1341]
                m2060507 = [944, 993, 945, 969]
                m3060507 = [959, 960, 1003, 1093]
                m2060607 = [753, 768, 720, 719]
                m3060607 = [1290, 270, 1207, 1309]
        if pressure == '1':
            if lobe == 'left':
                m1053007 = [1831, 1858]
                m2053007 = [2367, 2368]
                m1053107 = [3954, 3900]
                m2053107 = [2866, 2868]
                m1060107 = [190, 2722]
                m1060407 = [2272, 2273]
                m2060407 = [2776, 2774]
                m3060407 = [2473, 2472]
                m1060507 = [2121, 2040]
                m2060507 = [2257, 2258]
                m3060507 = [692, 2524]
                m2060607 = [1475, 1997]
                m3060607 = [53, 2576]
            if lobe == 'superior':
                m1053007 = [1836, 1835]
                m2053007 = [2464, 2406]
                m1053107 = [4071, 685]
                m2053107 = [2867, 2979]
                m1060107 = [2780, 2716, 2780, 2781]
                m1060407 = [2274, 2283]
                m2060407 = [2742, 2598]
                m3060407 = [2418, 2419, 2420, 2491]
                m1060507 = [1993, 1997]
                m2060507 = [2259, 2330]
                m3060507 = [2392, 2398]
                m2060607 = [1914, 1915]
                m3060607 = [2603, 2518]
            if lobe == 'inferior':
                m1053007 = [1841, 1839]
                m2053007 = [2502, 2501]
                m1053107 = [4018, 4019]
                m2053107 = [3150, 3170]
                m1060107 = [2992, 2895]
                m1060407 = [2455, 2394]
                m2060407 = [2964, 2727]
                m3060407 = [2507, 2508]
                m1060507 = [2030, 2027]
                m2060507 = [2391, 2423]
                m3060507 = [2537, 1693]
                m2060607 = [2057, 2056]
                m3060607 = [2677, 2676]
            if lobe == 'postcaval':
                m1053007 = [1864, 35]
                m2053007 = [2465, 2558]
                m1053107 = [4296, 1851]
                m2053107 = [2915, 2913]
                m1060107 = [2766, 2808]
                m1060407 = [2454, 2602]
                m2060407 = [2886, 2962]
                m3060407 = [2421, 2469]
                m1060507 = [1994, 2178]
                m2060507 = [2284, 2454]
                m3060507 = [2628, 2355]
                m2060607 = [1916, 1970]
                m3060607 = [2701, 2688]
        if pressure == '2':
            if lobe == 'left':
                m1053007 = [1075, 1960]
                m2053007 = [217, 1868]
                m1053107 = [1506, 2897]
                m2053107 = [44, 1616]
                m1060107 = [2089, 2169]
                m1060407 = [1225, 1265]
                m2060407 = [1371, 1392]
                m3060407 = [1882, 1924]
                m1060507 = [118, 1381]
                m2060507 = [1509, 1519]
                m3060507 = [51, 1591]
                m2060607 = [1495, 1419]
                m3060607 = [535, 1916]
            if lobe == 'superior':
                m1053007 = [1829, 1828]
                m2053007 = [1829, 1830]
                m1053107 = [2909, 2895]
                m2053107 = [1640, 1693]
                m1060107 = [2162, 2163]
                m1060407 = [1266, 1239]
                m2060407 = [1387, 1405]
                m3060407 = [1883, 1890, 1880, 1881]
                m1060507 = [1468, 1454]
                m2060507 = [1528, 1552]
                m3060507 = [1706, 1705]
                m2060607 = [1440, 1426]
                m3060607 = [1923, 1924]
            if lobe == 'inferior':
                m1053007 = [1842, 2020]
                m2053007 = [1908, 1979]
                m1053107 = [2959, 2957]
                m2053107 = [1642, 1755]
                m1060107 = [2153, 2151]
                m1060407 = [1253, 1252]
                m2060407 = [1465, 1515]
                m3060407 = [1939, 1982]
                m1060507 = [1320, 1424]
                m2060507 = [1620, 1602]
                m3060507 = [1749, 1809]
                m2060607 = [1460, 1375]
                m3060607 = [2002, 1903]
            if lobe == 'postcaval':
                m1053007 = [1862, 2003]
                m2053007 = [1907, 1909]
                m1053107 = [3056, 1829]
                m2053107 = [1724, 1722]
                m1060107 = [2164, 2190]
                m1060407 = [1290, 1323]
                m2060407 = [1386, 1410]
                m3060407 = [1907, 1796]
                m1060507 = [1319, 1321]
                m2060507 = [1578, 1650]
                m3060507 = [1739, 1689]
                m2060607 = [1441, 1509]
                m3060607 = [1951, 1950]
        if pressure == '3':
            if lobe == 'left':
                m1053007 = [1630, 1628]
                m2053007 = [1619, 1621]
                m1053107 = [246, 2864]
                m2053107 = [717, 1751]
                m1060107 = [318, 1645]
                m1060407 = [1122, 1121]
                m2060407 = [808, 814]
                m3060407 = [919, 921]
                m1060507 = [80, 1354]
                m2060507 = [1097, 1116]
                m3060507 = [191, 959]
                m2060607 = [653, 620]
                m3060607 = [142, 1171]
            if lobe == 'superior':
                m1053007 = [1657, 1656]
                m2053007 = [1620, 1703]
                m1053107 = [2982, 2955]
                m2053107 = [1612, 1610]
                m1060107 = [1644, 1650]
                m1060407 = [1126, 1095]
                m2060407 = [862, 832]
                m3060407 = [969, 970, 920, 924] 
                m1060507 = [1405, 1461]
                m2060507 = [1096, 1095]
                m3060507 = [952, 950]
                m2060607 = [615, 642]
                m3060607 = [1172, 1255]
            if lobe == 'inferior':
                m1053007 = [1721, 1624]
                m2053007 = [1682, 1683]
                m1053107 = [2892, 1313]
                m2053107 = [1714, 1715]
                m1060107 = [1611, 1757]
                m1060407 = [1140, 1157]
                m2060407 = [891, 844]
                m3060407 = [967, 975]
                m1060507 = [1469, 1488]
                m2060507 = [1040, 1043]
                m3060507 = [1062, 1063]
                m2060607 = [616, 637]
                m3060607 = [1257, 1284]
            if lobe == 'postcaval':
                m1053007 = [1755, 1636]
                m2053007 = [1704, 1732]
                m1053107 = [3039, 1474]
                m2053107 = [1655, 1653]
                m1060107 = [1651, 29]
                m1060407 = [1134, 1141]
                m2060407 = [851, 853]
                m3060407 = [968, 1020]
                m1060507 = [1381, 1379]
                m2060507 = [1042, 1136]
                m3060507 = [1009, 1057]
                m2060607 = [607, 664]
                m3060607 = [1250, 1249]
        if pressure == '4':
            if lobe == 'left':
                m1053007 = [1810, 1811]
                m2053007 = [1209, 1214]
                m1053107 = [1252, 2103]
                m2053107 = [1239, 1392]
                m1060107 = [215, 1664]
                m1060407 = [1065, 1013]
                m2060407 = [694, 701]
                m3060407 = [851, 926]
                m1060507 = [1202, 1169]
                m2060507 = [1006, 1064]
                m3060507 = [395, 975]
                m2060607 = [511, 737]
                m3060607 = [1223, 1222]
            if lobe == 'superior':
                m1053007 = [1819, 1818]
                m2053007 = [1213, 1259]
                m1053107 = [2125, 2154]
                m2053107 = [1391, 1477]
                m1060107 = [1574, 1522]
                m1060407 = [1053, 1051]
                m2060407 = [695, 727]
                m3060407 = [885, 894, 884, 864]
                m1060507 = [1232, 1233]
                m2060507 = [999, 983]
                m3060507 = [1027, 1026]
                m2060607 = [738, 741]
                m3060607 = [1268, 1267]
            if lobe == 'inferior':
                m1053007 = [1983, 1923]
                m2053007 = [1273, 1281]
                m1053107 = [2058, 1074]
                m2053107 = [1555, 1494]
                m1060107 = [1653, 1687]
                m1060407 = [1084, 1098]
                m2060407 = [690, 691]
                m3060407 = [878, 879]
                m1060507 = [1212, 1213]
                m2060507 = [943, 946]
                m3060507 = [1003, 1005]
                m2060607 = [720, 769]
                m3060607 = [1207, 1205]
            if lobe == 'postcaval':
                m1053007 = [1903, 511]
                m2053007 = [1224, 1361]
                m1053107 = [2057, 1273]
                m2053107 = [1550, 1556]
                m1060107 = [1652, 1713]
                m1060407 = [1059, 1083]
                m2060407 = [783, 739]
                m3060407 = [865, 920]
                m1060507 = [1234, 1326]
                m2060507 = [945, 1036]
                m3060507 = [986, 985]
                m2060607 = [753, 752]
                m3060607 = [1207, 1274]
        if lobe == 'middle':
            datasets = {
            'm1053007': twoInputVessels('m1p'+pressure+'_053007',m1053007[0], m1053007[1], m1053007[2], m1053007[3], 1),
            'm2053007': twoInputVessels('m2p'+pressure+'_053007',m2053007[0], m2053007[1], m2053007[2], m2053007[3], 1),
            'm1053107': twoInputVessels('m1p'+pressure+'_053107',m1053107[0], m1053107[1], m1053107[2], m1053107[3], 1),
            'm2053107': twoInputVessels('m2p'+pressure+'_053107',m2053107[0], m2053107[1], m2053107[2], m2053107[3], 1),
            'm1060107': twoInputVessels('m1p'+pressure+'_060107',m1060107[0], m1060107[1], m1060107[2], m1060107[3], 1),
            'm1060407': twoInputVessels('m1p'+pressure+'_060407',m1060407[0], m1060407[1], m1060407[2], m1060407[3], 0),
            'm2060407': twoInputVessels('m2p'+pressure+'_060407',m2060407[0], m2060407[1], m2060407[2], m2060407[3], 0),
            'm3060407': twoInputVessels('m3p'+pressure+'_060407',m3060407[0], m3060407[1], m3060407[2], m3060407[3], 0),
            'm1060507': twoInputVessels('m1p'+pressure+'_060507',m1060507[0], m1060507[1], m1060507[2], m1060507[3], 0),
            'm2060507': twoInputVessels('m2p'+pressure+'_060507',m2060507[0], m2060507[1], m2060507[2], m2060507[3], 0),
            'm3060507': twoInputVessels('m3p'+pressure+'_060507',m3060507[0], m3060507[1], m3060507[2], m3060507[3], 0),
            'm2060607': twoInputVessels('m2p'+pressure+'_060607',m2060607[0], m2060607[1], m2060607[2], m2060607[3], 0),
            'm3060607': twoInputVessels('m3p'+pressure+'_060607',m3060607[0], m3060607[1], m3060607[2], m3060607[3], 0),
            }
        elif lobe == 'superior':
            if pressure == 1:
                datasets = {
                'm1053007': lobeTermLoc('m1p'+pressure+'_053007',m1053007[0], m1053007[1], 1),
                'm2053007': lobeTermLoc('m2p'+pressure+'_053007',m2053007[0], m2053007[1], 1),
                'm1053107': lobeTermLoc('m1p'+pressure+'_053107',m1053107[0], m1053107[1], 1),
                'm2053107': lobeTermLoc('m2p'+pressure+'_053107',m2053107[0], m2053107[1], 1),
                'm1060107': twoInputVessels('m1p'+pressure+'_060107',m1060107[0], m1060107[1], m1060107[2], m1060107[3], 1),
                'm1060407': lobeTermLoc('m1p'+pressure+'_060407',m1060407[0], m1060407[1], 0),
                'm2060407': lobeTermLoc('m2p'+pressure+'_060407',m2060407[0], m2060407[1], 0),
                'm3060407': twoInputVessels('m3p'+pressure+'_060407',m3060407[0], m3060407[1], m3060407[2], m3060407[3], 0),
                'm1060507': lobeTermLoc('m1p'+pressure+'_060507',m1060507[0], m1060507[1], 0),
                'm2060507': lobeTermLoc('m2p'+pressure+'_060507',m2060507[0], m2060507[1], 0),
                'm3060507': lobeTermLoc('m3p'+pressure+'_060507',m3060507[0], m3060507[1], 0),
                'm2060607': lobeTermLoc('m2p'+pressure+'_060607',m2060607[0], m2060607[1], 0),
                'm3060607': lobeTermLoc('m3p'+pressure+'_060607',m3060607[0], m3060607[1], 0),
                }
            else: 
                datasets = {
                'm1053007': lobeTermLoc('m1p'+pressure+'_053007',m1053007[0], m1053007[1], 1),
                'm2053007': lobeTermLoc('m2p'+pressure+'_053007',m2053007[0], m2053007[1], 1),
                'm1053107': lobeTermLoc('m1p'+pressure+'_053107',m1053107[0], m1053107[1], 1),
                'm2053107': lobeTermLoc('m2p'+pressure+'_053107',m2053107[0], m2053107[1], 1),
                'm1060107': lobeTermLoc('m1p'+pressure+'_060107',m1060107[0], m1060107[1], 1),
                'm1060407': lobeTermLoc('m1p'+pressure+'_060407',m1060407[0], m1060407[1], 0),
                'm2060407': lobeTermLoc('m2p'+pressure+'_060407',m2060407[0], m2060407[1], 0),
                'm3060407': twoInputVessels('m3p'+pressure+'_060407',m3060407[0], m3060407[1], m3060407[2], m3060407[3], 0),
                'm1060507': lobeTermLoc('m1p'+pressure+'_060507',m1060507[0], m1060507[1], 0),
                'm2060507': lobeTermLoc('m2p'+pressure+'_060507',m2060507[0], m2060507[1], 0),
                'm3060507': lobeTermLoc('m3p'+pressure+'_060507',m3060507[0], m3060507[1], 0),
                'm2060607': lobeTermLoc('m2p'+pressure+'_060607',m2060607[0], m2060607[1], 0),
                'm3060607': lobeTermLoc('m3p'+pressure+'_060607',m3060607[0], m3060607[1], 0),
                }
        else:
            datasets = {
            'm1053007': lobeTermLoc('m1p'+pressure+'_053007',m1053007[0], m1053007[1], 1),
            'm2053007': lobeTermLoc('m2p'+pressure+'_053007',m2053007[0], m2053007[1], 1),
            'm1053107': lobeTermLoc('m1p'+pressure+'_053107',m1053107[0], m1053107[1], 1),
            'm2053107': lobeTermLoc('m2p'+pressure+'_053107',m2053107[0], m2053107[1], 1),
            'm1060107': lobeTermLoc('m1p'+pressure+'_060107',m1060107[0], m1060107[1], 1),
            'm1060407': lobeTermLoc('m1p'+pressure+'_060407',m1060407[0], m1060407[1], 0),
            'm2060407': lobeTermLoc('m2p'+pressure+'_060407',m2060407[0], m2060407[1], 0),
            'm3060407': lobeTermLoc('m3p'+pressure+'_060407',m3060407[0], m3060407[1], 0),
            'm1060507': lobeTermLoc('m1p'+pressure+'_060507',m1060507[0], m1060507[1], 0),
            'm2060507': lobeTermLoc('m2p'+pressure+'_060507',m2060507[0], m2060507[1], 0),
            'm3060507': lobeTermLoc('m3p'+pressure+'_060507',m3060507[0], m3060507[1], 0),
            'm2060607': lobeTermLoc('m2p'+pressure+'_060607',m2060607[0], m2060607[1], 0),
            'm3060607': lobeTermLoc('m3p'+pressure+'_060607',m3060607[0], m3060607[1], 0),
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
metrics_df.to_csv('PairwiseDistances/Lobe/curve_metrics_summary.csv', index=False)
print(f"\nSaved curve metrics summary to PairwiseDistances/Lobe/curve_metrics_summary.csv")
print(f"Total rows: {len(metrics_df)}")
print(f"Columns: {list(metrics_df.columns)}")