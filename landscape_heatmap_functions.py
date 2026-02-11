"""
Functions to add to your run_landscapes_withentropyloss_whole.py script
These work directly with the landscape_curves dictionary to create distance heatmaps
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os

def compute_landscape_metrics_from_curves(t_vals, landscape):
    """
    Compute metrics from a single landscape curve.
    
    Parameters:
    -----------
    t_vals : array
        Time/filtration values
    landscape : array
        Single landscape curve values
        
    Returns:
    --------
    dict with 'peak_x', 'peak_y', 'auc'
    """
    if len(landscape) == 0 or len(t_vals) == 0:
        return {'peak_x': 0, 'peak_y': 0, 'auc': 0}
    
    # Find peak
    peak_idx = np.argmax(landscape)
    peak_x = t_vals[peak_idx]
    peak_y = landscape[peak_idx]
    
    # Compute AUC
    auc = np.trapezoid(landscape, t_vals)
    
    return {'peak_x': peak_x, 'peak_y': peak_y, 'auc': auc}

def assign_mouse_groups():
    """
    Define which mice belong to which group.
    
    Returns:
    --------
    dict mapping mouse names to group labels
    """
    hyper_mice = ['m1053007', 'm2053007', 'm1053107', 'm2053107', 'm1060107']
    control_mice = ['m1060407', 'm2060407', 'm3060407', 'm1060507', 'm2060507', 
                    'm3060507', 'm2060607', 'm3060607']
    
    groups = {}
    for mouse in hyper_mice:
        groups[mouse] = 'Hyper'
    for mouse in control_mice:
        groups[mouse] = 'Control'
    
    return groups

def compute_pairwise_distances_from_landscapes(landscape_curves, dim, landscape_num):
    """
    Compute pairwise L2 distances between all mice for a specific landscape.
    
    Parameters:
    -----------
    landscape_curves : dict
        Dictionary of {mouse_name: {dim: (t_vals, landscapes)}}
    dim : int
        Dimension to analyze
    landscape_num : int
        Which landscape to analyze (0-indexed)
        
    Returns:
    --------
    distances : array
        Pairwise distance matrix
    mouse_labels : array
        Mouse names in order
    groups : array
        Group assignments in order
    """
    mouse_groups = assign_mouse_groups()
    
    # Extract features for all mice
    features_dict = {}
    
    for mouse_name, curves in landscape_curves.items():
        if dim not in curves:
            continue
            
        t_vals, landscapes = curves[dim]
        
        # Check if this landscape number exists
        if landscape_num >= landscapes.shape[0]:
            continue
            
        # Get the specific landscape
        landscape = landscapes[landscape_num]
        
        # Compute metrics
        metrics = compute_landscape_metrics_from_curves(t_vals, landscape)
        
        # Store as feature vector
        features_dict[mouse_name] = np.array([
            metrics['peak_x'],
            metrics['peak_y'],
            metrics['auc']
        ])
    
    # Convert to arrays
    mouse_labels = np.array(list(features_dict.keys()))
    features = np.array(list(features_dict.values()))
    groups = np.array([mouse_groups.get(m, 'Unknown') for m in mouse_labels])
    
    # Compute pairwise distances
    n = len(features)
    distances = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            distances[i, j] = np.linalg.norm(features[i] - features[j])
    
    # Verify diagonal is zero
    diagonal = np.diag(distances)
    if not np.allclose(diagonal, 0):
        print(f"WARNING: Diagonal is not zero!")
        print(f"  Diagonal values: {diagonal}")
        print(f"  Features shape: {features.shape}")
        for i in range(min(3, n)):
            print(f"  Mouse {mouse_labels[i]}: features={features[i]}")
            print(f"    Self-distance: {distances[i,i]}")
    
    return distances, mouse_labels, groups

def normalize_distances(distances):
    """Normalize distance matrix to [0, 1] range"""
    if distances is None or distances.size == 0:
        return None
    max_dist = np.max(distances)
    if max_dist > 0:
        return distances / max_dist
    return distances

def compute_group_statistics(distances, mouse_labels, groups):
    """
    Compute within-group and between-group statistics.
    
    Returns:
    --------
    dict with statistics
    """
    hyper_idx = np.where(groups == 'Hyper')[0]
    control_idx = np.where(groups == 'Control')[0]
    
    stats = {}
    
    # Within Group 1 (Hyper)
    if len(hyper_idx) > 1:
        within_hyper = []
        for i in range(len(hyper_idx)):
            for j in range(i+1, len(hyper_idx)):
                within_hyper.append(distances[hyper_idx[i], hyper_idx[j]])
        if within_hyper:
            stats['within_hyper_mean'] = np.mean(within_hyper)
            stats['within_hyper_std'] = np.std(within_hyper)
            stats['within_hyper_n'] = len(within_hyper)
        else:
            stats['within_hyper_mean'] = np.nan
            stats['within_hyper_std'] = np.nan
            stats['within_hyper_n'] = 0
    else:
        stats['within_hyper_mean'] = np.nan
        stats['within_hyper_std'] = np.nan
        stats['within_hyper_n'] = 0
    
    # Within Group 2 (Control)
    if len(control_idx) > 1:
        within_control = []
        for i in range(len(control_idx)):
            for j in range(i+1, len(control_idx)):
                within_control.append(distances[control_idx[i], control_idx[j]])
        if within_control:
            stats['within_control_mean'] = np.mean(within_control)
            stats['within_control_std'] = np.std(within_control)
            stats['within_control_n'] = len(within_control)
        else:
            stats['within_control_mean'] = np.nan
            stats['within_control_std'] = np.nan
            stats['within_control_n'] = 0
    else:
        stats['within_control_mean'] = np.nan
        stats['within_control_std'] = np.nan
        stats['within_control_n'] = 0
    
    # Between groups
    if len(hyper_idx) > 0 and len(control_idx) > 0:
        between = []
        for i in hyper_idx:
            for j in control_idx:
                between.append(distances[i, j])
        stats['between_mean'] = np.mean(between)
        stats['between_std'] = np.std(between)
        stats['between_n'] = len(between)
    else:
        stats['between_mean'] = np.nan
        stats['between_std'] = np.nan
        stats['between_n'] = 0
    
    # Overall average
    upper_tri = []
    for i in range(len(distances)):
        for j in range(i+1, len(distances)):
            upper_tri.append(distances[i, j])
    stats['overall_mean'] = np.mean(upper_tri) if upper_tri else np.nan
    stats['overall_std'] = np.std(upper_tri) if upper_tri else np.nan
    
    return stats

def compute_distance_to_group_averages(distances, mouse_labels, groups):
    """
    Compute distance of each mouse to the average of each group.
    """
    hyper_idx = np.where(groups == 'Hyper')[0]
    control_idx = np.where(groups == 'Control')[0]
    
    n = len(mouse_labels)
    dist_to_hyper_avg = np.zeros(n)
    dist_to_control_avg = np.zeros(n)
    
    for i in range(n):
        if len(hyper_idx) > 0:
            # Average distance to all hyper mice (excluding self if applicable)
            hyper_dists = [distances[i, j] for j in hyper_idx if i != j]
            dist_to_hyper_avg[i] = np.mean(hyper_dists) if hyper_dists else 0
        else:
            dist_to_hyper_avg[i] = np.nan
            
        if len(control_idx) > 0:
            control_dists = [distances[i, j] for j in control_idx if i != j]
            dist_to_control_avg[i] = np.mean(control_dists) if control_dists else 0
        else:
            dist_to_control_avg[i] = np.nan
    
    return dist_to_hyper_avg, dist_to_control_avg

def plot_landscape_distance_heatmap(distances, mouse_labels, groups, stats, 
                                    dim, landscape_num, pressure, lobe, output_dir,
                                    normalize=True):
    """
    Create and save a comprehensive distance heatmap with BOTH group statistics AND average columns.
    
    This merged version includes:
    - Pairwise distances between all mice
    - Three additional columns: Hyper_Avg, Control_Avg, Overall_Avg
    - Statistics panel showing within-group and between-group metrics
    
    Parameters:
    -----------
    distances : array
        Pairwise distance matrix
    mouse_labels : array
        Mouse names
    groups : array
        Group assignments
    stats : dict
        Group statistics
    dim : int
        Dimension number
    landscape_num : int
        Landscape number (1-indexed for display)
    pressure : str or int
        Pressure level
    lobe : str or None
        Lobe identifier or None
    output_dir : str
        Directory to save the plot
    normalize : bool
        Whether to normalize distances
    """
    if distances is None or distances.size == 0:
        lobe_str = f", Lobe {lobe}" if lobe else ""
        print(f"Skipping Dim {dim}, Landscape {landscape_num}, Pressure {pressure}{lobe_str}: No data")
        return
    
    # VERIFY: Check that diagonal is zero
    n = len(mouse_labels)
    diagonal_orig = np.diag(distances)
    if not np.allclose(diagonal_orig, 0, atol=1e-10):
        print(f"ERROR: Input distances diagonal is not zero!")
        print(f"  Diagonal: {diagonal_orig}")
        return
    
    # Compute distance to group averages
    dist_to_hyper, dist_to_control = compute_distance_to_group_averages(
        distances, mouse_labels, groups
    )
    
    # Compute overall average distance for each mouse
    n = len(mouse_labels)
    overall_avg = np.array([np.mean([distances[i, j] for j in range(n) if i != j]) 
                           for i in range(n)])
    
    # Normalize BEFORE creating extended matrix to preserve structure
    if normalize:
        # Find max distance for normalization (excluding diagonal which is 0)
        max_dist = np.max(distances)
        
        if max_dist > 0:
            # Normalize pairwise distances
            normalized_distances = distances / max_dist
            # Normalize average columns using same scale
            normalized_hyper_avg = dist_to_hyper / max_dist
            normalized_control_avg = dist_to_control / max_dist
            normalized_overall_avg = overall_avg / max_dist
        else:
            normalized_distances = distances
            normalized_hyper_avg = dist_to_hyper
            normalized_control_avg = dist_to_control
            normalized_overall_avg = overall_avg
        
        # Create extended matrix from normalized components
        plot_distances = np.column_stack([
            normalized_distances,
            normalized_hyper_avg.reshape(-1, 1),
            normalized_control_avg.reshape(-1, 1),
            normalized_overall_avg.reshape(-1, 1)
        ])
        
        title_prefix = "Norm_"
        value_label = "Normalized L² Distance"
    else:
        # Create extended matrix from raw distances
        plot_distances = np.column_stack([
            distances,
            dist_to_hyper.reshape(-1, 1),
            dist_to_control.reshape(-1, 1),
            overall_avg.reshape(-1, 1)
        ])
        
        title_prefix = ""
        value_label = "L² Distance"
    
    # Sort by group for better visualization
    sort_idx = np.argsort(groups)[::-1]  # Hyper first, then Control
    sorted_labels = mouse_labels[sort_idx]
    sorted_groups = groups[sort_idx]
    
    # The extended matrix has structure: [n×n pairwise | n×3 averages]
    # We need to sort rows for everything, but columns only for the pairwise part
    pairwise_block = plot_distances[:, :n]  # First n columns (pairwise)
    averages_block = plot_distances[:, n:]  # Last 3 columns (averages)
    
    # Sort pairwise block: both rows AND columns
    pairwise_sorted = pairwise_block[sort_idx][:, sort_idx]
    
    # Sort averages block: only rows
    averages_sorted = averages_block[sort_idx]
    
    # Recombine
    plot_distances = np.column_stack([pairwise_sorted, averages_sorted])
    
    # VERIFY: Final check that diagonal is still zero in the pairwise section
    diagonal_final = np.array([plot_distances[i, i] for i in range(n)])
    if not np.allclose(diagonal_final, 0, atol=1e-10):
        print(f"ERROR: Final diagonal is not zero after sorting!")
        print(f"  Diagonal: {diagonal_final}")
        print(f"  Expected all zeros")
    
    # Create figure with more width for the extra columns
    fig, ax = plt.subplots(figsize=(20, 12))
    
    # Create heatmap
    im = ax.imshow(plot_distances, cmap='viridis', aspect='auto')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(value_label, fontsize=14)
    
    # X-axis labels (mice + 3 average columns)
    x_labels = list(sorted_labels) + ['Hyper_Avg', 'Control_Avg', 'Overall_Avg']
    ax.set_xticks(np.arange(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=10)
    
    # Y-axis labels
    ax.set_yticks(np.arange(len(sorted_labels)))
    ax.set_yticklabels(sorted_labels, fontsize=10)
    
    # Add values in cells
    for i in range(len(sorted_labels)):
        for j in range(len(x_labels)):
            value = plot_distances[i, j]
            text = ax.text(j, i, f'{value:.3f}',
                          ha="center", va="center", 
                          color="white" if value > 0.5 else "black",
                          fontsize=7)
    
    # Draw group separation lines (horizontal and vertical through mice)
    hyper_count = np.sum(sorted_groups == 'Hyper')
    if 0 < hyper_count < len(sorted_groups):
        ax.axhline(y=hyper_count - 0.5, color='red', linewidth=3, linestyle='--')
        ax.axvline(x=hyper_count - 0.5, color='red', linewidth=3, linestyle='--')
    
    # Draw vertical line before group average columns
    ax.axvline(x=n - 0.5, color='red', linewidth=2, linestyle=':')
    
    # Title (include lobe if present)
    lobe_str = f", {lobe}" if lobe else ""
    title = f'{title_prefix}Landscape{landscape_num} Curve Pairwise L² Distances (Dimension {dim}, Pressure {pressure}{lobe_str})'
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Curves and Group Averages', fontsize=12)
    ax.set_ylabel('Curves', fontsize=12)
    
    # Add statistics box
    textstr = 'Group Statistics\n' + '='*25 + '\n\n'
    
    if not np.isnan(stats['within_hyper_mean']):
        textstr += f'Within Group 1 (Hyper):\n'
        textstr += f'  Mean:  {stats["within_hyper_mean"]:.4f}\n'
        textstr += f'  Std:   {stats["within_hyper_std"]:.4f}\n'
        textstr += f'  N:     {stats["within_hyper_n"]}\n\n'
    
    if not np.isnan(stats['within_control_mean']):
        textstr += f'Within Group 2 (Control):\n'
        textstr += f'  Mean:  {stats["within_control_mean"]:.4f}\n'
        textstr += f'  Std:   {stats["within_control_std"]:.4f}\n'
        textstr += f'  N:     {stats["within_control_n"]}\n\n'
    
    if not np.isnan(stats['between_mean']):
        textstr += f'Between Groups:\n'
        textstr += f'  Mean:  {stats["between_mean"]:.4f}\n'
        textstr += f'  Std:   {stats["between_std"]:.4f}\n'
        textstr += f'  N:     {stats["between_n"]}\n'
    
    # Position the text box on the right side
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(1.12, 0.5, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='center', bbox=props, family='monospace')
    
    plt.tight_layout()
    
    # Save figure (include lobe in filename if present)
    lobe_suffix = f"_{lobe}" if lobe else ""
    filename = f'{title_prefix}Pressure{pressure}{lobe_suffix}_Dim{dim}_Landscape{landscape_num}.png'
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved: {filename}")
    return filepath

def generate_heatmaps_for_pressure_and_lobe(landscape_curves, pressure, lobe, global_k, output_dir):
    """
    Generate all heatmaps for a single pressure and lobe combination.
    
    Parameters:
    -----------
    landscape_curves : dict
        Dictionary of {mouse_name: {dim: (t_vals, landscapes)}}
    pressure : str or int
        Pressure level
    lobe : str or None
        Lobe identifier (e.g., 'lobe1', 'RUL', etc.) or None if not using lobes
    global_k : dict
        Dictionary of {dim: k_value} for this pressure/lobe
    output_dir : str
        Directory to save the plots
        
    Returns:
    --------
    list of generated file paths
    """
    os.makedirs(output_dir, exist_ok=True)
    
    generated_files = []
    
    lobe_str = f", Lobe {lobe}" if lobe else ""
    print(f"\n{'='*70}")
    print(f"Generating heatmaps for Pressure {pressure}{lobe_str}")
    print(f"  Mice available: {list(landscape_curves.keys())}")
    print(f"  global_k: {global_k}")
    print(f"{'='*70}")
    
    # For each dimension
    for dim in [0, 1, 2]:
        if dim not in global_k:
            print(f"  Dimension {dim} not in global_k, skipping")
            continue
            
        k = global_k[dim]
        print(f"  Dimension {dim}: Processing {k} landscapes")
        
        # For each landscape (0-indexed, but display as 1-indexed)
        for landscape_idx in range(k):
            landscape_num_display = landscape_idx + 1
            
            # Compute distances
            distances, mouse_labels, groups = compute_pairwise_distances_from_landscapes(
                landscape_curves, dim, landscape_idx
            )
            
            if distances is None or distances.size == 0:
                print(f"    Landscape {landscape_num_display}: No data (returned None or empty)")
                continue
            
            print(f"    Landscape {landscape_num_display}: Found {len(mouse_labels)} mice with data")
            
            # Compute statistics
            stats = compute_group_statistics(distances, mouse_labels, groups)
            
            # Generate merged heatmap (with both stats panel AND average columns)
            filepath = plot_landscape_distance_heatmap(
                distances, mouse_labels, groups, stats,
                dim, landscape_num_display, pressure, lobe, output_dir,
                normalize=True
            )
            if filepath:
                generated_files.append(filepath)
    
    print(f"Generated {len(generated_files)} heatmaps for Pressure {pressure}{lobe_str}")
    return generated_files


# ============================================================================
# MAIN FUNCTION TO CALL FROM YOUR SCRIPT
# ============================================================================

def generate_all_landscape_heatmaps(all_landscape_curves, global_k_by_pressure, 
                                   output_base_dir='landscape_heatmaps'):
    """
    Generate all heatmaps for all pressures and lobes.
    
    CALL THIS FUNCTION after your main loop completes.
    
    Each heatmap includes:
    - Pairwise distances between all mice
    - Three additional columns: Hyper_Avg, Control_Avg, Overall_Avg
    - Statistics panel with within-group and between-group metrics
    
    Parameters:
    -----------
    all_landscape_curves : dict
        Dictionary with keys like 'm1053007_p1_lobe1', 'm2053007_p1_lobe1', etc.
        OR 'm1053007_p1', 'm2053007_p1' (if no lobe info)
        Each value is {dim: (t_vals, landscapes)}
    global_k_by_pressure : dict
        Dictionary of {pressure: {dim: k_value}}
        OR {(pressure, lobe): {dim: k_value}} if using lobes
    output_base_dir : str
        Base directory for all outputs
        
    Usage in your script:
    ---------------------
    # After all your pressure/lobe loops complete:
    generate_all_landscape_heatmaps(all_landscape_curves, global_k_by_pressure)
    """
    print("\n" + "="*70)
    print("GENERATING ALL LANDSCAPE DISTANCE HEATMAPS")
    print("="*70)
    
    all_files = []
    
    # Group landscapes by (pressure, lobe)
    grouped_landscapes = {}
    
    for full_name, curves in all_landscape_curves.items():
        # Parse the full_name to extract pressure and optionally lobe
        # Format: 'mousename_p{pressure}_{lobe}' or 'mousename_p{pressure}'
        import re
        match = re.match(r'^(.+)_p(\d)_(.+)$', full_name)
        if match:
            # Has lobe: mousename_p{pressure}_{lobe}
            mouse_name = match.group(1)
            pressure = match.group(2)
            lobe = match.group(3)
        else:
            # Try without lobe: mousename_p{pressure}
            match = re.match(r'^(.+)_p(\d)$', full_name)
            if match:
                mouse_name = match.group(1)
                pressure = match.group(2)
                lobe = None
            else:
                print(f"Warning: Skipping {full_name} - cannot parse format")
                continue
        
        key = (pressure, lobe)
        
        # Initialize group if needed
        if key not in grouped_landscapes:
            grouped_landscapes[key] = {}
        
        # Add mouse to this group
        grouped_landscapes[key][mouse_name] = curves
    
    # Generate heatmaps for each (pressure, lobe) combination
    print(f"\nFound {len(grouped_landscapes)} unique (pressure, lobe) combinations:")
    for key in sorted(grouped_landscapes.keys()):
        print(f"  {key}: {len(grouped_landscapes[key])} mice")
    
    for (pressure, lobe), pressure_landscapes in grouped_landscapes.items():
        # Determine the correct global_k
        # Try (pressure, lobe) first, then just pressure
        if (pressure, lobe) in global_k_by_pressure:
            global_k = global_k_by_pressure[(pressure, lobe)]
        elif pressure in global_k_by_pressure:
            global_k = global_k_by_pressure[pressure]
        else:
            print(f"WARNING: No global_k found for Pressure {pressure}, Lobe {lobe}")
            print(f"  Available keys in global_k_by_pressure: {list(global_k_by_pressure.keys())}")
            continue
        
        if not pressure_landscapes:
            print(f"WARNING: No landscape data for Pressure {pressure}, Lobe {lobe}")
            continue
        
        # Generate heatmaps for this pressure/lobe combination
        # Create subdirectory structure: Pressure{X}/lobe{Y}/ or just Pressure{X}/
        if lobe is not None:
            output_subdir = f'{output_base_dir}/Complement/Pressure{pressure}/{lobe}'
        else:
            output_subdir = f'{output_base_dir}/Full/Pressure{pressure}'
        
        os.makedirs(output_subdir, exist_ok=True)
        
        print(f"\nProcessing Pressure {pressure}, Lobe {lobe if lobe else 'N/A'}")
        files = generate_heatmaps_for_pressure_and_lobe(
            pressure_landscapes,
            pressure,
            lobe,
            global_k,
            output_subdir
        )
        all_files.extend(files)
    
    print("\n" + "="*70)
    print(f"COMPLETE! Generated {len(all_files)} total heatmaps")
    print(f"Each heatmap includes pairwise distances + group average columns + statistics")
    print(f"Output directory: {output_base_dir}")
    print("="*70)
    
    return all_files