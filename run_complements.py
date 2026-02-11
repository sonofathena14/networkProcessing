import gudhi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from gudhi.representations import Landscape
from scipy.interpolate import interp1d
import os

# Global functions
def nodesExtractor(name):
    """Extracts nodes and their corresponding information"""
    file_path = 'Complement/' + name + '.csv'
    points = pd.read_csv(file_path)
    points = points.iloc[:, 1:]
    return points

def nodesToArray(name):
    """Convert nodes to numpy array"""
    nodes_df = nodesExtractor(name)
    nodes_loc = nodes_df.loc[:, ['x', 'y', 'z']]
    loc_array = nodes_loc.to_numpy()
    return loc_array

def compute_persistence(points):
    """
    Compute persistence diagram and track the last death in each dimension.
    Excludes infinite values in H1 and H2.
    """
    alpha_complex = gudhi.AlphaComplex(points=points)
    simplex_tree = alpha_complex.create_simplex_tree()
    simplex_tree.compute_persistence()
    
    persistence_pairs = simplex_tree.persistence()
    
    # Check for infinite bars
    for dim, (birth, death) in persistence_pairs:
        if np.isinf(death):
            print(f"WARNING: Infinite bar in H{dim}: birth={birth}, death=∞")
    
    diag = [(dim, (birth, death)) for dim, (birth, death) in persistence_pairs]
    
    # Track last death in each dimension
    last_deaths = {}
    for dim, (birth, death) in diag:
        if np.isfinite(death):
            if dim not in last_deaths:
                last_deaths[dim] = death
            else:
                last_deaths[dim] = max(last_deaths[dim], death)
    
    return diag, last_deaths

def average_curves(curve_dicts, resolution=100000):
    """Generic averaging function for any type of curve dicts"""
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

def save_curves(curves_dicts, labels, dimension, filename, title='Curves Comparison', limit=[0, 100], dpi=300):
    """Save curves to file"""
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
            print(f"Warning: {label} does not contain Dimension-{dimension}")
    
    plt.title(f'{title} (Dimension-{dimension})', fontsize=28)
    plt.xlabel('Filtration Value', fontsize=20)
    plt.ylabel(f'Value', fontsize=20)
    plt.xlim(limit)
    plt.legend(fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved plot to {filename}")

def compute_betti_curve(diag, max_dim=2, resolution=100000, cutoff=None):
    """Compute Betti curves"""
    betti_curves = {}
    
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
        
        grid = np.linspace(global_min, global_max, resolution)
        curve = np.zeros_like(grid)
        
        for birth, death in diag_dim:
            if np.isinf(death):
                curve += (grid >= birth)
            else:
                curve += (grid >= birth) & (grid <= death)
        
        betti_curves[dim] = (grid, curve)
    
    return betti_curves

def extract_intervals(persistence_diagram, dim):
    """Extract finite intervals for a given dimension"""
    return np.array([
        (b, d)
        for d_dim, (b, d) in persistence_diagram
        if d_dim == dim and np.isfinite(d)
    ])

def compute_lifespan_curves(persistence_diagram, t_min, t_max, max_dim=2, resolution=100000):
    """Compute lifespan curves"""
    lifespan_curves = {}
    
    for dim in range(max_dim + 1):
        intervals = extract_intervals(persistence_diagram, dim)
        
        if intervals.size == 0:
            lifespan_curves[dim] = (np.zeros(resolution), np.zeros(resolution))
            continue
        
        births = intervals[:, 0]
        deaths = intervals[:, 1]
        lifespans = deaths - births
        
        t_vals = np.linspace(t_min, t_max, resolution)
        LC = np.zeros_like(t_vals)
        
        for i, t in enumerate(t_vals):
            alive = (births <= t) & (t < deaths)
            LC[i] = lifespans[alive].sum()
        
        lifespan_curves[dim] = (t_vals, LC)
    
    return lifespan_curves

def compute_norm_lifespan_curves(persistence_diagram, t_min, t_max, max_dim=2, resolution=100000):
    """Compute normalized lifespan curves"""
    lifespan_curves = {}
    
    for dim in range(max_dim + 1):
        intervals = extract_intervals(persistence_diagram, dim)
        
        if intervals.size == 0:
            lifespan_curves[dim] = (np.zeros(resolution), np.zeros(resolution))
            continue
        
        births = intervals[:, 0]
        deaths = intervals[:, 1]
        lifespans = deaths - births
        total_lifespan = lifespans.sum()
        
        t_vals = np.linspace(t_min, t_max, resolution)
        LC = np.zeros_like(t_vals)
        
        for i, t in enumerate(t_vals):
            alive = (births <= t) & (t < deaths)
            if total_lifespan > 0:
                LC[i] = lifespans[alive].sum() / total_lifespan
            else:
                LC[i] = 0
        
        lifespan_curves[dim] = (t_vals, LC)
    
    return lifespan_curves

def compute_landscape_curves(persistence_diagram, t_min, t_max, max_dim=2, k=3, resolution=500):
    """Compute persistence landscapes"""
    landscape_curves = {}
    
    for dim in range(max_dim + 1):
        diagram = np.array([
            (birth, death)
            for d, (birth, death) in persistence_diagram
            if d == dim and np.isfinite(death)
        ])
        
        if diagram.shape[0] == 0:
            t_vals = np.linspace(t_min, t_max, resolution)
            landscapes = np.zeros((k, resolution))
            landscape_curves[dim] = (t_vals, landscapes)
            continue
        
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
    """Save persistence landscape plot"""
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

def process_pressure_lobe(pressure, lobe):
    """Process analysis for a given pressure and lobe"""
    print(f"\n{'='*60}")
    print(f"Processing Pressure {pressure}, Lobe {lobe}")
    print(f"{'='*60}\n")
    
    # Create output directory
    output_dir = f'2DComplementGraphs/Pressure{pressure}/{lobe}'
    os.makedirs(output_dir, exist_ok=True)
    
    # Load datasets
    datasets = {
        'm1053007': nodesToArray(f'm1p{pressure}_053007_{lobe}'),
        'm2053007': nodesToArray(f'm2p{pressure}_053007_{lobe}'),
        'm1053107': nodesToArray(f'm1p{pressure}_053107_{lobe}'),
        'm2053107': nodesToArray(f'm2p{pressure}_053107_{lobe}'),
        'm1060107': nodesToArray(f'm1p{pressure}_060107_{lobe}'),
        'm1060407': nodesToArray(f'm1p{pressure}_060407_{lobe}'),
        'm2060407': nodesToArray(f'm2p{pressure}_060407_{lobe}'),
        'm3060407': nodesToArray(f'm3p{pressure}_060407_{lobe}'),
        'm1060507': nodesToArray(f'm1p{pressure}_060507_{lobe}'),
        'm2060507': nodesToArray(f'm2p{pressure}_060507_{lobe}'),
        'm3060507': nodesToArray(f'm3p{pressure}_060507_{lobe}'),
        'm2060607': nodesToArray(f'm2p{pressure}_060607_{lobe}'),
        'm3060607': nodesToArray(f'm3p{pressure}_060607_{lobe}')
    }
    
    # Compute persistence
    print("Computing persistence diagrams...")
    persistence_results = {
        name: compute_persistence(points)
        for name, points in datasets.items()
    }
    
    persistence_diagrams = {name: res[0] for name, res in persistence_results.items()}
    last_deaths = {name: res[1] for name, res in persistence_results.items()}
    
    # Compute global max deaths
    global_last_deaths = {}
    for name, deaths_by_dim in last_deaths.items():
        for dim, death in deaths_by_dim.items():
            if dim not in global_last_deaths:
                global_last_deaths[dim] = death
            else:
                global_last_deaths[dim] = max(global_last_deaths[dim], death)
    
    print("Global max deaths by dimension:")
    for dim in sorted(global_last_deaths.keys()):
        print(f"  H{dim}: {global_last_deaths[dim]:.4f}")
    
    # Compute all curves
    print("\nComputing Betti curves...")
    betti_curves = {}
    for name, diag in persistence_diagrams.items():
        betti_curves[name] = compute_betti_curve(diag, cutoff=max(global_last_deaths.values()))
    
    print("Computing lifespan curves...")
    lifespan_curves = {}
    for name, diag in persistence_diagrams.items():
        lifespan_curves[name] = compute_lifespan_curves(diag, t_min=0, t_max=max(global_last_deaths.values()))
    
    print("Computing normalized lifespan curves...")
    norm_lifespan_curves = {}
    for name, diag in persistence_diagrams.items():
        norm_lifespan_curves[name] = compute_norm_lifespan_curves(diag, t_min=0, t_max=max(global_last_deaths.values()))
    
    print("Computing persistence landscapes...")
    landscape_curves = {}
    for name, diag in persistence_diagrams.items():
        landscape_curves[name] = compute_landscape_curves(
            diag, t_min=0, t_max=max(global_last_deaths.values()), k=3, resolution=500
        )
    
    # Separate hyper and control groups
    keys = list(betti_curves.keys())
    hyper_betti = [betti_curves[k] for k in keys[:5]]
    control_betti = [betti_curves[k] for k in keys[5:]]
    
    hyper_lifespan = [lifespan_curves[k] for k in keys[:5]]
    control_lifespan = [lifespan_curves[k] for k in keys[5:]]
    
    hyper_norm_lifespan = [norm_lifespan_curves[k] for k in keys[:5]]
    control_norm_lifespan = [norm_lifespan_curves[k] for k in keys[5:]]
    
    # Compute averages
    print("\nComputing average curves...")
    avg_betti_hyper = average_curves(hyper_betti)
    avg_betti_control = average_curves(control_betti)
    
    avg_lifespan_hyper = average_curves(hyper_lifespan)
    avg_lifespan_control = average_curves(control_lifespan)
    
    avg_norm_lifespan_hyper = average_curves(hyper_norm_lifespan)
    avg_norm_lifespan_control = average_curves(control_norm_lifespan)
    
    # Save average curves
    print("\nSaving average curves...")
    for dim in range(3):
        # Betti curves
        filename = f'{output_dir}/Betti_{dim}_Average_P{pressure}_{lobe}'
        save_curves([avg_betti_hyper, avg_betti_control], ['Hyper', 'Control'], dim, filename,
                   f"Comparison of Average Betti Curves Lobe {lobe}", [0, global_last_deaths[dim]])
        
        # Lifespan curves
        filename = f'{output_dir}/Lifespan_{dim}_Average_P{pressure}_{lobe}'
        save_curves([avg_lifespan_hyper, avg_lifespan_control], ['Hyper', 'Control'], dim, filename,
                   f"Comparison of Average Lifespan Curves Lobe {lobe}", [0, global_last_deaths[dim]])
        
        # Normalized lifespan curves
        filename = f'{output_dir}/Norm_Lifespan_{dim}_Average_P{pressure}_{lobe}'
        save_curves([avg_norm_lifespan_hyper, avg_norm_lifespan_control], ['Hyper', 'Control'], dim, filename,
                   f"Comparison of Average Norm Lifespan Curves Lobe {lobe}", [0, global_last_deaths[dim]])
    
    # Save individual curves
    print("\nSaving individual Betti curves...")
    for dim in range(3):
        for i, mouse_name in enumerate(keys):
            group = 'Hyper' if i < 5 else 'Control'
            filename = f'{output_dir}/Betti_{dim}_{group}_{mouse_name}_P{pressure}_{lobe}'
            save_curves([betti_curves[mouse_name]], [mouse_name], dim, filename,
                       f"Betti Curve Lobe {lobe}", [0, global_last_deaths[dim]])
    
    print("\nSaving individual lifespan curves...")
    for dim in range(3):
        for i, mouse_name in enumerate(keys):
            group = 'Hyper' if i < 5 else 'Control'
            filename = f'{output_dir}/Lifespan_{dim}_{group}_{mouse_name}_P{pressure}_{lobe}'
            save_curves([lifespan_curves[mouse_name]], [mouse_name], dim, filename,
                       f"Lifespan Curve Lobe {lobe}", [0, global_last_deaths[dim]])
    
    print("\nSaving individual normalized lifespan curves...")
    for dim in range(3):
        for i, mouse_name in enumerate(keys):
            group = 'Hyper' if i < 5 else 'Control'
            filename = f'{output_dir}/Norm_Lifespan_{dim}_{group}_{mouse_name}_P{pressure}_{lobe}'
            save_curves([norm_lifespan_curves[mouse_name]], [mouse_name], dim, filename,
                       f"Norm Lifespan Curve Lobe {lobe}", [0, global_last_deaths[dim]])
    
    print("\nSaving persistence landscapes...")
    for i, mouse_name in enumerate(keys):
        group = 'Hyper' if i < 5 else 'Control'
        for dim in range(3):
            filename = f'{output_dir}/Persistence_Landscape_{dim}_{group}_{mouse_name}_P{pressure}_{lobe}'
            t_vals, landscapes = landscape_curves[mouse_name][dim]
            save_landscape(t_vals, landscapes, dim, filename,
                         f'Persistence Landscape {mouse_name} Dimension {dim} Lobe {lobe}',
                         [0, global_last_deaths[dim]], dpi=150)
    
    print(f"\nCompleted processing for Pressure {pressure}, Lobe {lobe}")

# Main execution
if __name__ == "__main__":
    pressures = ['1', '2', '3', '4']  # Add or modify pressure values as needed
    lobes = ['left', 'superior', 'middle', 'inferior', 'postcaval']  # Add or modify lobes as needed
    
    for pressure in pressures:
        for lobe in lobes:
            try:
                process_pressure_lobe(pressure, lobe)
            except Exception as e:
                print(f"\nERROR processing Pressure {pressure}, Lobe {lobe}:")
                print(f"  {str(e)}")
                print("Continuing with next combination...\n")
    
    print("\n" + "="*60)
    print("All processing complete!")
    print("="*60)