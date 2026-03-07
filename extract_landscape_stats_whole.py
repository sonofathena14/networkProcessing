import gudhi
import numpy as np
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt
from gudhi.representations import Landscape


def nodesExtractorC(name):
    file_path = 'Networks/Network_Vessels_' + name + '.mat'
    matlab_data = scipy.io.loadmat(file_path)
    data_structure = matlab_data['nodesC2']
    nodes_data = data_structure.squeeze()
    nodes_df = pd.DataFrame(nodes_data, columns=['NodeID', 'X', 'Y', 'Z', 'Degree'])
    return nodes_df

def nodesExtractorH(name):
    file_path = 'Pruned/Pruned_Network_' + name + '.mat'
    matlab_data = scipy.io.loadmat(file_path)
    data_structure = matlab_data['nodesC3']
    nodes_data = data_structure.squeeze()
    nodes_df = pd.DataFrame(nodes_data, columns=['NodeID', 'X', 'Y', 'Z', 'Degree'])
    return nodes_df

def nodesToArrayC(name):
    nodes_df = nodesExtractorC(name)
    nodes_loc = nodes_df.loc[:, ['X', 'Y', 'Z']] / 1000
    return nodes_loc.to_numpy()

def nodesToArrayH(name):
    nodes_df = nodesExtractorH(name)
    nodes_loc = nodes_df.loc[:, ['X', 'Y', 'Z']] / 1000
    return nodes_loc.to_numpy()


def compute_persistence(points):
    """
    Compute persistence diagram and track the last finite death in each dimension.
    """
    alpha_complex = gudhi.AlphaComplex(points=points)
    simplex_tree = alpha_complex.create_simplex_tree()
    simplex_tree.compute_persistence()
    persistence_pairs = simplex_tree.persistence()

    for dim, (birth, death) in persistence_pairs:
        if np.isinf(death):
            print(f"WARNING: Infinite bar in H{dim}: birth={birth}, death=∞")

    diag = list(persistence_pairs)

    last_deaths = {}
    for dim, (birth, death) in diag:
        if np.isfinite(death):
            last_deaths[dim] = max(last_deaths.get(dim, death), death)

    return diag, last_deaths


def compute_landscape_curves(persistence_diagram, t_min, t_max_per_dim, max_dim=2, resolution=500):
    """
    Compute the first persistence landscape (lambda_1) for each dimension.
    t_max_per_dim: dict {dim: t_max} so each dimension is sampled over its own filtration range,
    preventing low-death dimensions (e.g. H0) from being compressed into a window sized for H2.

    Returns:
        dict of {dim: (t_vals, landscapes)} where landscapes has shape (1, resolution)
    """
    landscape_curves = {}

    for dim in range(max_dim + 1):
        t_max = t_max_per_dim.get(dim, max(t_max_per_dim.values()))
        t_vals = np.linspace(t_min, t_max, resolution)

        diagram = np.array([
            (birth, death)
            for d, (birth, death) in persistence_diagram
            if d == dim and np.isfinite(death)
        ])

        if diagram.shape[0] == 0:
            landscape_curves[dim] = (t_vals, np.zeros((1, resolution)))
            continue

        landscape = Landscape(num_landscapes=1, resolution=resolution, sample_range=(t_min, t_max))
        lam1 = landscape.fit_transform([diagram])[0].reshape(1, resolution)
        landscape_curves[dim] = (t_vals, lam1)

    return landscape_curves



def save_landscape(t_vals, landscapes, dim, filename, title=None, limit=[0, 100], dpi=300):
    """
    Save persistence landscape plot to file.
    landscapes: shape (k, resolution)
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


def create_landscape_summary_csv(
    landscape_curves_dict,
    persistence_diagrams_by_dataset,
    output_filename='landscape_metrics_summary.csv'
):
    """
    Create a CSV with one row per (mouse, pressure) containing 4 landscape metrics
    for each of the 3 homology dimensions (H0, H1, H2).

    Metrics per dimension:
      1. landscape1_max          : Maximum value of lambda_1
      2. total_persistence       : Sum of (death - birth) over all finite persistence pairs
      3. landscape1_L1_norm      : Area under lambda_1 (trapezoidal integral)
      4. landscape1_active_range : Filtration span where lambda_1 > 0
                                   (last nonzero t - first nonzero t)
    """
    dims = [0, 1, 2]

    columns = ['mouse_name', 'pressure']
    for dim in dims:
        prefix = f'H{dim}'
        columns += [
            f'{prefix}_landscape1_max',
            f'{prefix}_total_persistence',
            f'{prefix}_landscape1_L1_norm',
            f'{prefix}_landscape1_active_range',
        ]

    rows = []

    for dataset_name, landscape_data in landscape_curves_dict.items():
        if '_p' in dataset_name:
            mouse_name, pressure = dataset_name.rsplit('_p', 1)
        else:
            mouse_name = dataset_name
            pressure = ''

        row = {'mouse_name': mouse_name, 'pressure': pressure}
        pers_diag = persistence_diagrams_by_dataset.get(dataset_name, [])

        for dim in dims:
            prefix = f'H{dim}'

            # Total persistence from raw diagram
            finite_pairs = [
                (birth, death)
                for d, (birth, death) in pers_diag
                if d == dim and np.isfinite(death)
            ]
            row[f'{prefix}_total_persistence'] = (
                sum(d - b for b, d in finite_pairs) if finite_pairs else np.nan
            )

            # Landscape-based metrics
            if dim in landscape_data:
                t_vals, landscapes = landscape_data[dim]
                lam1 = landscapes[0]

                row[f'{prefix}_landscape1_max'] = float(np.max(lam1))
                row[f'{prefix}_landscape1_L1_norm'] = float(np.trapz(np.abs(lam1), t_vals))

                nonzero_idx = np.where(lam1 > 0)[0]
                if len(nonzero_idx) > 0:
                    active_range = float(t_vals[nonzero_idx[-1]] - t_vals[nonzero_idx[0]])
                else:
                    active_range = np.nan
                row[f'{prefix}_landscape1_active_range'] = active_range
            else:
                row[f'{prefix}_landscape1_max'] = np.nan
                row[f'{prefix}_landscape1_L1_norm'] = np.nan
                row[f'{prefix}_landscape1_active_range'] = np.nan

        rows.append(row)

    df = pd.DataFrame(rows, columns=columns)
    df = df.sort_values(['mouse_name', 'pressure']).reset_index(drop=True)
    df.to_csv(output_filename, index=False)
    print(f"\nSaved landscape metrics summary to {output_filename}")
    return df


# ── Main ──────────────────────────────────────────────────────────────────────

all_landscape_curves = {}
all_persistence_diagrams = {}

for pressure in ['1', '2', '3', '4']:
    datasets = {
        'm1053007': nodesToArrayH('m1p' + pressure + '_053007'),
        'm2053007': nodesToArrayH('m2p' + pressure + '_053007'),
        'm1053107': nodesToArrayH('m1p' + pressure + '_053107'),
        'm2053107': nodesToArrayH('m2p' + pressure + '_053107'),
        'm1060107': nodesToArrayH('m1p' + pressure + '_060107'),
        'm1060407': nodesToArrayC('m1p' + pressure + '_060407'),
        'm2060407': nodesToArrayC('m2p' + pressure + '_060407'),
        'm3060407': nodesToArrayC('m3p' + pressure + '_060407'),
        'm1060507': nodesToArrayC('m1p' + pressure + '_060507'),
        'm2060507': nodesToArrayC('m2p' + pressure + '_060507'),
        'm3060507': nodesToArrayC('m3p' + pressure + '_060507'),
        'm2060607': nodesToArrayC('m2p' + pressure + '_060607'),
        'm3060607': nodesToArrayC('m3p' + pressure + '_060607'),
    }

    persistence_results = {name: compute_persistence(points) for name, points in datasets.items()}
    persistence_diagrams = {name: res[0] for name, res in persistence_results.items()}
    last_deaths = {name: res[1] for name, res in persistence_results.items()}

    global_last_deaths = {}
    for deaths_by_dim in last_deaths.values():
        for dim, death in deaths_by_dim.items():
            global_last_deaths[dim] = max(global_last_deaths.get(dim, death), death)

    landscape_curves = {
        name: compute_landscape_curves(diag, t_min=0, t_max_per_dim=global_last_deaths, max_dim=2, resolution=100000)
        for name, diag in persistence_diagrams.items()
    }

    """# Save landscape plots
    for i, mouse_name in enumerate(landscape_curves.keys()):
        group = 'Hyper' if i < 5 else 'Control'
        for dim in [0, 1, 2]:
            filename = (
                f'fullLungGraphs/Pressure{pressure}/'
                f'Persistence_Landscape_{dim}_{group}_{mouse_name}_P{pressure}_fullLung'
            )
            t_vals, landscapes = landscape_curves[mouse_name][dim]
            save_landscape(
                t_vals, landscapes, dim=dim, filename=filename,
                title=f'Persistence Landscape {mouse_name} Dimension {dim} Full Lung',
                limit=[0, global_last_deaths[dim]], dpi=150
            )
    """
    # Accumulate results
    for name, curves in landscape_curves.items():
        full_name = name + '_p' + pressure
        all_landscape_curves[full_name] = curves
        all_persistence_diagrams[full_name] = persistence_diagrams[name]

    print("Completed Pressure " + pressure)

# Create summary CSV
print("\n=== Creating Landscape Metrics Summary CSV ===")
summary_df = create_landscape_summary_csv(
    all_landscape_curves,
    all_persistence_diagrams,
    output_filename='landscape_metrics_summary_full.csv'
)

print(f"\nTotal datasets processed: {len(all_landscape_curves)}")
print(f"Rows in summary: {len(summary_df)}")
print(summary_df.head())