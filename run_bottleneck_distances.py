import gudhi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from gudhi.hera import bottleneck_distance
from itertools import combinations

# Global functions
def nodesExtractor(name):
    """Extracts nodes and their corresponding information"""
    file_path = 'fullComplement/' + name + '.csv'
    points = pd.read_csv(file_path)
    points = points.iloc[:, 1:]
    return points

def nodesToArray(name):
    """Converts nodes to numpy array of coordinates"""
    nodes_df = nodesExtractor(name)
    nodes_loc = nodes_df.loc[:, ['x', 'y', 'z']]
    loc_array = nodes_loc.to_numpy()
    return loc_array

def compute_persistence(points):
    """
    Compute persistence diagram and track the last death in each dimension.
    Excludes infinite values in H1 and H2.
    """
    # Build alpha complex
    alpha_complex = gudhi.AlphaComplex(points=points)
    simplex_tree = alpha_complex.create_simplex_tree()
    simplex_tree.compute_persistence()
    
    # Collect persistence pairs
    persistence_pairs = simplex_tree.persistence()
    
    # Check for infinite bars
    for dim, (birth, death) in persistence_pairs:
        if np.isinf(death):
            print(f"WARNING: Infinite bar in H{dim}: birth={birth}, death=∞")
    
    # Store persistence diagram
    diag = []
    for dim, (birth, death) in persistence_pairs:
        diag.append((dim, (birth, death)))
    
    # Track last death in each dimension
    last_deaths = {}
    for dim, (birth, death) in diag:
        if np.isfinite(death):
            if dim not in last_deaths:
                last_deaths[dim] = death
            else:
                last_deaths[dim] = max(last_deaths[dim], death)
    
    return diag, last_deaths

def extract_dimension_diagrams(persistence_diagram, dimension):
    """
    Extract persistence pairs for a specific dimension.
    Returns a numpy array of shape (n, 2) with birth-death pairs.
    Filters out infinite death values.
    """
    pairs = []
    for dim, (birth, death) in persistence_diagram:
        if dim == dimension and np.isfinite(death):
            pairs.append([birth, death])
    
    if len(pairs) == 0:
        return np.array([]).reshape(0, 2)
    
    return np.array(pairs)

def compute_bottleneck_distance_matrix(persistence_diagrams, dimension):
    """
    Compute pairwise bottleneck distances for a specific dimension.
    
    Parameters:
    -----------
    persistence_diagrams : dict
        Dictionary of {name: persistence_diagram}
    dimension : int
        Homology dimension (0, 1, or 2)
    
    Returns:
    --------
    distance_matrix : numpy.ndarray
        Symmetric matrix of pairwise bottleneck distances
    names : list
        List of dataset names corresponding to matrix rows/columns
    """
    names = list(persistence_diagrams.keys())
    n = len(names)
    distance_matrix = np.zeros((n, n))
    
    # Extract diagrams for the specific dimension
    dim_diagrams = {}
    for name, diag in persistence_diagrams.items():
        dim_diagrams[name] = extract_dimension_diagrams(diag, dimension)
    
    # Compute pairwise distances
    for i in range(n):
        for j in range(i + 1, n):
            name_i, name_j = names[i], names[j]
            diag_i = dim_diagrams[name_i]
            diag_j = dim_diagrams[name_j]
            
            # Handle empty diagrams
            if len(diag_i) == 0 and len(diag_j) == 0:
                dist = 0.0
            elif len(diag_i) == 0 or len(diag_j) == 0:
                # Distance is the largest persistence in the non-empty diagram divided by 2
                non_empty = diag_i if len(diag_i) > 0 else diag_j
                persistences = non_empty[:, 1] - non_empty[:, 0]
                dist = np.max(persistences) / 2.0 if len(persistences) > 0 else 0.0
            else:
                # Compute bottleneck distance
                print('Start: '+str(i)+', '+str(j))
                dist = bottleneck_distance(diag_i, diag_j, delta=1e-4)
                print('Finish: '+str(i)+', '+str(j))

            
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist
    
    return distance_matrix, names

def plot_bottleneck_heatmap(distance_matrix, names, dimension, pressure, filename=None):
    """
    Create and save a heatmap of bottleneck distances.
    
    Parameters:
    -----------
    distance_matrix : numpy.ndarray
        Pairwise distance matrix
    names : list
        List of dataset names
    dimension : int
        Homology dimension
    pressure : str
        Pressure level
    filename : str, optional
        Path to save the figure
    """
    plt.figure(figsize=(12, 10))
    
    # Create heatmap
    sns.heatmap(
        distance_matrix,
        xticklabels=names,
        yticklabels=names,
        cmap='viridis',
        annot=True,
        fmt='.3f',
        square=True,
        cbar_kws={'label': 'Bottleneck Distance'}
    )
    
    plt.title(f'Bottleneck Distance Heatmap - Dimension {dimension}, Pressure {pressure}', 
              fontsize=16, pad=20)
    plt.xlabel('Dataset', fontsize=12)
    plt.ylabel('Dataset', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Saved heatmap to {filename}")
    else:
        plt.show()
    
    plt.close()

def save_distance_matrix_to_csv(distance_matrix, names, dimension, pressure, filename):
    """
    Save distance matrix to CSV file.
    
    Parameters:
    -----------
    distance_matrix : numpy.ndarray
        Pairwise distance matrix
    names : list
        List of dataset names
    dimension : int
        Homology dimension
    pressure : str
        Pressure level
    filename : str
        Path to save CSV file
    """
    df = pd.DataFrame(distance_matrix, index=names, columns=names)
    df.to_csv(filename)
    print(f"Saved distance matrix to {filename}")

def create_summary_statistics(all_distance_matrices, output_filename='bottleneck_summary_stats.csv'):
    """
    Create a summary CSV with statistics for each dimension and pressure.
    
    Parameters:
    -----------
    all_distance_matrices : dict
        Dictionary of {(pressure, dimension): distance_matrix}
    output_filename : str
        Path to save the summary CSV
    """
    rows = []
    
    for (pressure, dimension), distance_matrix in all_distance_matrices.items():
        # Get upper triangle (excluding diagonal)
        upper_triangle = distance_matrix[np.triu_indices_from(distance_matrix, k=1)]
        
        row = {
            'pressure': pressure,
            'dimension': dimension,
            'mean_distance': np.mean(upper_triangle),
            'std_distance': np.std(upper_triangle),
            'min_distance': np.min(upper_triangle),
            'max_distance': np.max(upper_triangle),
            'median_distance': np.median(upper_triangle),
            'q25_distance': np.percentile(upper_triangle, 25),
            'q75_distance': np.percentile(upper_triangle, 75)
        }
        rows.append(row)
    
    df = pd.DataFrame(rows)
    df = df.sort_values(['dimension', 'pressure'])
    df.to_csv(output_filename, index=False)
    print(f"\nSaved summary statistics to {output_filename}")
    
    return df

# Main execution
print("="*80)
print("COMPUTING PAIRWISE BOTTLENECK DISTANCES")
print("="*80)

# Dictionary to store all distance matrices
all_distance_matrices = {}
all_names_by_pressure = {}

for pressure in ['1', '2', '3', '4']:
    print(f"\n{'='*80}")
    print(f"Processing Pressure {pressure}")
    print(f"{'='*80}")
    
    # Load datasets
    datasets = {
        'm1053007': nodesToArray('m1p' + pressure + '_053007_Full'),
        'm2053007': nodesToArray('m2p' + pressure + '_053007_Full'),
        'm1053107': nodesToArray('m1p' + pressure + '_053107_Full'),
        'm2053107': nodesToArray('m2p' + pressure + '_053107_Full'),
        'm1060107': nodesToArray('m1p' + pressure + '_060107_Full'),
        'm1060407': nodesToArray('m1p' + pressure + '_060407_Full'),
        'm2060407': nodesToArray('m2p' + pressure + '_060407_Full'),
        'm3060407': nodesToArray('m3p' + pressure + '_060407_Full'),
        'm1060507': nodesToArray('m1p' + pressure + '_060507_Full'),
        'm2060507': nodesToArray('m2p' + pressure + '_060507_Full'),
        'm3060507': nodesToArray('m3p' + pressure + '_060507_Full'),
        'm2060607': nodesToArray('m2p' + pressure + '_060607_Full'),
        'm3060607': nodesToArray('m3p' + pressure + '_060607_Full')
    }
    
    # Compute persistence diagrams
    print(f"\nComputing persistence diagrams...")
    persistence_results = {
        name: compute_persistence(points)
        for name, points in datasets.items()
    }
    
    # Separate results
    persistence_diagrams = {name: res[0] for name, res in persistence_results.items()}
    last_deaths = {name: res[1] for name, res in persistence_results.items()}
    
    # Store names for this pressure
    all_names_by_pressure[pressure] = list(datasets.keys())
    
    # Compute bottleneck distances for each dimension
    # uses 1e-4 relative error threshold
    for dimension in [0, 1, 2]:
        print(f"\n  Computing bottleneck distances for dimension {dimension}...")
        
        distance_matrix, names = compute_bottleneck_distance_matrix(
            persistence_diagrams, 
            dimension
        )
        
        # Store the distance matrix
        all_distance_matrices[(pressure, dimension)] = distance_matrix
        
        # Create output directory if it doesn't exist
        import os
        os.makedirs(f'bottleneck_results/pressure_{pressure}', exist_ok=True)
        
        # Save heatmap
        heatmap_filename = f'bottleneck_results/pressure_{pressure}/heatmap_dim{dimension}_p{pressure}_fullComplement.png'
        plot_bottleneck_heatmap(
            distance_matrix, 
            names, 
            dimension, 
            pressure, 
            filename=heatmap_filename
        )
        
        # Save distance matrix to CSV
        csv_filename = f'bottleneck_results/pressure_{pressure}/distances_dim{dimension}_p{pressure}_fullComplement.csv'
        save_distance_matrix_to_csv(
            distance_matrix, 
            names, 
            dimension, 
            pressure, 
            filename=csv_filename
        )
    
    print(f"\nCompleted Pressure {pressure}")

# Create summary statistics across all pressures and dimensions
print(f"\n{'='*80}")
print("Creating Summary Statistics")
print(f"{'='*80}")

summary_df = create_summary_statistics(
    all_distance_matrices,
    output_filename='bottleneck_results/bottleneck_summary_stats.csv'
)

print("\n" + "="*80)
print("SUMMARY STATISTICS")
print("="*80)
print(summary_df.to_string(index=False))

print("\n" + "="*80)
print("PROCESSING COMPLETE")
print("="*80)
print(f"\nResults saved in 'bottleneck_results/' directory:")
print(f"  - Heatmaps: bottleneck_results/pressure_X/heatmap_dimY_pX.png")
print(f"  - Distance matrices: bottleneck_results/pressure_X/distances_dimY_pX.csv")
print(f"  - Summary statistics: bottleneck_results/bottleneck_summary_stats.csv")
