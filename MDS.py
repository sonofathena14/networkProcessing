import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
import seaborn as sns

def load_distance_matrix(csv_file):
    """
    Load distance matrix from CSV file.
    
    Parameters:
    -----------
    csv_file : str
        Path to CSV file containing distance matrix
    
    Returns:
    --------
    distance_matrix : numpy.ndarray
        Distance matrix as numpy array
    labels : list
        Row/column labels from the CSV
    """
    df = pd.read_csv(csv_file, index_col=0)
    labels = df.index.tolist()
    distance_matrix = df.values
    
    return distance_matrix, labels

def perform_mds(distance_matrix, n_components=2, random_state=42):
    """
    Perform Multidimensional Scaling on distance matrix.
    
    Parameters:
    -----------
    distance_matrix : numpy.ndarray
        Symmetric distance matrix
    n_components : int
        Number of dimensions for the embedding (default: 2)
    random_state : int
        Random seed for reproducibility
    
    Returns:
    --------
    embedding : numpy.ndarray
        MDS coordinates in n_components dimensions
    stress : float
        Final stress value (lower is better)
    """
    mds = MDS(
        n_components=n_components,
        dissimilarity='precomputed',
        random_state=random_state,
        max_iter=300,
        n_init=4
    )
    
    embedding = mds.fit_transform(distance_matrix)
    stress = mds.stress_
    
    return embedding, stress

def plot_mds_2d(embedding, labels, title='MDS Visualization', 
                first_n_blue=5, filename=None, show_labels=True):
    """
    Plot 2D MDS embedding with color coding.
    
    Parameters:
    -----------
    embedding : numpy.ndarray
        2D coordinates from MDS
    labels : list
        Labels for each point
    title : str
        Plot title
    first_n_blue : int
        Number of first samples to color blue (rest will be red)
    filename : str, optional
        Path to save the figure
    show_labels : bool
        Whether to show labels on points
    """
    plt.figure(figsize=(12, 10))
    
    # Create color array: first n are blue, rest are red
    colors = ['blue' if i < first_n_blue else 'red' for i in range(len(labels))]
    
    # Create scatter plot
    for i, (x, y) in enumerate(embedding):
        plt.scatter(x, y, c=colors[i], s=200, alpha=0.7, edgecolors='black', linewidth=1.5)
        
        if show_labels:
            # Add label next to each point
            plt.annotate(
                labels[i], 
                (x, y), 
                xytext=(5, 5), 
                textcoords='offset points',
                fontsize=10,
                fontweight='bold'
            )
    
    # Create custom legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='blue', edgecolor='black', label='Hypertensive'),
        Patch(facecolor='red', edgecolor='black', label='Control')
    ]
    plt.legend(handles=legend_elements, loc='best', fontsize=12)
    
    plt.xlabel('MDS Dimension 1', fontsize=14)
    plt.ylabel('MDS Dimension 2', fontsize=14)
    plt.title(title, fontsize=16, pad=20)
    plt.grid(True, alpha=0.3)
    plt.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
    plt.axvline(x=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Saved plot to {filename}")
    else:
        plt.show()
    
    plt.close()

def create_mds_summary(embedding, labels, stress, filename=None):
    """
    Create a summary CSV with MDS coordinates.
    
    Parameters:
    -----------
    embedding : numpy.ndarray
        2D coordinates from MDS
    labels : list
        Labels for each point
    stress : float
        MDS stress value
    filename : str, optional
        Path to save CSV
    
    Returns:
    --------
    df : pandas.DataFrame
        DataFrame with MDS coordinates
    """
    df = pd.DataFrame({
        'dataset': labels,
        'mds_dim1': embedding[:, 0],
        'mds_dim2': embedding[:, 1],
        'group': ['hypertensive' if i < 5 else 'control' for i in range(len(labels))]
    })
    
    if filename:
        df.to_csv(filename, index=False)
        print(f"Saved MDS coordinates to {filename}")
    
    print(f"\nMDS Stress: {stress:.6f}")
    print("(Lower stress indicates better preservation of original distances)")
    
    return df

def process_single_csv(csv_file, output_prefix, dimension=None, pressure=None):
    """
    Process a single distance matrix CSV file.
    
    Parameters:
    -----------
    csv_file : str
        Path to input CSV file
    output_prefix : str
        Prefix for output files
    dimension : int, optional
        Dimension number (for labeling)
    pressure : str, optional
        Pressure level (for labeling)
    """
    print(f"\n{'='*80}")
    print(f"Processing: {csv_file}")
    print(f"{'='*80}")
    
    # Load distance matrix
    distance_matrix, labels = load_distance_matrix(csv_file)
    print(f"Loaded distance matrix: {distance_matrix.shape}")
    print(f"Datasets: {labels}")
    
    # Perform MDS
    embedding, stress = perform_mds(distance_matrix)
    
    # Create title
    title = 'MDS Visualization of Bottleneck Distances'
    if dimension is not None and pressure is not None:
        title += f' (Dimension {dimension}, Pressure {pressure})'
    elif dimension is not None:
        title += f' (Dimension {dimension})'
    
    # Plot MDS
    plot_filename = f'{output_prefix}_mds_plot.png'
    plot_mds_2d(embedding, labels, title=title, filename=plot_filename)
    
    # Save coordinates
    csv_filename = f'{output_prefix}_mds_coordinates.csv'
    df = create_mds_summary(embedding, labels, stress, filename=csv_filename)
    
    return embedding, stress, df

# Example usage
if __name__ == "__main__":
    import os
    
    # Example: Process a single CSV file
    # Modify this path to point to your actual CSV file
    csv_file = 'distances_dim2_p4_fullComplement.csv'
    
    if os.path.exists(csv_file):
        # Process single file
        embedding, stress, df = process_single_csv(
            csv_file,
            output_prefix='mds_results/dim2_p4',
            dimension=2,
            pressure='4'
        )
        
        print("\n" + "="*80)
        print("MDS Analysis Complete!")
        print("="*80)
        print("\nMDS Coordinates:")
        print(df.to_string(index=False))
    else:
        print(f"File not found: {csv_file}")
        print("\nTo use this script:")
        print("1. Update the csv_file variable to point to your distance matrix CSV")
        print("2. Run the script")
        print("\nOr import the functions and use them directly:")
        print("  from mds_analysis import process_single_csv")
        print("  process_single_csv('your_file.csv', 'output_prefix', dimension=0, pressure='1')")
    
    # Example: Process all CSV files in a directory
    print("\n" + "="*80)
    print("To process all dimension/pressure combinations:")
    print("="*80)
    
    #process_all_example = 
    """
# Create output directory
os.makedirs('mds_results', exist_ok=True)

# Process all combinations
for pressure in ['1', '2', '3', '4']:
    for dimension in [0, 1, 2]:
        csv_file = f'bottleneck_results/pressure_{pressure}/distances_dim{dimension}_p{pressure}.csv'
        
        if os.path.exists(csv_file):
            output_prefix = f'mds_results/dim{dimension}_p{pressure}'
            
            embedding, stress, df = process_single_csv(
                csv_file,
                output_prefix=output_prefix,
                dimension=dimension,
                pressure=pressure
            )
"""
    
    #print(process_all_example)