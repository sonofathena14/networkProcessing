import re
import gudhi
import numpy as np
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt
from gudhi.representations import Landscape
from collections import defaultdict, deque

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

    columns = ['mouse_name', 'pressure','lobe']
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
        # Key format: '<mouse>_p<pressure>_<lobe>'
        # Use regex to split on _p<digit>_ to avoid false matches in mouse names
        match = re.search(r'_p(\d+)_(.+)$', dataset_name)
        if match:
            mouse_name = dataset_name[:match.start()]
            pressure = match.group(1)
            lobe = match.group(2)
        else:
            mouse_name = dataset_name
            pressure = ''
            lobe = ''

        row = {'mouse_name': mouse_name, 'pressure': pressure, 'lobe': lobe}
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
                row[f'{prefix}_landscape1_L1_norm'] = float(np.trapezoid(np.abs(lam1), t_vals))

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
    df = df.sort_values(['mouse_name', 'pressure', 'lobe']).reset_index(drop=True)
    df.to_csv(output_filename, index=False)
    print(f"\nSaved landscape metrics summary to {output_filename}")
    return df


# ── Main ──────────────────────────────────────────────────────────────────────

all_landscape_curves = {}
all_persistence_diagrams = {}

for pressure in ['1', '2', '3', '4']:
    for lobe in ['left','superior','inferior','middle','postcaval']:
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

        persistence_results = {name: compute_persistence(points) for name, points in datasets.items()}
        persistence_diagrams = {name: res[0] for name, res in persistence_results.items()}
        last_deaths = {name: res[1] for name, res in persistence_results.items()}

        global_last_deaths = {}
        for deaths_by_dim in last_deaths.values():
            for dim, death in deaths_by_dim.items():
                global_last_deaths[dim] = max(global_last_deaths.get(dim, death), death)

        landscape_curves = {
            name: compute_landscape_curves(diag, t_min=0, t_max_per_dim=global_last_deaths, max_dim=2, resolution=500)
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
            full_name = name + '_p' + pressure + '_' + lobe
            all_landscape_curves[full_name] = curves
            all_persistence_diagrams[full_name] = persistence_diagrams[name]

        print("Completed Pressure " + pressure+ " lobe " + lobe)

# Create summary CSV
print("\n=== Creating Landscape Metrics Summary CSV ===")
summary_df = create_landscape_summary_csv(
    all_landscape_curves,
    all_persistence_diagrams,
    output_filename='landscape_metrics_summary_lobe.csv'
)

print(f"\nTotal datasets processed: {len(all_landscape_curves)}")
print(f"Rows in summary: {len(summary_df)}")
print(summary_df.head())