import gudhi
import numpy as np
import pandas as pd
import scipy.io
from collections import defaultdict, deque

# ── Data extractor functions ──────────────────────────────────────────────────

def connectivityExtractor(name, pruned):
    file_path = ('Networks/Network_Vessels_' if pruned == 0 else 'Pruned/Pruned_Network_') + name + '.mat'
    matlab_data = scipy.io.loadmat(file_path)
    data_structure = matlab_data['Data']
    connectivity_data = data_structure['connectivity'][0, 0].squeeze()
    connectivity_df = pd.DataFrame(connectivity_data, columns=['Parent', 'Daughter1', 'Daughter2', 'Daughter3'])
    connectivity_df.replace(0, np.nan, inplace=True)
    connectivity_df.at[0, 'Parent'] = 0
    return connectivity_df

def nodesExtractor(name, pruned):
    file_path = ('Networks/Network_Vessels_' if pruned == 0 else 'Pruned/Pruned_Network_') + name + '.mat'
    matlab_data = scipy.io.loadmat(file_path)
    data_structure = matlab_data['nodesC2'] if pruned == 0 else matlab_data['nodesC3']
    nodes_df = pd.DataFrame(data_structure.squeeze(), columns=['NodeID', 'X', 'Y', 'Z', 'Degree'])
    return nodes_df

def edgesExtractor(name, pruned):
    file_path = ('Networks/Network_Vessels_' if pruned == 0 else 'Pruned/Pruned_Network_') + name + '.mat'
    matlab_data = scipy.io.loadmat(file_path)
    edge_df = pd.DataFrame(matlab_data['segments'].squeeze(), columns=['ID', 'From', 'To'])
    return edge_df

def findInputVessel(segments, fromnode, to):
    vessel = segments[
        ((segments['From'] == fromnode) & (segments['To'] == to)) |
        ((segments['From'] == to) & (segments['To'] == fromnode))
    ]
    return int(vessel['ID'].iloc[0])

def mapIDExtractor(name, pruned):
    file_path = ('Networks/Network_Vessels_' if pruned == 0 else 'Pruned/Pruned_Network_') + name + '.mat'
    matlab_data = scipy.io.loadmat(file_path)
    data_structure = matlab_data['Data']
    map_df = pd.DataFrame(data_structure['mapIDs'][0, 0].squeeze(), columns=['New', 'Old'])
    return map_df

def lobeExtractor(name, vesID, pruned):
    data = connectivityExtractor(name, pruned)
    tree = defaultdict(list)
    for _, row in data.iterrows():
        parent = row['Parent']
        for daughter_col in ['Daughter1', 'Daughter2', 'Daughter3']:
            daughter = row[daughter_col]
            if pd.notna(daughter):
                tree[parent].append(daughter)
    visited = set()
    queue = deque([vesID])
    while queue:
        current = queue.popleft()
        if current not in visited:
            visited.add(current)
            queue.extend(tree.get(current, []))
    visited.discard(vesID)
    downstream_df = data[data['Parent'].isin(visited)]
    return downstream_df

def node_loc(name, lobe_nodes, pruned):
    nodes = nodesExtractor(name, pruned)
    lobe = nodes[nodes['NodeID'].isin(lobe_nodes)]
    return lobe[['X', 'Y', 'Z']]

def lobeTermLoc_with_edges(name, fromnode, tonode, pruned):
    """Returns (points Nx3, edges list-of-tuples) with 0-based indices."""
    segments = edgesExtractor(name, pruned)
    maps = mapIDExtractor(name, pruned)
    vesID = findInputVessel(segments, fromnode, tonode)
    newID = int(maps[maps['Old'] == vesID]['New'].iloc[0])
    lobe_ves = lobeExtractor(name, newID, pruned)
    new_lobe_ves_ID = lobe_ves['Parent'].to_numpy()
    oldID = maps[maps['New'].isin(new_lobe_ves_ID)]['Old'].to_numpy()
    lobe_segs = segments[segments['ID'].isin(oldID)]
    fromnodes = lobe_segs['From'].to_numpy()
    tonodes   = lobe_segs['To'].to_numpy()
    lobe_node_ids = np.unique(np.concatenate((fromnodes, tonodes))).astype(int)
    
    # Get coordinates
    node_coords = node_loc(name, lobe_node_ids, pruned).to_numpy() / 1000
    
    # Build a map from original node ID -> 0-based row index
    id_to_idx = {nid: i for i, nid in enumerate(lobe_node_ids)}
    
    # Build edge list in 0-based indices
    edges = [(id_to_idx[f], id_to_idx[t])
             for f, t in zip(fromnodes.astype(int), tonodes.astype(int))
             if f in id_to_idx and t in id_to_idx]
    
    return node_coords, edges

def twoInputVessels_with_edges(name, from1, to1, from2, to2, pruned):
    pts1, edg1 = lobeTermLoc_with_edges(name, from1, to1, pruned)
    pts2, edg2 = lobeTermLoc_with_edges(name, from2, to2, pruned)
    # Combine, re-index
    all_pts = np.concatenate((pts1, pts2), axis=0)
    unique_pts, inv = np.unique(all_pts, axis=0, return_inverse=True)
    offset = len(pts1)
    edges_combined = [(inv[a], inv[b]) for a, b in edg1] + \
                     [(inv[a + offset], inv[b + offset]) for a, b in edg2]
    return unique_pts, edges_combined

def get_points(name, coords, pruned):
    if len(coords) == 2:
        return lobeTermLoc_with_edges(name, coords[0], coords[1], pruned)
    elif len(coords) == 4:
        return twoInputVessels_with_edges(name, coords[0], coords[1], coords[2], coords[3], pruned)
    else:
        raise ValueError(f"coords must have 2 or 4 values, got {len(coords)}")

def run_height_persistence(points, edges, direction, height_col=2):
    height_values = points[:, height_col].copy().astype(float)
    if direction == 'neg':
        height_values = height_values.max() + height_values.min() - height_values

    st = gudhi.SimplexTree()
    # Insert vertices with their height as filtration value
    for i, h in enumerate(height_values):
        st.insert([i], filtration=float(h))
    # Insert edges with filtration = max of the two endpoint heights
    for u, v in edges:
        st.insert([u, v], filtration=max(height_values[u], height_values[v]))

    st.make_filtration_non_decreasing()
    st.compute_persistence()
    return st.persistence()

def analyze_persistence(persistence, dimension=0):
    lifespans = [death - birth
                 for dim, (birth, death) in persistence
                 if dim == dimension and death != float('inf')]
    if not lifespans:
        return {'num_bars': 0, 'complexity': 0.0}
    return {'num_bars': len(lifespans), 'complexity': sum(lifespans)}

# ── Pruned flag per dataset ───────────────────────────────────────────────────

pruned_flag = {
    'm1053007': 1, 'm2053007': 1, 'm1053107': 1, 'm2053107': 1,
    'm1060107': 1, 'm1060407': 0, 'm2060407': 0, 'm3060407': 0,
    'm1060507': 0, 'm2060507': 0, 'm3060507': 0,
    'm2060607': 0, 'm3060607': 0,
}

# ── Lobe input coordinates ────────────────────────────────────────────────────
# [from, to] -> lobeTermLoc | [from1, to1, from2, to2] -> twoInputVessels

lobe_coords = {
    '1': {
        'left': {
            'm1053007': [1831, 1858], 'm2053007': [2367, 2368], 'm1053107': [3954, 3900],
            'm2053107': [2866, 2868], 'm1060107': [190, 2722],  'm1060407': [2272, 2273],
            'm2060407': [2776, 2774], 'm3060407': [2473, 2472], 'm1060507': [2121, 2040],
            'm2060507': [2257, 2258], 'm3060507': [692, 2524],  'm2060607': [1475, 1997],
            'm3060607': [53, 2576],
        },
        'middle': {
            'm1053007': [1841, 1945, 1841, 1964], 'm2053007': [2463, 2497, 2502, 2501],
            'm1053107': [4179, 4304, 4018, 4379], 'm2053107': [2920, 2921, 3150, 3149],
            'm1060107': [2766, 2744, 2992, 1332], 'm1060407': [2456, 2459, 2455, 1061],
            'm2060407': [2963, 2568, 2964, 2995], 'm3060407': [2509, 2613, 2507, 408],
            'm1060507': [1992, 1929, 2030, 2065], 'm2060507': [2390, 2227, 2284, 2283],
            'm3060507': [2624, 2377, 2568, 1365], 'm2060607': [1916, 2052, 2057, 1991],
            'm3060607': [2604, 2466, 2677, 2607],
        },
        'superior': {
            'm1053007': [1836, 1835],             'm2053007': [2464, 2406],
            'm1053107': [4071, 685],               'm2053107': [2867, 2979],
            'm1060107': [2780, 2716, 2780, 2781],  'm1060407': [2274, 2283],
            'm2060407': [2742, 2598],              'm3060407': [2418, 2419, 2420, 2491],
            'm1060507': [1993, 1997],              'm2060507': [2259, 2330],
            'm3060507': [2392, 2398],              'm2060607': [1914, 1915],
            'm3060607': [2603, 2518],
        },
        'inferior': {
            'm1053007': [1841, 1839], 'm2053007': [2502, 2501], 'm1053107': [4018, 4019],
            'm2053107': [3150, 3170], 'm1060107': [2992, 2895], 'm1060407': [2455, 2394],
            'm2060407': [2964, 2727], 'm3060407': [2507, 2508], 'm1060507': [2030, 2027],
            'm2060507': [2391, 2423], 'm3060507': [2537, 1693], 'm2060607': [2057, 2056],
            'm3060607': [2677, 2676],
        },
        'postcaval': {
            'm1053007': [1864, 35],   'm2053007': [2465, 2558], 'm1053107': [4296, 1851],
            'm2053107': [2915, 2913], 'm1060107': [2766, 2808], 'm1060407': [2454, 2602],
            'm2060407': [2886, 2962], 'm3060407': [2421, 2469], 'm1060507': [1994, 2178],
            'm2060507': [2284, 2454], 'm3060507': [2628, 2355], 'm2060607': [1916, 1970],
            'm3060607': [2701, 2688],
        },
    },
    '2': {
        'left': {
            'm1053007': [1075, 1960], 'm2053007': [217, 1868],  'm1053107': [1506, 2897],
            'm2053107': [44, 1616],   'm1060107': [2089, 2169], 'm1060407': [1225, 1265],
            'm2060407': [1371, 1392], 'm3060407': [1882, 1924], 'm1060507': [118, 1381],
            'm2060507': [1509, 1519], 'm3060507': [51, 1591],   'm2060607': [1495, 1419],
            'm3060607': [535, 1916],
        },
        'middle': {
            'm1053007': [1841, 1957, 1841, 2025], 'm2053007': [1831, 1854, 1908, 1919],
            'm1053107': [2738, 2739, 2959, 1567], 'm2053107': [1645, 1644, 1642, 1641],
            'm1060107': [2164, 2098, 2153, 927],  'm1060407': [1289, 1269, 1253, 1252],
            'm2060407': [1464, 1466, 1465, 592],  'm3060407': [1907, 1962, 1939, 337],
            'm1060507': [1322, 1400, 1320, 1425], 'm2060507': [1565, 1566, 1578, 1577],
            'm3060507': [1732, 1589, 1749, 1109], 'm2060607': [1441, 1221, 1460, 221],
            'm3060607': [1925, 1946, 2002, 2094],
        },
        'superior': {
            'm1053007': [1829, 1828],             'm2053007': [1829, 1830],
            'm1053107': [2909, 2895],             'm2053107': [1640, 1693],
            'm1060107': [2162, 2163],             'm1060407': [1266, 1239],
            'm2060407': [1387, 1405],             'm3060407': [1883, 1890, 1880, 1881],
            'm1060507': [1468, 1454],             'm2060507': [1528, 1552],
            'm3060507': [1706, 1705],             'm2060607': [1440, 1426],
            'm3060607': [1923, 1924],
        },
        'inferior': {
            'm1053007': [1842, 2020], 'm2053007': [1908, 1979], 'm1053107': [2959, 2957],
            'm2053107': [1642, 1755], 'm1060107': [2153, 2151], 'm1060407': [1253, 1252],
            'm2060407': [1465, 1515], 'm3060407': [1939, 1982], 'm1060507': [1320, 1424],
            'm2060507': [1620, 1602], 'm3060507': [1749, 1809], 'm2060607': [1460, 1375],
            'm3060607': [2002, 1903],
        },
        'postcaval': {
            'm1053007': [1862, 2003], 'm2053007': [1907, 1909], 'm1053107': [3056, 1829],
            'm2053107': [1724, 1722], 'm1060107': [2164, 2190], 'm1060407': [1290, 1323],
            'm2060407': [1386, 1410], 'm3060407': [1907, 1796], 'm1060507': [1319, 1321],
            'm2060507': [1578, 1650], 'm3060507': [1739, 1689], 'm2060607': [1441, 1509],
            'm3060607': [1951, 1950],
        },
    },
    '3': {
        'left': {
            'm1053007': [1630, 1628], 'm2053007': [1619, 1621], 'm1053107': [246, 2864],
            'm2053107': [717, 1751],  'm1060107': [318, 1645],  'm1060407': [1122, 1121],
            'm2060407': [808, 814],   'm3060407': [919, 921],   'm1060507': [80, 1354],
            'm2060507': [1097, 1116], 'm3060507': [191, 959],   'm2060607': [653, 620],
            'm3060607': [142, 1171],
        },
        'middle': {
            'm1053007': [1721, 1720, 1721, 95],   'm2053007': [1731, 1737, 1682, 1733],
            'm1053107': [3103, 2859, 2892, 844],  'm2053107': [1621, 1622, 1714, 1709],
            'm1060107': [1651, 1727, 1611, 1610], 'm1060407': [1110, 1109, 1140, 426],
            'm2060407': [854, 283, 891, 50],       'm3060407': [968, 938, 967, 976],
            'm1060507': [1470, 1475, 1469, 1370], 'm2060507': [1041, 1121, 1042, 1163],
            'm3060507': [1000, 451, 1062, 1016],  'm2060607': [607, 595, 616, 640],
            'm3060607': [1261, 1202, 1257, 1268],
        },
        'superior': {
            'm1053007': [1657, 1656],            'm2053007': [1620, 1703],
            'm1053107': [2982, 2955],            'm2053107': [1612, 1610],
            'm1060107': [1644, 1650],            'm1060407': [1126, 1095],
            'm2060407': [862, 832],              'm3060407': [969, 970, 920, 924],
            'm1060507': [1405, 1461],            'm2060507': [1096, 1095],
            'm3060507': [952, 950],              'm2060607': [615, 642],
            'm3060607': [1172, 1255],
        },
        'inferior': {
            'm1053007': [1721, 1624], 'm2053007': [1682, 1683], 'm1053107': [2892, 1313],
            'm2053107': [1714, 1715], 'm1060107': [1611, 1757], 'm1060407': [1140, 1157],
            'm2060407': [891, 844],   'm3060407': [967, 975],   'm1060507': [1469, 1488],
            'm2060507': [1040, 1043], 'm3060507': [1062, 1063], 'm2060607': [616, 637],
            'm3060607': [1257, 1284],
        },
        'postcaval': {
            'm1053007': [1755, 1636], 'm2053007': [1704, 1732], 'm1053107': [3039, 1474],
            'm2053107': [1655, 1653], 'm1060107': [1651, 29],   'm1060407': [1134, 1141],
            'm2060407': [851, 853],   'm3060407': [968, 1020],  'm1060507': [1381, 1379],
            'm2060507': [1042, 1136], 'm3060507': [1009, 1057], 'm2060607': [607, 664],
            'm3060607': [1250, 1249],
        },
    },
    '4': {
        'left': {
            'm1053007': [1810, 1811], 'm2053007': [1209, 1214], 'm1053107': [1252, 2103],
            'm2053107': [1239, 1392], 'm1060107': [215, 1664],  'm1060407': [1065, 1013],
            'm2060407': [694, 701],   'm3060407': [851, 926],   'm1060507': [1202, 1169],
            'm2060507': [1006, 1064], 'm3060507': [395, 975],   'm2060607': [511, 737],
            'm3060607': [1223, 1222],
        },
        'middle': {
            'm1053007': [1954, 1803, 1954, 1788], 'm2053007': [1223, 161, 1273, 1342],
            'm1053107': [2059, 2186, 2058, 1054], 'm2053107': [1520, 1532, 1555, 1565],
            'm1060107': [1652, 1619, 1653, 1529], 'm1060407': [1055, 1054, 1084, 516],
            'm2060407': [692, 127, 690, 313],      'm3060407': [865, 61, 878, 880],
            'm1060507': [1243, 1244, 1212, 1341], 'm2060507': [944, 993, 945, 969],
            'm3060507': [959, 960, 1003, 1093],   'm2060607': [753, 768, 720, 719],
            'm3060607': [1290, 270, 1207, 1309],
        },
        'superior': {
            'm1053007': [1819, 1818],            'm2053007': [1213, 1259],
            'm1053107': [2125, 2154],            'm2053107': [1391, 1477],
            'm1060107': [1574, 1522],            'm1060407': [1053, 1051],
            'm2060407': [695, 727],              'm3060407': [885, 894, 884, 864],
            'm1060507': [1232, 1233],            'm2060507': [999, 983],
            'm3060507': [1027, 1026],            'm2060607': [738, 741],
            'm3060607': [1268, 1267],
        },
        'inferior': {
            'm1053007': [1983, 1923], 'm2053007': [1273, 1281], 'm1053107': [2058, 1074],
            'm2053107': [1555, 1494], 'm1060107': [1653, 1687], 'm1060407': [1084, 1098],
            'm2060407': [690, 691],   'm3060407': [878, 879],   'm1060507': [1212, 1213],
            'm2060507': [943, 946],   'm3060507': [1003, 1005], 'm2060607': [720, 769],
            'm3060607': [1207, 1205],
        },
        'postcaval': {
            'm1053007': [1903, 511],  'm2053007': [1224, 1361], 'm1053107': [2057, 1273],
            'm2053107': [1550, 1556], 'm1060107': [1652, 1713], 'm1060407': [1059, 1083],
            'm2060407': [783, 739],   'm3060407': [865, 920],   'm1060507': [1234, 1326],
            'm2060507': [945, 1036],  'm3060507': [986, 985],   'm2060607': [753, 752],
            'm3060607': [1207, 1274],
        },
    },
}

# ── Main loop ─────────────────────────────────────────────────────────────────

all_metrics = []

for pressure in ['1', '2', '3', '4']:
    for lobe in ['left', 'middle', 'superior', 'inferior', 'postcaval']:
        for mouse, coords in lobe_coords[pressure][lobe].items():

            pruned = pruned_flag[mouse]
            name = f'{mouse[:2]}p{pressure}_{mouse[2:]}'  # e.g. m1053007 -> m1p1_053007

            points, edges = get_points(name, coords, pruned)

            post_persistence = run_height_persistence(points, edges, direction='pos', height_col=2)
            post_results     = analyze_persistence(post_persistence, dimension=0)

            ant_persistence  = run_height_persistence(points, edges, direction='neg', height_col=2)
            ant_results      = analyze_persistence(ant_persistence, dimension=0)

            all_metrics.append({
                'mouse_name':             mouse,
                'pressure':               pressure,
                'lobe':                   lobe,
                'anterior_branch_count':  ant_results['num_bars'],
                'anterior_complexity':    ant_results['complexity'],
                'posterior_branch_count': post_results['num_bars'],
                'posterior_complexity':   post_results['complexity'],
            })

            print(f"  Done: {name} | {lobe} | ant={ant_results['num_bars']} | post={post_results['num_bars']}")

# ── Save ──────────────────────────────────────────────────────────────────────

metrics_df = pd.DataFrame(all_metrics, columns=[
    'mouse_name', 'pressure', 'lobe',
    'anterior_branch_count', 'anterior_complexity',
    'posterior_branch_count', 'posterior_complexity',
])
metrics_df.to_csv('directional_complexity_lobe.csv', index=False)
print(f"\nSaved to directional_complexity_lobe.csv — {len(metrics_df)} rows")
print(metrics_df.head(10))