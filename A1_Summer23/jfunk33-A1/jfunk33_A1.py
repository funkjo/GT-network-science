def practice_graphs():
    
    """
    Returns:
    
    G: NetworkX graph object
    """
    
    # This line is a placeholder
    G = nx.Graph()
    G.add_nodes_from([0,1,2,3,4,5,6,7,8])
    G.add_edges_from([(1, 2), 
                      (1, 3),
                      (2, 3),
                      (4, 5),
                      (4, 6),
                      (4, 7),
                      (4, 8)])
    pos=nx.spring_layout(G, k=3, center=[4,1])
    
    return G, pos

def create_toy_graphs():
    
    """
    
    Returns:
    cycle: a networkx graph object meeting the requirements of a cycle
    clique: a networkx graph object meeting the requirements of a clique
    star: a networkx graph object meeting the requirements of a star network
    
    """
    
    #These lines are placeholders
    cycle = nx.cycle_graph(10)
    clique = nx.complete_graph(10)
    star = nx.star_graph(10)
    
    return cycle, clique, star

def calculate_leading_eigenvalue(G):
    """
    Inputs:
    G: NetworkX graph object
    
    Returns:
    eig (float): leading eigenvalue of the adjacency matrix 
    
    """
    
    # This is a placeholder
    A = nx.to_numpy_array(G)
    eig = np.linalg.eig(A)
    eig = list(eig[0])
    eig = [abs(ele) for ele in eig]
    
    return max(eig)

def cities_map_analysis(G):
    """
    Input
    G: NetworkX Graph object
    
    Output
    node_count (int): number of nodes
    edge_count (int): number of edges
    city_count (int): number of city pairs within 50 miles of each other
    """
    
    # these are placeholders
    node_count = G.number_of_nodes()
    edge_count = G.number_of_edges()
    city_count = len([(u,v) for u,v,e in G.edges(data=True) if e['weight'] > 50])
    
    return node_count, edge_count, city_count

def cities_within_100(G, city_list):
    '''
    Input
    G: NetworkX graph object
    city_list: list of strings (names of cities in G)

    Output
    S: subgraph of G that only contains edges between cities in “city_list” and directly neighboring cities that are less than 100 miles away
    '''

    # This is a placeholder
#     S = nx.Graph()
    print(city_list)
    S = nx.Graph(((u, v, e) for u,v,e in G.edges(data=True) if (u in city_list or v in city_list) and e['weight'] < 100))
    print(S.nodes(data=True))
    
    return S

def les_mis_connected(G):
    """
    Returns:
    les_mis_connected (Boolean): whether the graph is connected or not
    
    """
    
    return nx.is_connected(G)

def calculate_average_paths(G):
    """
    Inputs
    G: NetworkX graph object
    
    Returns:
    aspl (float): average shortest path length
    mspl (int): maximum shortest path length
    """
    
    combos = set(itertools.combinations(G.nodes(), 2))
    shortest_paths = []
    for combo in combos:
        shortest_path = nx.shortest_path(G, source=combo[0], target=combo[1])
        shortest_paths.append(len(shortest_path)-1)
    
    
    # These lines are placeholders
    aspl = sum(shortest_paths) / len(shortest_paths)
    mspl = max(shortest_paths)
    
    sns.countplot(x=shortest_paths)
    
    return aspl, mspl

def popular_characters(G):
    
    """
    Inputs:
    G: NetworkX graph object
    
    Returns:
    top_3: List[str]: a list of the string names of the 3 most commonly mentioned characters in the novel
    
    """
    
    A = nx.adjacency_matrix(G, weight='value')
    A = A.todense()
#     print(A)
    A = A.transpose()
    A = np.array(A)
    T = A
#     print(T)
    q0 = [1/len(G.nodes())] * len(G.nodes())
#     print(q0)
    q = np.array(q0)
    for i in range(150):
        q = np.matmul(T, q)
#         print(q)

    popular = pd.DataFrame()
    popular['Character'] = np.array(G.nodes())
    popular['StationaryDistribution'] = q
    popular = popular.sort_values(by=['StationaryDistribution'])
    top_3 = list(popular['Character'].head(3))
        
    
    
    
    
    # This line is a placeholder
#     top_3 = ["George", "P.", "Burdell"]
    
    return top_3

def weakly_connected(G):
    """
    Inputs
    G: NetworkX graph object
    
    Returns:
    wccs (int): number of weakly connected components in the graph
    wpct (float): percent of nodes in the graph that belong to the largest weakly connected component
    
    """
    
    largest = max(nx.weakly_connected_components(G), key=len)
    
    # These lines are placeholders
    wccs = nx.number_weakly_connected_components(G)
    wpct = len(largest) / len(G.nodes())
    
    return wccs, wpct

def strongly_connected(G):
    """
    Inputs
    G: NetworkX graph object
    
    Returns:
    sccs (int): number of strongly connected components in the graph
    spct (float): percent of nodes in the graph that belong to the largest strongly connected component
    
    """
    
    largest = max(nx.strongly_connected_components(G), key=len)
    
    #These lines are placeholders
    sccs = nx.number_strongly_connected_components(G)
    spct = len(largest) / len(G.nodes())
    
    return sccs, spct

def scc_shortest_paths(G):
    """
    Inputs
    G: NetworkX graph object
    
    Returns:
    aspl (float): average shortest path length of the largest strongly connected component
    mspl (int): maximum shortest path length of the largest strongly connected component
    
    """
    
    largest = max(nx.strongly_connected_components(G), key=len)
    
    combos = set(itertools.combinations(G.subgraph(largest).nodes(), 2))
    shortest_paths = []
    for combo in combos:
        shortest_path = nx.shortest_path(G.subgraph(largest), source=combo[0], target=combo[1])
        shortest_paths.append(len(shortest_path)-1)
    
    
    # These lines are placeholders
    aspl = sum(shortest_paths) / len(shortest_paths)
    mspl = max(shortest_paths)
    
    sns.countplot(x=shortest_paths)
    
    return aspl, mspl

def is_graph_dag(G):
    """
    Input:
    G: NetworkX graph object
    
    Returns:
    is_dag (Boolean): is the graph directed or not
    """
    
    # This is a placeholder
    is_dag = True
    is_dag = False
    
    return nx.is_directed_acyclic_graph(G)

def remove_cycles(G):
    """
    Input:
    G: NetworkX graph object
    
    Returns:
    G2: NetworkX graph object with cycles removed 
    
    """
    G2 = G
    
    loop = True
    while loop:
        try:
            cycle = nx.find_cycle(G2, orientation="original")
            edge = cycle[0][:2]
            G2.remove_edge(edge[0], edge[1])
        except Exception as e:
            loop = False
    
    return G2

def get_sources(G):
    """
    Inputs:
    G: NetworkX graph object
    
    Returns:
    num_sources (int): number of source nodes
    top_source (str): the name of the node with the highest level of influence
    """
    
    #These are placeholders
    num_sources = len([x for x in G.nodes() if G.in_degree(x)==0])
    
    _max = 0
    top_source = "George P. Burdell"
    for source in [x for x in G.nodes() if G.in_degree(x)==0]:
        if len(nx.descendants(G, source)) > _max:
            _max = len(nx.descendants(G, source))
            top_source = source
    
    
    return num_sources, top_source

def get_sources(G):
    """
    Inputs:
    G: NetworkX graph object
    
    Returns:
    num_sources (int): number of source nodes
    top_source (str): the name of the node with the highest level of influence
    """
    
    #These are placeholders
    num_sources = len([x for x in G.nodes() if G.in_degree(x)==0])
    
    _max = 0
    top_source = "George P. Burdell"
    for source in [x for x in G.nodes() if G.in_degree(x)==0]:
        if len(nx.descendants(G, source)) > _max:
            _max = len(nx.descendants(G, source))
            top_source = source
    
    
    return num_sources, top_source

def calculated_projections(G):
    """
    Inputs:
    G: NetworkX graph object
    
    Returns:
    user_matrix: one mode projection for users
    project_matrix: one mode projection for projects
    """
    
    user_nodes = [x for x,y in G.nodes(data=True) if y['node_type']=='user']
    project_nodes = [x for x,y in G.nodes(data=True) if y['node_type']=='project']
    u = nx.bipartite.biadjacency_matrix(G, row_order = user_nodes, column_order = project_nodes)
    p = nx.bipartite.biadjacency_matrix(G, row_order = project_nodes, column_order = user_nodes)
    
    user_matrix = u @ u.T
    project_matrix = p @ p.T
    
    return user_matrix, project_matrix

def get_user_pair(M):
    """
    Inputs:
    M: projected matrix
    
    Returns:
    
    u1 (str) - first user
    u2 (str) - second user
    """
    
    u1 = "George123"
    u2 = "Burdell456"
    
    return u1, u2

def get_project_pair(G):
    """
    Inputs:
    G: projected matrix
    
    Returns:
    
    p1 (str) - first project
    p2 (str) - second project
    """
    
    p1 = "George123"
    p2 = "Burdell456"
    
    return p1, p2
