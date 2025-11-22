def blog_graph_degree_adjacent_nodes(filename="blog.txt"):
    """
    Parse the blog.txt file into an undirected graph and then return the degree of adjacent nodes

    Input: filename

    Output: x, y node degree-degree pairs for edges in G
    """
    # This is a placeholder
    G = nx.read_edgelist(filename,create_using=nx.DiGraph())
    U = G.to_undirected()
    x, y = [], []
    for edge in U.edges():
        x.append(U.degree(edge[0]))  # Degree of the first node in the edge
        y.append(U.degree(edge[1]))  # Degree of the second node in the edge

    return x, y

def pearson_correlation_coefficient(x: List[float], y: List[float]) -> Tuple[float, float]:
    """
    Calculates the Pearson correlation coefficient between two vectors.

    Input: undirected graph G

    Output: Pearson correlation coefficient (r) and p-value
    """
    
    result = pearsonr(x, y)
    r, p = result.statistic, result.pvalue

    # your code here

    return r, p

def largest_scc(filename="blog.txt") -> nx.Graph:
    """
    Returns the largest strongly connected component of a directed graph.

    Input: filename for graph data

    Output: largest strongly connected component of G as an undirected graph.
    """
    # This is a placeholder
    # your code here
    G = nx.read_edgelist(filename,create_using=nx.DiGraph())
    largest_scc = max(nx.strongly_connected_components(G), key=len)
    largest_scc = G.subgraph(largest_scc)
    largest_scc = largest_scc.to_undirected()

    return largest_scc

def clustering_coeff(G: nx.Graph) -> Dict[int, float]:
    """
    Calculates the average clustering coefficient of a graph.

    Input: undirected graph G

    Output: average clustering coefficient of G
    """
    # This is a placeholder. It _can_ be deleted or modified.
    acc_dict = nx.clustering(G)

    # your code here

    # the return object should be a dict where each key is a node in the graph and the value is the clustering coefficient of that node.
    return acc_dict

def define_n_p(G: nx.Graph):
    """
    Determine n and p values for a graph

    Input: Undirected graph, G

    Output: n (number of nodes), p (probability of an edge)
    """
    # This is a placeholder. It should be modified.
    n, p = len(G.nodes()), nx.density(G)

    # your code here

    return n, p

def create_gnp(n: int, p: float) -> nx.Graph:
    """
    Creates a random graph with n nodes and edge probability p.

    Input: n (number of nodes), p (probability of an edge)

    Output: a random graph with n nodes and edge probability p
    """
    # This is a placeholder. It should be modified.
    G_r = nx.fast_gnp_random_graph(n,p)

    # your code here

    return G_r

def blog_graph_scc(filename="blog.txt") -> nx.DiGraph():
    """
    Parse the blog.txt file and find the largest strongly connected component and convert to a directed graph.

    Input: filename

    Output: scc, your largest strongly connected component as a directed graph
    """
    # This is a placeholder
    G = nx.read_edgelist(filename,create_using=nx.DiGraph())
    scc = max(nx.strongly_connected_components(G), key=len)
    scc = G.subgraph(scc)
    
    return scc

def count_type_5(scc: nx.DiGraph) -> int:
    """
    Counts the number of type 5 motifs in a graph.
    A type 5 motif is a feed-forward loop: A->B->C with the additional edge A->C

    Input: directed graph scc

    Output: number of type 5 nodes in scc
    """

    # your code here
    triads = nx.triadic_census(scc)
    type_5 = triads['030T']

    return type_5

def count_type_9(scc: nx.DiGraph) -> int:
    """
    Counts the number of type 9 motifs in a graph.
    A type 9 motif is a directed cycle: A->B->C->A

    Input: directed graph scc

    Output: number of type 9 nodes in scc
    """

    # your code here
    triads = nx.triadic_census(scc)
    type_9 = triads['030C']

    return type_9
