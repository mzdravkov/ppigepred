import networkx as nx


def remove_disconnected_components(graph, references):
    all_connected_nodes = set()
    for reference_node in references:
        connected_nodes = nx.node_connected_component(graph, reference_node)
        all_connected_nodes.update(connected_nodes)
    return nx.induced_subgraph(graph, all_connected_nodes)


def get_protein_graph(protein_interactions, references):
    graph = nx.from_pandas_edgelist(protein_interactions,
                                    'protein1',
                                    'protein2',
                                    'combined_score')
    print("Initial number of nodes:", len(graph.nodes))
    filtered_graph = remove_disconnected_components(graph, references)
    print("Nodes after removing disconnected components:", len(filtered_graph.nodes))
    return filtered_graph
    
    