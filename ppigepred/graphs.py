import logging

import networkx as nx


def remove_disconnected_components(graph, references):
    """Takes a graph and an iterable of reference nodes
    and removes components of the graph that are disconnected from
    the reference nodes."""
    all_connected_nodes = set()
    for reference_node in references:
        connected_nodes = nx.node_connected_component(graph, reference_node)
        all_connected_nodes.update(connected_nodes)
    return nx.induced_subgraph(graph, all_connected_nodes)


def get_protein_graph(protein_interactions, references):
    """Creates a graph from a dataframe with protein interactions.
    The graph won't contain components disconnected from the
    provided reference nodes."""
    graph = nx.from_pandas_edgelist(protein_interactions,
                                    'protein1',
                                    'protein2',
                                    'score')

    # If minimum score filtering is performed and
    # a reference protein is connected only to edges
    # that are removed, it'll be missing in the edgelist.
    # We make sure that reference nodes are added to the graph.
    for reference in references:
        if not graph.has_node(reference):
            graph.add_node(reference)
            # Connect the reference node to itself to
            # make sure there's at least one edge
            # and prevent a division by zero
            # when inverting the adjacency matrix.
            graph.add_edge(reference, reference)

    logging.info("Initial number of nodes: %d", len(graph.nodes))
    filtered_graph = remove_disconnected_components(graph, references)
    logging.info("Nodes after removing disconnected components: %d", len(filtered_graph.nodes))
    return filtered_graph


def filter_graph(graph, probabilities):
    """Returns a subgraph containing only nodes with probability > min_score"""
    subgraph = nx.subgraph(graph, probabilities.keys())
    logging.info(subgraph)
    return subgraph


def export_to_dot(graph, filename):
    dot = nx.nx_pydot.to_pydot(graph)
    with open(filename, 'w') as file:
        file.write(str(dot))