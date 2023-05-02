import random
import logging

from collections import defaultdict

import pandas as pd
import numpy as np
import networkx as nx
import scipy


def predict(graph,
            references,
            candidates,
            walks=1000,
            return_prob=0.05,
            min_score=None):
    # Get the adjacency matrix A.
    adj_matrix = nx.adjacency_matrix(graph)
    node_index = {v: i for i, v in enumerate(graph.nodes)}

    # Create a probability vector.
    # The initial vector equal probability for all reference nodes
    # and 0 probability for all others.
    init_prob = np.zeros(len(graph.nodes))
    for node in references:
        init_prob[node_index[node]] = 1 / len(references)

    # Get a diagonal matrix D with the inverted sums of
    # the adjacency matrix rows.
    inv_sum = np.sum(adj_matrix, axis=1).astype(float) ** -1
    deg_inverse = np.diag(inv_sum)

    # The transition matrix is D.A
    transition = scipy.sparse.csr_matrix(deg_inverse).dot(adj_matrix)
    transition_transp = np.transpose(transition)

    # Perform <iterations> random walks to get the final probability vector.
    # Each step of every random walk has <return_prob> chance of ending the walk.
    prob = np.zeros(len(graph.nodes))
    for _ in range(walks):
        current_walk_prob = init_prob
        while random.random() >= return_prob:
            current_walk_prob = transition_transp.dot(current_walk_prob)
        prob += current_walk_prob/walks

    logging.info('Probability vector: %s', pd.DataFrame(prob).describe())

    # if a minimum score is specified we filter
    # out all proteins with score < minimum_score
    selected_nodes = []
    if min_score:
        selected_nodes = [n for n in graph.nodes if prob[node_index[n]] >= min_score]

    subgraph = nx.subgraph(graph, selected_nodes)
    logging.info(subgraph)

    # prepare data for visalization
    selected_node_index = {n: i for i, n in enumerate(subgraph.nodes)}

    node_data = {i: {
        'name': n,
        'density': prob[node_index[n]],
        'is_reference': n in references,
        } for n, i in selected_node_index.items()}

    subgraph_data = []
    for edge in subgraph.edges():
         node1_index = selected_node_index[edge[0]]
         node2_index = selected_node_index[edge[1]]
         score = subgraph.get_edge_data(*edge)['combined_score']
         subgraph_data.append((node1_index, node2_index, score))

    return node_data, subgraph_data


def get_edge_value(dictionary, edge):
    return dictionary.get(edge, dictionary.get((edge[1], edge[0])))


def predict_iterative(graph,
                      references,
                      candidates,
                      walks=10000,
                      return_prob=0.05,
                      min_score=None):
    density = defaultdict(lambda: 0)
    scores = nx.get_edge_attributes(graph, 'combined_score')
    for start in references:
        current_node = start
        density[current_node] += 1
        for _ in range(walks):
            while random.random() >= return_prob:
                edges = list(nx.edges(graph, current_node))
                edge_scores = [get_edge_value(scores, e) for e in edges]
                neighbors = [e[1] for e in edges]
                current_node = random.choices(neighbors, weights=edge_scores, k=1)[0]
                density[current_node] += 1
            current_node = start

    # if a minimum score is specified we filter
    # out all proteins with score < minimum_score
    if min_score:
        density = {prot: score for prot, score in density.items() if score >= min_score}

    # for gene, score in sorted(density.items(), key=lambda x: x[1], reverse=True):
    #     print(gene, score/walks)

    subgraph = nx.subgraph(graph, density.keys())
    logging.info(subgraph)

    node_index = {n: i for n, i in zip(subgraph.nodes(), range(len(subgraph.nodes())))}

    node_data = {i: {
        'name': n,
        'density': density[n]/walks,
        'is_reference': n in references,
        } for n, i in node_index.items()}

    subgraph_data = []
    for edge in subgraph.edges():
         node1_index = node_index[edge[0]]
         node2_index = node_index[edge[1]]
         score = subgraph.get_edge_data(*edge)['combined_score']
         subgraph_data.append((node1_index, node2_index, score))

    return node_data, subgraph_data

    # with open('subgraph.csv', 'w') as f:
    #     for edge in subgraph.edges():
    #         score = subgraph.get_edge_data(*edge)['combined_score']
    #         line = "{},{},{}\n".format(edge[0], edge[1], score)
    #         f.write(line)



    # dot = nx.nx_pydot.to_pydot(subgraph)
    # with open('subgraph.dot', 'w') as f:
    #     f.write(str(dot))





