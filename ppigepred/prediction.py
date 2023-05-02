import logging
import random

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import partial

import networkx as nx
import numpy as np
import pandas as pd
import scipy


def simulate_random_walks(return_prob, transition_transp, init_prob, num_walks):
    """Performs <num_walks> number of random walks (with a restart)
    and returns the probability vector obtained."""
    prob = np.zeros(len(init_prob))
    for _ in range(num_walks):
        current_walk_prob = init_prob
        while random.random() >= return_prob:
            current_walk_prob = transition_transp.dot(current_walk_prob)
        prob += current_walk_prob / num_walks
    return prob


def predict(graph,
            references,
            walks=1000,
            return_prob=0.05,
            min_score=None,
            processes=1):
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
    # The total number of random walks is split into chunks that are executed
    # in parallel.
    prob = np.zeros(len(graph.nodes))
    random_walk_fn = partial(simulate_random_walks,
                             return_prob,
                             transition_transp,
                             init_prob)
    chunks = [walks // processes] * processes
    if walks % processes != 0:
        chunks[-1] += walks % processes
    with ProcessPoolExecutor(max_workers=processes) as executor:
        for walks_chunk_prob in executor.map(random_walk_fn, chunks):
            prob += walks_chunk_prob / processes

    logging.info('Probability vector: %s', pd.DataFrame(prob).describe())

    # if a minimum score is specified we filter
    # out all proteins with score < minimum_score
    selected_nodes = graph.nodes
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


    # with open('subgraph.csv', 'w') as f:
    #     for edge in subgraph.edges():
    #         score = subgraph.get_edge_data(*edge)['combined_score']
    #         line = "{},{},{}\n".format(edge[0], edge[1], score)
    #         f.write(line)



    # dot = nx.nx_pydot.to_pydot(subgraph)
    # with open('subgraph.dot', 'w') as f:
    #     f.write(str(dot))

