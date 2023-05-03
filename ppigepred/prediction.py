import logging
import random

from collections import OrderedDict
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
            processes=1):
    """Simulates multiple random walks with restart from a set
    of reference nodes to calculate the probability densitiy
    for all nodes in the graph"""
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

    return {v: prob[i] for v, i in node_index.items()}

def filter_probabilities(probabilities, min_score, top):
    """Takes a dict of node probabilities, a min score, and requested number
    of top nodes and returns an ordered dictionary of the top nodes with score
    higher than the provided value."""
    filtered_probabilities = probabilities
    if min_score:
        filtered_probabilities = {v: p for v, p in probabilities.items()
                                  if p >= min_score}

    ordered_proteins = sorted(filtered_probabilities.items(),
                              key=lambda x: x[1],
                              reverse=True)
    limit = top if top else len(filtered_probabilities)
    return OrderedDict(ordered_proteins[:limit])
