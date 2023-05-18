import logging
import random

from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
from functools import partial

import networkx as nx
import numpy as np
import pandas as pd
import scipy

from .bioseqs import protein_ids_to_gene_ids

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


def filter_probabilities(probabilities, candidates, min_score):
    """Takes a dict of node probabilities and returns
    only the probabilities that match all provided filters.
    Filters out all probabilities that:
    - Have a score lower than the provided one.
    - Are not in the candidates set."""
    filtered_probabilities = probabilities
    if min_score:
        filtered_probabilities = {v: p for v, p in probabilities.items()
                                  if p >= min_score}

    if candidates is not None:
        filtered_probabilities = {v: p for v, p in filtered_probabilities.items()
                                  if v in candidates}

    return filtered_probabilities


def convert_to_gene_probabilities(probabilities, ensembl, graph, top=None):
    """Converts a dictionary of protein probabilities to
    a new dictionary of gene probabilities.
    If the "top=N" parameter is provided it would
    only take the top N genes"""
    prot_to_gene = protein_ids_to_gene_ids(ensembl, probabilities)
    graph = nx.relabel_nodes(graph, prot_to_gene)
    gene_probabilities = OrderedDict()
    for prot_id, prob in probabilities.items():
        gene_id = prot_to_gene[prot_id]
        if gene_id not in gene_probabilities:
            gene_probabilities[gene_id] = prob
        elif prob > gene_probabilities[gene_id]:
            gene_probabilities[gene_id] = prob
        if top and len(gene_probabilities) >= top:
            break
    return gene_probabilities