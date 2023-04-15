import random
import pandas as pd
import numpy as np
import networkx as nx

from collections import defaultdict

PORT = 8000


def predict(protein_interactions, references, candidates, iterations=1000):
    graph = nx.from_pandas_edgelist(protein_interactions, 'protein1', 'protein2', 'combined_score')
    #adj_matrix = nx.adjacency_matrix(graph).todense()

    # protein_columns = protein_interactions[['protein1', 'protein2']]
    # proteins = pd.unique(protein_columns.values.ravel('K')).sort()

    # adj_matrix = pd.crosstab(protein_interactions['protein1'],
    #                          protein_interactions['protein2'],
    #                          protein_interactions['combined_score'],
    #                          aggfunc=np.mean).fillna(0)

    #print(adj_matrix)

    node_index = {n: i for n, i in zip(graph.nodes(), range(len(graph.nodes())))}
    num_nodes = len(graph.nodes())

    transition_matrix = np.zeros((num_nodes, num_nodes))

    for edge in graph.edges():
        i = node_index[edge[0]]
        j = node_index[edge[1]]
        transition_matrix[i, j] = 1 / graph.degree(edge[0])

    print(transition_matrix)

    power = np.linalg.matrix_power(transition_matrix, 2)

    print(power)


    # deg_matrix = np.diag(adj_matrix.gt(0).sum(axis=0))
    # deg_matrix = np.diag(adj_matrix.sum(axis=1).A1)

    #deg_inverse = np.linalg.inv(deg_matrix)
    #print(deg_inverse)

    #transition_matrix = np.multiply(deg_matrix, adj_matrix)
    # transition_matrix = adj_matrix / deg_matrix
    #transition_matrix = np.matmul(deg_matrix, adj_matrix)

    density = defaultdict(lambda: 0)

    for reference in references:
        prob_vector = np.zeros(transition_matrix.shape[1]).reshape(-1, 1)
        prob_vector[node_index[reference]] = 1
        print(prob_vector)
        for i in range(iterations):
            if i % 100 == 0:
                print(i)
            prob_vector = np.dot(transition_matrix, prob_vector)
            # density[np.argmax(prob_vector)] += 1

    for gene, score in sorted(density.items(), key=lambda x: x[1], reverse=True):
        print(gene, score)

    subgraph = nx.subgraph(graph, density.keys())
    print(subgraph)
    with open('subgraph.csv', 'w') as f:
        for edge in subgraph.edges():
            score = subgraph.get_edge_data(*edge)['combined_score']
            line = "{},{},{}\n".format(edge[0], edge[1], score)
            f.write(line)


    # print(protein_interactions)
    

def get_edge_value(dictionary, edge):
    return dictionary.get(edge, dictionary.get((edge[1], edge[0])))


def predict_iterative(protein_interactions,
                      references,
                      candidates,
                      iterations=10000,
                      return_prob=0.2,
                      min_score=None):
    density = defaultdict(lambda: 0)
    graph = nx.from_pandas_edgelist(protein_interactions, 'protein1', 'protein2', 'combined_score')
    scores = nx.get_edge_attributes(graph, 'combined_score')
    for start in references:
        current_node = start
        for i in range(iterations):
            if random.random() <= return_prob:
                current_node = start
            else:
                # neighbors = list(nx.neighbors(graph, current_node))
                edges = list(nx.edges(graph, current_node))
                edge_scores = [get_edge_value(scores, e) for e in edges]
                neighbors = [e[1] for e in edges]
                current_node = random.choices(neighbors, weights=edge_scores, k=1)[0]
            density[current_node] += 1
            
    # if a minimum score is specified we filter
    # out all proteins with score < minimum_score
    if min_score:
        density = {prot: score for prot, score in density.items() if score >= min_score}

    for gene, score in sorted(density.items(), key=lambda x: x[1], reverse=True):
        print(gene, score/iterations)

    subgraph = nx.subgraph(graph, density.keys())

    node_index = {n: i for n, i in zip(subgraph.nodes(), range(len(subgraph.nodes())))}
    
    node_data = {i: {'name': n, 'density': density[n]/iterations} for n, i in node_index.items()}

    print(subgraph)
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





