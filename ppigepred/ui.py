import argparse
import logging
import pandas as pd
import numpy as np
import time

from .prediction import predict, predict_iterative
from .graphs import get_protein_graph

logging.basicConfig(filename='logs.log', level=logging.DEBUG)

class UI:
    parser = argparse.ArgumentParser(
            description='Predicts genes related to disease based ' + \
                    'on the protein-protein interactions they participate in.')

    parser.add_argument('--db',
                        help='Protein-protein interaction db. CSV with PROT1,PROT2,SCORE',
                        required=True)
    parser.add_argument('-rp',
                        '--restart-probability',
                        help='Probability of returning back to the initial node when doing the random walk simulation',
                        default=0.2)
    parser.add_argument('-i',
                        '--iterations',
                        help='Number of random walk simulations for each reference node',
                        default=1000)
    parser.add_argument('-mis',
                        '--min-interaction-score',
                        help='Will ignore protein interactions with combined score below the specified value. ' + \
                             '(this filters the input network on which the random walk is performed)',
                        default=10)
    parser.add_argument('-m',
                        '--min-score',
                        help='Minimum score to filter out insignificant results.',
                        default=10)

    ref_group = parser.add_mutually_exclusive_group(required=True)
    ref_group.add_argument("-r",
                           "--references",
                           help='List of disease-related gene symbols in comma-separated form')
    ref_group.add_argument("-rf",
                           "--references-file",
                           help='File that contains one disease-related gene symbol per line')

    candidate_group = parser.add_mutually_exclusive_group(required=False)
    candidate_group.add_argument("-c",
                                 "--candidates",
                                 help='List of candidate gene symbols in comma-separated form')
    candidate_group.add_argument("-cf",
                                 "--candidates-file",
                                 help='File that contains one candidate-gene symbol per line')


    @classmethod
    def run(cls):
        args = cls.parser.parse_args()

        if args.references:
            references = args.references.split(',')
        else:
            references = pd.read_csv(args.references_file)

        if args.candidates:
            candidates = args.candidates.split(',')
        elif args.candidates_file:
            candidates = pd.read_csv(args.candidates_file)
        else:
            candidates = set()

        protein_interactions = pd.read_csv(args.db)
        
        protein_graph = get_protein_graph(protein_interactions, references)

        # print(len(protein_interactions))
        # if args.min_interaction_score:
        #     min_interaction_score = int(args.min_interaction_score)
        #     relevant_interactions = protein_interactions['combined_score'] > min_interaction_score
        #     protein_interactions = protein_interactions[relevant_interactions]
        # print(len(protein_interactions))

        # print(np.min(protein_interactions['combined_score']))
        # print(np.mean(protein_interactions['combined_score']))
        # print(np.median(protein_interactions['combined_score']))
        # print(np.max(protein_interactions['combined_score']))

        t = time.time()
        node_index, edges = predict(protein_graph,
                          references,
                          candidates,
                          iterations=int(args.iterations),
                          min_score=float(args.min_score))
        print(time.time() - t)
        return node_index, edges

        # t = time.time()
        # node_index, edges = predict_iterative(protein_graph,
        #                             references,
        #                             candidates,
        #                             iterations=int(args.iterations),
        #                             return_prob=float(args.restart_probability),
        #                             min_score=float(args.min_score))
        # print(time.time() - t)
        # return node_index, edges

