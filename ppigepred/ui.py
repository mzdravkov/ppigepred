import argparse
import logging
import pandas as pd
import numpy as np
import json
import time

from flask import Flask

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
    parser.add_argument('-w',
                        '--walks',
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
    parser.add_argument('-i',
                        '--interactive',
                        help='Run interactive graph visualization server.',
                        action='store_true')

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

        t = time.time()
        node_index, edges = predict(protein_graph,
                          references,
                          candidates,
                          walks=int(args.walks),
                          min_score=float(args.min_score))
        print(time.time() - t)

        t = time.time()
        node_index, edges = predict_iterative(protein_graph,
                                    references,
                                    candidates,
                                    walks=int(args.walks),
                                    return_prob=float(args.restart_probability),
                                    min_score=float(args.min_score))
        print(time.time() - t)

        if args.interactive:
            cls.run_interactive(node_index, edges)
        

    def run_interactive(node_index, edges):
        """Runs a flask web server with an interactive
        visualization of the obtained protein-protein
        interaction graph."""
        app = Flask(__name__)

        @app.route('/')
        def visualize():
            f = open('ppigepred/visualize.html')
            return f.read()

        @app.route('/graphjs')
        def graphjs():
            f = open('ppigepred/graph.js')
            return f.read()

        @app.route('/subgraph')
        def serve_graph():
            data = {
                'nodes': node_index,
                'edges': edges,
            }
            return json.dumps(data)

        print('Visualization available at http://localhost:5000/')
        app.run(debug=True, use_reloader=False)

