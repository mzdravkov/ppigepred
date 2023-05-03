import argparse
import csv
import json
import logging
import multiprocessing

import pandas as pd

from flask import Flask

from .graphs import get_protein_graph
from .graphs import filter_graph
from .graphs import export_to_dot
from .prediction import predict
from .prediction import filter_probabilities


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
                        help='Probability of returning back to the initial node ' + \
                             'when doing the random walk simulation',
                        type=float,
                        default=0.2)
    parser.add_argument('-w',
                        '--walks',
                        help='Number of random walk simulations for each reference node',
                        type=int,
                        default=1000)
    parser.add_argument('-mis',
                        '--min-interaction-score',
                        help='Will ignore protein interactions with combined ' + \
                             'score below the specified value. ' + \
                             '(this filters the input network on which ' + \
                             'the random walk is performed)',
                        type=int,
                        default=10)
    parser.add_argument('-m',
                        '--min-score',
                        help='Minimum score to filter out insignificant results.',
                        type=float,
                        default=0.0005)
    parser.add_argument('-i',
                        '--interactive',
                        help='Run interactive graph visualization server.',
                        action='store_true')
    parser.add_argument('-p',
                        '--processes',
                        help='Number of processes that will be used.',
                        default=1,
                        type=int,
                        choices=range(1, multiprocessing.cpu_count()))
    parser.add_argument('-t',
                        '--top',
                        help='Get only the top N proteins',
                        type=int)
    parser.add_argument('-o',
                        '--output',
                        help='Output file')
    parser.add_argument('-gv',
                        '--graphviz',
                        help='Export the graph induced by the relevant ' + \
                             'nodes to a GraphViz dot file.')

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
        
        protein_graph = get_protein_graph(protein_interactions,
                                          references,
                                          candidates)

        probabilities = predict(protein_graph,
                                references,
                                walks=args.walks,
                                return_prob=args.restart_probability,
                                processes=args.processes)

        # Filter probabilities based on the user's input
        # (remove low score nodes or lower-ranked ones).
        probabilities = filter_probabilities(probabilities,
                                             args.min_score,
                                             args.top)
        
        # Get the subgraph induced by the relevant nodes.
        graph = filter_graph(protein_graph, probabilities)

        if args.output:
            with open(args.output, 'w') as file:
                writer = csv.writer(file)
                writer.writerow(('protein', 'probability'))
                writer.writerows(probabilities.items())
        else:
            for prot, prob in probabilities.items():
                print("{},{:.8f}".format(prot, prob))
                
        if args.graphviz:
            export_to_dot(graph, args.graphviz)

        if args.interactive:
            node_data, edges = cls.prepare_interactive_data(graph,
                                                            probabilities,
                                                            references)
            cls.run_interactive(node_data, edges)
                        

    def prepare_interactive_data(graph, probabilities, references):
        selected_node_index = {n: i for i, n in enumerate(graph.nodes)}

        node_data = {i: {
            'name': n,
            'density': probabilities[n],
            'is_reference': n in references,
            } for n, i in selected_node_index.items()}

        subgraph_data = []
        for edge in graph.edges():
             node1_index = selected_node_index[edge[0]]
             node2_index = selected_node_index[edge[1]]
             score = graph.get_edge_data(*edge)['combined_score']
             subgraph_data.append((node1_index, node2_index, score))

        return node_data, subgraph_data
    

    def run_interactive(node_data, edges):
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
                'nodes': node_data,
                'edges': edges,
            }
            return json.dumps(data)

        print('Visualization available at http://localhost:5000/')
        app.run(debug=True, use_reloader=False)

