import argparse
import csv
import json
import logging
import multiprocessing

import networkx as nx
import pandas as pd

from flask import Flask
from flask import render_template

from .bioseqs import get_ensembl_release
from .bioseqs import gene_ids_to_protein_ids
from .bioseqs import protein_ids_to_gene_ids
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
    parser.add_argument('-g',
                        '--gene_ids',
                        action='store_true',
                        help='Use ENSEMBL gene ids instead of ENSEMBL protein ids.')
    parser.add_argument('-a',
                        '--assembly',
                        type=str,
                        default='GRCh38',
                        help='Allows to specify the assembly that will be used if --gene_ids flag is up.')

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
            references = cls.read_csv(args.references_file, header=None)
            cls.validate_sequence_list_file(references, args.references_file)
            references = set(references.iloc[:, 0])

        # If gene ids are used, convert the provided gene ids to
        # protein ids for internal use.
        if args.gene_ids:
            ensembl = get_ensembl_release(args.assembly)
            references = gene_ids_to_protein_ids(ensembl, references)

        candidates = None
        if args.candidates:
            candidates = set(args.candidates.split(','))
        elif args.candidates_file:
            candidates = cls.read_csv(args.candidates_file, header=None)
            cls.validate_sequence_list_file(candidates, args.candidates_file)
            candidates = set(candidates.iloc[:, 0])

        protein_interactions = cls.read_csv(args.db)
        cls.validate_protein_db_file(protein_interactions)
        protein_interactions.columns = ['protein1', 'protein2', 'score']
        if args.min_interaction_score:
            selected_rows = protein_interactions.score > args.min_interaction_score
            protein_interactions = protein_interactions[selected_rows]

        protein_graph = get_protein_graph(protein_interactions, references)

        probabilities = predict(protein_graph,
                                references,
                                walks=args.walks,
                                return_prob=args.restart_probability,
                                processes=args.processes)

        # Filter probabilities based on the user's input. It removes:
        # - nodes with low score
        # - nodes that are ranked too low
        # - nodes not included in the candidates set
        probabilities = filter_probabilities(probabilities,
                                             candidates,
                                             args.min_score,
                                             args.top)

        # Get the subgraph induced by the relevant nodes.
        graph = filter_graph(protein_graph, probabilities)

        # If gene ids are used, convert the internally
        # used protein ids to gene ids for output.
        if args.gene_ids:
            prot_to_gene_mapping = protein_ids_to_gene_ids(ensembl, probabilities)
            graph = nx.relabel_nodes(graph, prot_to_gene_mapping)
            probabilities = {prot_to_gene_mapping[prot_id]: prob
                             for prot_id, prob in probabilities.items()}

        # If output file is provided, write the results there.
        # Else, just print them to stdout.
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
            cls.run_interactive(node_data,
                                edges,
                                min_interaction_score=args.min_interaction_score)


    def read_csv(file, header='infer'):
        try:
            return pd.read_csv(file, header=header)
        except Exception as err:
            wrapper_err = RuntimeError("Cannot read file {}. Please provide a csv file.".format(file))
            raise wrapper_err from err


    def validate_protein_db_file(protein_interactions):
        """Validates the proper formatting of the
        provided protein interaction file."""
        col_types = protein_interactions.dtypes
        if (protein_interactions.shape[1] != 3 or
            not pd.api.types.is_string_dtype(col_types[0]) or
            not pd.api.types.is_string_dtype(col_types[1]) or
            not pd.api.types.is_numeric_dtype(col_types[2])):
            raise ValueError("Protein db file should contain three columns: " \
                             "(protein1:str, protein2:str, score:numeric)")


    def validate_sequence_list_file(seq_list, file):
        """Validates the proper formatting of the
        provided protein/gene sequence file."""
        if (seq_list.shape[1] != 1 or
            not pd.api.types.is_string_dtype(seq_list.iloc[:, 0].dtype)):
            raise ValueError("{} file should contain a single headerless column " \
                             "of protein or gene names.".format(file))


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
             score = graph.get_edge_data(*edge)['score']
             subgraph_data.append((node1_index, node2_index, score))

        return node_data, subgraph_data


    def run_interactive(node_data, edges, min_interaction_score=0):
        """Runs a flask web server with an interactive
        visualization of the obtained protein-protein
        interaction graph."""
        app = Flask(__name__)

        @app.route('/')
        def visualize():
            return render_template('visualize.html',
                                   min_interaction_score=min_interaction_score)

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

