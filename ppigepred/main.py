import argparse
import logging

logging.basicConfig(filename='logs.log', level=logging.DEBUG)

parser = argparse.ArgumentParser(description="Predicts genes related to disease based on the protein-protein interactions they participate in.")

parser.add_argument('--db', help='Protein-protein interaction db. CSV with PROT1,PROT2,SCORE')
parser.add_argument("-r", "--ref", help='List of disease-related gene symbols in comma-separated form')
parser.add_argument("-rf", "--ref-file", help='File that contains one disease-related gene symbol per line')
parser.add_argument("-c", "--candidate", help='List of candidate gene symbols in comma-separated form')
parser.add_argument("-cf", "--candidate-file", help='File that contains one candidate-gene symbol per line')


if __name__ == '__main__':
    args = parser.parse_args()

    #if args.subcommand == 'align':
    #    align(args)
