# "python src/philharmonic_cytoscape.py -s {input.styles} -o {output} {input.clusters}

import argparse
from pathlib import Path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a Cytoscape style file from a cluster file')
    parser.add_argument('-s', '--styles', type=str, help='Cytoscape style file')
    parser.add_argument('-o', '--output', type=str, help='Output file')
    parser.add_argument('clusters', type=str, help='Cluster file')

    args = parser.parse_args()
    Path(args.output).touch()