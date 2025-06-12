import argparse
import logging
import sys

import networkx as nx
import pandas

import gfa_io

logger = logging.getLogger(__name__)


def main() -> None:

    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true',
                               help='Verbose output')

    parser = argparse.ArgumentParser('GFA utils', parents=[global_parser], add_help=True)
    subparsers = parser.add_subparsers(title='commands', dest='command', description='Valid commands',
                                       help='choose an analysis stage for further options')
    subparsers.required = True

    cmd_update = subparsers.add_parser('update-segments', 
                                       description='Update segment sequences from FASTA')
    cmd_update.add_argument('GFA_INPUT', help='Input GFA file')
    cmd_update.add_argument('FASTA_INPUT',
                            help='Input FASTA file containing updated segment sequences (ie. post-polishing')
    cmd_update.add_argument('GFA_OUTPUT', help='Output GFA file [stdout]', nargs='?', default=sys.stdout,
                            type=argparse.FileType('wt'))

    cmd_isolates = subparsers.add_parser('isolates',
                                         description='Find isolated segments in GFA file')
    cmd_isolates.add_argument('--reciprocal', default=False, action='store_true',
                              help='Require reciprocal edges between segments')
    cmd_isolates.add_argument('--min-len', type=int, default=2000, help='Minimum segment length')
    cmd_isolates.add_argument('--max-len', type=int, default=None, help='Maximum segment length')
    cmd_isolates.add_argument('--circular', default=False, action='store_true', help='Maximum segment length')
    cmd_isolates.add_argument('GFA_INPUT', help='Input GFA file')
    cmd_isolates.add_argument('CSV_OUTPUT', help='Output CSV file [stdout]', nargs='?', default=sys.stdout,
                              type=argparse.FileType('wt'))

    cmd_segments = subparsers.add_parser('dump-segments',
                                         description='Dump segments to FASTA')
    cmd_segments.add_argument('GFA_INPUT', help='Input GFA file')
    cmd_segments.add_argument('FASTA_OUTPUT', help='Output FASTA file [stdout]', nargs='?', default=sys.stdout,
                              type=argparse.FileType('wt'))

    cmd_convert = subparsers.add_parser('convert',
                                        description='Convert GFA to general graph formats')
    cmd_convert.add_argument('-f', '--format', choices=['graphml', 'gml', 'edgelist'],
                             default='graphml', help='Graph format')
    cmd_convert.add_argument('--node-labels', default=False, action='store_true',
                             help='Add additional label attributes to nodes')
    cmd_convert.add_argument('--seqs', default=False, action='store_true', help='Include sequence data')
    cmd_convert.add_argument('--paths', default=False, action='store_true', help='Annotate paths')
    cmd_convert.add_argument('--ignore', default=False, action='store_true', help='Ignore empty paths')
    cmd_convert.add_argument('GFA_INPUT', help='Input GFA file')
    cmd_convert.add_argument('GRAPH_OUTPUT', help='Output graph file')

    args = parser.parse_args()

    logging.captureWarnings(True)
    logger = logging.getLogger('main')

    # root log listens to everything
    root = logging.getLogger('')
    root.setLevel(logging.DEBUG)

    # log message format
    formatter = logging.Formatter(fmt='%(levelname)-8s | %(asctime)s | %(name)7s | %(message)s')

    # Runtime console listens to INFO by default
    ch = logging.StreamHandler()
    if args.verbose:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    if args.command == 'update-segments':
        gfa_io.update_segments(args.GFA_INPUT, args.FASTA_INPUT, args.GFA_OUTPUT)

    elif args.command == 'isolates':

        suspects = gfa_io.find_isolated_segments(args.GFA_INPUT, args.min_len, args.max_len,
                                                 args.circular, args.reciprocal)
        if len(suspects) <= 0:
            logger.info('No isolated segments found')
        else:
            suspects = pandas.DataFrame(suspects).sort_values('length', ascending=False)
            logger.info('Found {} segments'.format(len(suspects)))
            suspects.to_csv(args.CSV_OUTPUT, index=False)

    elif args.command == 'dump-segments':

        gfa_io.GFA(args.GFA_INPUT).to_fasta(args.FASTA_OUTPUT)

    elif args.command == 'convert':

        logger.debug('Reading GFA')
        gfa = gfa_io.GFA(args.GFA_INPUT, args.ignore)

        logger.debug('Converting to GraphML')
        g = gfa.to_networkx(include_seq=args.seqs, annotate_paths=args.paths,
                            collections_to_str=True, add_label=args.node_labels)
        logger.debug(f'Graph size: {g.size()}, order: {g.order()}')

        write_fmt = {
            'graphml': nx.write_graphml,
            'gml': nx.write_gml,
            'edgelist': nx.write_edgelist
        }

        try:
            write_fmt[args.format](G=g, path=args.GRAPH_OUTPUT)
        except KeyError:
            raise RuntimeError('Unsupported format {}'.format(args.format))
