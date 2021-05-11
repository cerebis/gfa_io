from gfa_io import gfa_io
import argparse
import logging
import pandas

logger = logging.getLogger(__name__)

def main():

    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')

    parser = argparse.ArgumentParser('GFA utils', parents=[global_parser], add_help=True)
    subparsers = parser.add_subparsers(title='commands', dest='command', description='Valid commands',
                                       help='choose an analysis stage for further options')
    subparsers.required = True

    cmd_isolates = subparsers.add_parser('isolates',
                                         description='Find isolated segments in GFA file')
    cmd_isolates.add_argument('--reciprocal', default=False, action='store_true',
                              help='Require reciprocal edges between segments')
    cmd_isolates.add_argument('--min-len', type=int, default=2000, help='Minimum segment length')
    cmd_isolates.add_argument('--max-len', type=int, default=None, help='Maximum segment length')
    cmd_isolates.add_argument('--circular', default=False, action='store_true', help='Maximum segment length')
    cmd_isolates.add_argument('GFA_INPUT', help='Input GFA file')
    cmd_isolates.add_argument('CSV_OUTPUT', help='Output CSV file')

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

    if args.command == 'isolates':

        suspects = gfa_io.find_isolated_segments(args.GFA_INPUT, args.min_len, args.max_len,
                                                 args.circular, args.reciprocal)
        suspects = pandas.DataFrame(suspects).sort_values('length', ascending=False)
        logger.info('Found {} segments'.format(len(suspects)))
        suspects.to_csv(args.CSV_OUTPUT, index=False)
