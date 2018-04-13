from collections import Iterable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO as SeqIO
import networkx as nx
import re
import bz2
import gzip


def open_input(fname):
    """
    Open a text file for input. Only the ends of filenames are inspected to determine
    if the stream requires decompression, rather than reading the start of file.

    Recognises gzip and bz2.

    :param fname: the name of the input file
    :return: open file handle, possibly wrapped in a decompressor
    """
    suffix = fname.split('.')[-1].lower()
    if suffix == 'bz2':
        return bz2.BZ2File(fname, 'r')
    elif suffix == 'gz':
        return gzip.GzipFile(fname, 'r')
    else:
        return open(fname, 'r')


class MissingFieldsError(IOError):
    """
    Thrown records of fewer fields that expected by the spec.
    """
    def __init__(self, min_fields, sep, str_val):
        super(MissingFieldsError, self).__init__(
            'less than {} fields delimited by "{}" were found for record: {}'.format(
                min_fields, sep, str_val))


class BaseHelper(object):
    """
    A simple base, which simply removes the initial common task of parsing.
    """
    def __init__(self, min_fields, _type=None, sep='\t'):
        self._type = _type
        self.min_fields = min_fields
        self.sep = sep

    def _id(self):
        raise RuntimeError('_id() has not been implemented')

    def __repr__(self):
        o = self._id()
        if isinstance(o, (int, float, basestring)):
            return o
        elif isinstance(o, Iterable):
            return ','.join(str(v) for v in o)
        return o

    def __hash__(self):
        return hash(self._id())

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self._id() == other._id()

    def __ne__(self, other):
        return not self.__eq__(other)

    def split(self, str_val):
        """
        Split the delimited string representation into fields.
        :param str_val: string representation of a record
        :return: array of fields
        """
        str_val = str_val.strip()
        if not str_val:
            raise IOError('string representation was empty or only whitespace')
        fields = str_val.split(self.sep)
        if len(fields) < self.min_fields:
            raise MissingFieldsError(self.min_fields, self.sep.encode('string_escape'), str_val)
        return fields

    def read_optionals(self, fields):
        """
        If a record has optional fields, beyond what are required, these will be read
        into a dictionary, keyed by their tags. Therefore, it is assumed that tags are
        unique within each record.
        :param fields: the array of fields for a given record
        :return: dict of OptionalField objects
        """
        optionals = {}
        if len(fields) > self.min_fields:
            for fi in fields[self.min_fields:]:
                o = OptionalField(fi)
                assert o not in optionals, 'optional fields with duplicate tag {}'.format(o.tag)
                optionals[o.tag] = o
        return optionals


class OptionalField(BaseHelper):
    """
    A catch all object for optional fields. As these fields declare their
    value type, this is taken into account when returning it with .get(),
    however only one type of int is returned.
    """
    def __init__(self, str_val):
        super(OptionalField, self).__init__(3, sep=':')

        tokens = self.split(str_val)

        # check that tag adheres to the specification for naming
        if not re.match('^[a-zA-Z][a-zA-Z0-9]$', tokens[0]):
            raise IOError('invalid field tag [{}]'.format(tokens[0]))

        self.tag = tokens[0]
        self.type = tokens[1]
        # we're assuming ":" can appear in value, therefore we don't use the
        # last token split on ":"
        self.val = str_val[str_val.find(tokens[1]) + 2:]

        # set the function for returning a correctly typed value.
        # added all cases as some software does not respect specification.
        if self.type in 'AZJazj':
            self.get = self.get_str
        elif self.type in 'iI':
            self.get = self.get_int
        elif self.type in 'fF':
            self.get = self.get_float
        elif self.type in 'Hh':
            self.get = self.get_bytearray
        elif self.type in 'Bb':
            self.get = self.get_numarray
        else:
            raise IOError('unknown field type [{}] on line [{}]'.format(self.type, str_val))

    def get_int(self):
        return int(self.val)

    def get_float(self):
        return float(self.val)

    def get_str(self):
        return self.val

    def get_bytearray(self):
        return bytearray.fromhex(self.val)

    def get_numarray(self):
        if self.val[0] in 'cCsSiI':
            return [int(v) for v in self.val[1:].split(',')]
        elif self.val[0] == 'f':
            return [float(v) for v in self.val[1:].split(',')]
        else:
            raise IOError('unknown number array type [{}]'.format(self.val[0]))

    def _id(self):
        return self.tag


class Segment(BaseHelper):
    """
    A Segment record
    """
    def __init__(self, str_val):
        super(Segment, self).__init__(3, _type='S')

        tokens = self.split(str_val)

        self.name = tokens[1]
        self.seq = None
        # '*' marks a dummy sequence
        if tokens[2] != '*':
            self.seq = tokens[2]

        self.optionals = self.read_optionals(tokens)

    def length(self):
        """
        :return: sequence length
        """
        return len(self.seq)

    def _id(self):
        return self.name


class Link(BaseHelper):
    """
    A Link record.
    """
    def __init__(self, str_val):
        super(Link, self).__init__(6, _type='L')

        tokens = self.split(str_val)

        self._from = tokens[1]
        self.from_orient = tokens[2]
        self.to = tokens[3]
        self.to_orient = tokens[4]

        self.overlap = None
        # overlaps can have dummy entries
        if tokens[5] != '*':
            self.overlap = tokens[5]

        self.optionals = self.read_optionals(tokens)

    def _id(self):
        return self._from, self.from_orient, self.to, self.to_orient


class Containment(BaseHelper):
    """
    A Containment record
    """
    def __init__(self, str_val):
        super(Containment, self).__init__(7, _type='C')

        tokens = self.split(str_val)

        self.container = tokens[1]
        self.container_orient = tokens[2]
        self.contained = tokens[3]
        self.contained_orient = tokens[4]
        self.pos = int(tokens[5])

        self.overlap = None
        # overlap can have dummy entries
        if tokens[6] != '*':
            self.overlap = tokens[6]

        self.optionals = self.read_optionals(tokens)

    def _id(self):
        return self.container, self.container_orient, self.contained, self.contained_orient


class Path(BaseHelper):
    """
    A Path record.
    """
    def __init__(self, str_val):
        super(Path, self).__init__(4, _type='P')

        tokens = self.split(str_val)

        self.name = tokens[1]

        self.segment_names = []
        for ti in tokens[2].split(','):
            sid, so = ti[:-1], ti[-1]
            self.segment_names.append((sid, so))
        self.overlaps = []

        # overlap can have dummy entries
        if tokens[3] != '*':
            for ti in tokens[3].split(','):
                self.overlaps.append(ti)

    def _id(self):
        return self.name


class GFA(object):

    def to_graph(self, include_seq=False, annotate_paths=False, collections_to_str=True, add_label=False):
        """
        Convert the instance of a Networkx DiGraph.
        :param include_seq: include the sequencce as a node attribute
        :param annotate_paths: add node path membership as node attribute
        :param collections_to_str: convert attributes which are collections to strings.
        :return: networkx.DiGraph
        """
        g = nx.DiGraph()

        # add nodes. We do this first to include unconnected nodes
        for si in self.segments.values():
            attrs = {'length': si.length()}
            for k, v in si.optionals.iteritems():
                attrs[k] = v.get()
            if include_seq:
                attrs['seq'] = si.seq
            if add_label:
                attrs['label'] = si.name
            g.add_node(si.name, attrs)

        # add edges
        for li in self.links.values():
            attrs = {
                'udir': True if li.from_orient == '+' else False,
                'vdir': True if li.to_orient == '+' else False,
                'ovlp': li.overlap
            }
            g.add_edge(li._from, li.to, attrs)

        if annotate_paths:
            # g.graph['paths'] = {}
            # collect all the path steps for each node
            for n, p in enumerate(gfa.paths):
                # g.graph['paths'][p.name] = p.segment_names
                last = len(p.segment_names)
                for m, (sid, so) in enumerate(p.segment_names):
                    if m == last-1:
                        # mark the last element so its easy to discern
                        g.node[sid].setdefault('paths', []).append('p{}.{}|'.format(n, m))
                    else:
                        g.node[sid].setdefault('paths', []).append('p{}.{}'.format(n, m))

        if collections_to_str:
            # stringify collection attributes.
            # as for instance networkx graphml does not support collections.
            for n, d in g.nodes_iter(data=True):
                for k, v in d.iteritems():
                    if not isinstance(v, basestring) and isinstance(v, Iterable):
                        d[k] = ' '.join(vv for vv in sorted(v))

        return g

    def to_fasta(self, output, verbose=False):
        """
        Write all segments to Fasta
        :param output: output file name or handle
        :param verbose: enable some runtime information
        """
        if isinstance(output, basestring):
            out_hndl = open(output, 'w')
        else:
            out_hndl = output

        try:
            n = 0
            for _id, _s in self.segments.iteritems():
                SeqIO.write(SeqRecord(Seq(_s.seq), id=_id, description=''), out_hndl, 'fasta')
                n += 1
            if verbose:
                print 'Wrote {} segments'.format(n)

        finally:
            if isinstance(output, basestring):
                out_hndl.close()

    def __repr__(self):
        return 'Segments: {}, Paths: {}, Links: {}, Containments: {}'.format(
            len(self.segments), len(self.paths), len(self.links), len(self.containments))

    def __str__(self):
        return self.__repr__()

    def __init__(self, filename, ignore_isolate_paths=False):
        """
        Instantiate from a file. Only a single GFA record is expected to be found.
        :param filename: the gfa file to read
        :param ignore_isolate_paths: some assemblers specify non-compliant, single node paths.
        """
        self.filename = filename
        self.ignore_isolate_paths = ignore_isolate_paths
        self.comments = []
        self.header = {}
        self.segments = {}
        self.links = {}
        self.containments = []
        self.paths = []

        with open_input(filename) as input_h:

            for line in input_h:

                line = line.strip()
                if not line:
                    break

                if line.startswith('#'):
                    # a comment, could possibly read in
                    self.comments.append(line[1:])

                elif line.startswith('H'):
                    # header
                    tokens = line.split('\t')
                    if len(tokens) < 2:
                        raise IOError('Unspecified version number')

                    found_ver = False
                    for ti in tokens[1:]:
                        o = OptionalField(ti)
                        self.header[o.tag] = o
                        if o.tag == 'VN':
                            found_ver = True
                            self.version = o.get()

                    # check version? Why bother
                    if found_ver:
                        if not re.match('^1[\.0-9]+', self.version):
                            raise IOError('unsupported GFA version [{}]'.format(self.version))

                elif line.startswith('S'):
                    s = Segment(line)
                    self.segments[s.name] = s

                elif line.startswith('L'):
                    l = Link(line)
                    self.links[l] = l

                elif line.startswith('C'):
                    self.containments.append(Containment(line))

                elif line.startswith('P'):
                    try:
                        self.paths.append(Path(line))
                    except MissingFieldsError as e:
                        if not ignore_isolate_paths:
                            raise e


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--node-labels', default=False, action='store_true', help='Add additional label attributes to nodes')
    parser.add_argument('--seqs', default=False, action='store_true', help='Include sequence data')
    parser.add_argument('--paths', default=False, action='store_true', help='Annotate paths')
    parser.add_argument('--ignore', default=False, action='store_true', help='Ignore empty paths')
    parser.add_argument('--dedup', default=False, action='store_true', help='Remove what appear to be redundant paths')
    parser.add_argument('input', help='Input GFA file')
    parser.add_argument('output', help='Output GraphML file')
    args = parser.parse_args()

    print 'Reading GFA'
    gfa = GFA(args.input, args.ignore)

    print 'Converting to GraphML'
    g = gfa.to_graph(include_seq=args.seqs, annotate_paths=args.paths, collections_to_str=True,
                     add_label=args.node_labels)
    print nx.info(g)

    if args.dedup:
        print '\n\nRemoving redundant edges'
        to_del = set()
        for u, v in g.edges_iter():
            if g.has_edge(v, u):
                if g[u][v]['udir'] == g[v][u]['udir'] and g[u][v]['vdir'] == g[v][u]['vdir']:
                    if v < u:
                        u, v = v, u
                    to_del.add((u, v))
        g.remove_edges_from(to_del)
        print nx.info(g)

    nx.write_graphml(g, args.output)
    # nx.write_yaml(g, args.output)
