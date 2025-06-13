import bz2
import gzip
import io
import logging
import re
from abc import ABC, abstractmethod
from collections import OrderedDict, namedtuple
from collections.abc import Iterable
from typing import Hashable, List, Optional, Self, TextIO

import Bio.SeqIO as SeqIO
import networkx as nx
import tqdm
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

SegmentInfo = namedtuple('SegmentInfo', ['name', 'length', 'depth'])


def update_segments(gfa_file: str, fasta_file: str, output_stream: io.BufferedWriter) -> None:
    """
    Update the segment sequences using the supplied fasta file. 

    :param gfa_file: Input gfa to update
    :param fasta_file: Input fasta for to use as sequence source
    :param output_stream: Output stream for an updated gfa file
    """

    gfa = GFA(gfa_file, skip_sequence_data=False, progress=True)

    changed = set()
    for new_seq in tqdm.tqdm(SeqIO.parse(fasta_file, 'fasta'), desc='Updating nodes'):

        if new_seq.id not in gfa.segments:
            raise KeyError('Fasta sequence {} was not found as a segment'.format(new_seq.id))

        si = gfa.segments[new_seq.id]
        si.length = len(new_seq)
        si.seq = str(new_seq.seq)
        logger.debug('Updated {}'.format(si.name))
        
        if si.name in changed:
            raise ValueError('Duplicate fasta entries for segment {}'.format(si.name))
        changed.add(si.name)

    logger.info('Updated {} of {} segments'.format(len(changed), len(gfa.segments)))

    gfa.write(output_stream)


def find_isolated_segments(gfa_file: str, min_len: int=2000, max_len: Optional[int]=None,
                           require_circular: bool=True, require_reciprocal: bool=False) -> List[SegmentInfo]:
    """
    Check a GFA file for putative circular, isolated segments. The test for circularity is simply
    the existence of a self-loop.

    Reciprocal edges is a strong constraint, leading to the removal of many edges. This may or may
    not lead to a substantial increase of false positives relative to any increase in true positives.

    :param gfa_file: The filename of the GFA.
    :param min_len: Minimum length below which segments are ignored.
    :param max_len: Maximum length above which segments are ignored.
    :param require_circular: Requires self-loop indicating circular.
    :param require_reciprocal: In conversion of the GFA to an undirected graph, require that there be
    reciprocal edges between u and v. i.e. (u,v) and (v,u). If not true, then no edge is added to
    the undirected graph.
    :return: A list of suspected circular isolated segments.
    """

    g = GFA(gfa_file, skip_sequence_data=True, progress=True)\
        .to_networkx(include_seq=False, progress=True)

    logger.debug(f'Graph details: order={g.order()}, size={g.size()}')

    if require_reciprocal:
        raise NotImplementedError('Checking for reciprocal edges is not currently implemented')

    suspects = []
    for i, gi in enumerate(nx.connected_components(g)):

        # single isolated segment
        if len(gi) != 1:
            continue

        u = gi.pop()

        # a simple test for circularity, the existence of a self-loop
        if require_circular and u not in g[u]:
            continue

        # skip very short or long segments
        len_u = g.nodes[u]['length']
        if len_u < min_len or (max_len is not None and len_u > max_len):
            continue

        dp_u = None
        if 'dp' in g.nodes[u]:
            dp_u = g.nodes[u]['dp']

        suspects.append(SegmentInfo(u, len_u, dp_u))

    return suspects


def open_input(filename: str) -> TextIO:
    """
    Open a text file for input. Only the ends of filenames are inspected to determine
    if the stream requires decompression, rather than reading the start of file.

    Recognises gzip and bz2.

    :param filename: The name of the input file.
    :return: An open file handle, possibly wrapped in a decompressor.
    """
    suffix = filename.split('.')[-1].lower()
    if suffix == 'bz2':
        return bz2.open(filename, 'rt')
    elif suffix == 'gz':
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'rt')


class MissingFieldsError(IOError):
    """
    Thrown records of fewer fields that expected by the spec.
    """
    def __init__(self, min_fields: int, sep: str, str_val: str) -> None:
        super(MissingFieldsError, self).__init__(
            'less than {} fields delimited by "{}" were found for record: {}'.format(
                min_fields, sep, str_val))


class BaseHelper(ABC):
    """
    A simple base, which simply removes the initial common task of parsing.
    """
    def __init__(self, min_fields: int, _type: Optional[str]=None, sep: str='\t') -> None:
        self._type = _type
        self.min_fields = min_fields
        self.sep = sep

    @abstractmethod
    def _id(self) -> Hashable:
        raise NotImplementedError('_id() has not been implemented')

    def __repr__(self) -> str:
        o = self._id()
        if isinstance(o, str):
            return o
        if isinstance(o, (int, float)):
            return str(o)
        elif isinstance(o, Iterable):
            return ','.join(str(v) for v in o)
        else:
            raise RuntimeError(f'Unknown type for __repr__: {o} which is type: {type(o)}')

    def __hash__(self) -> int:
        return hash(self._id())

    def __eq__(self, other: Self) -> bool:
        if not isinstance(other, self.__class__):
            return False
        return self._id() == other._id()

    def __ne__(self, other: Self) -> bool:
        return not self.__eq__(other)

    def split(self, str_val: str) -> List[str]:
        """
        Split the delimited string representation into fields.
        :param str_val: A string representation of a record
        :return: An array of fields
        """
        str_val = str_val.strip()
        if not str_val:
            raise IOError('string representation was empty or only whitespace')
        fields = str_val.split(self.sep)
        if len(fields) < self.min_fields:
            raise MissingFieldsError(self.min_fields,
                                     self.sep.encode(errors='backslashreplace').decode(),
                                     str_val)
        return fields

    def read_optionals(self, fields: List[str]) -> dict:
        """
        If a record has optional fields, beyond what are required, these will be read
        into a dictionary, keyed by their tags. Therefore, it is assumed that tags are
        unique within each record.
        :param fields: The array of fields for a given record
        :return: A dict of OptionalField objects
        """
        optionals = OrderedDict()
        if len(fields) > self.min_fields:
            for fi in fields[self.min_fields:]:
                o = OptionalField(fi)
                if o.tag in optionals:
                    raise KeyError(f'optional fields with duplicate tag {o.tag}')
                optionals[o.tag] = o
        return optionals


class OptionalField(BaseHelper):
    """
    A catch-all object for optional fields. As these fields declare their
    value type. This is taken into account when returning it with .get(), but
    only one type of int is returned.
    """
    def __init__(self, str_val: str) -> None:
        super(OptionalField, self).__init__(3, sep=':')

        tokens = self.split(str_val)

        # check that tag adheres to the specification for naming
        if not re.match('^[a-zA-Z][a-zA-Z0-9]$', tokens[0]):
            raise IOError('invalid field tag [{}]'.format(tokens[0]))

        self.tag = tokens[0]
        self.type = tokens[1]
        # we're assuming ":" can appear in val, therefore, we don't use the
        # last token split on ":"
        self.val = str_val[str_val.find(tokens[1]) + 2:]

        # set the function for returning a correctly typed value.
        # added all cases as some software does not respect specification.
        if self.type in 'AZJ':
            self.get = self.get_str
        elif self.type == 'i':
            self.get = self.get_int
        elif self.type == 'f':
            self.get = self.get_float
        elif self.type == 'H':
            self.get = self.get_bytearray
        elif self.type == 'B':
            self.get = self.get_numarray
        else:
            raise ValueError('unknown field type [{}] on line [{}]'.format(self.type, str_val))

    def __str__(self) -> str:
        if self.type in 'AZJif':
            return f'{self.tag}:{self.type}:{self.val}'
        else:
            raise NotImplementedError('Writing byte or num arrays (types H or B) is not implemented')

    def get_int(self) -> int:
        return int(self.val)

    def get_float(self) -> float:
        return float(self.val)

    def get_str(self) -> str:
        return self.val

    def get_bytearray(self) -> bytearray:
        return bytearray.fromhex(self.val)

    def get_numarray(self) -> List[int|float]:
        if self.val[0] in 'cCsSiI':
            return [int(v) for v in self.val[1:].split(',') if v]
        elif self.val[0] == 'f':
            return [float(v) for v in self.val[1:].split(',') if v]
        else:
            raise ValueError('unknown number array type [{}]'.format(self.val[0]))

    def _id(self) -> str:
        return self.tag


class Segment(BaseHelper):
    """
    A Segment record
    """
    def __init__(self, str_val: str, skip_seq: bool) -> None:
        super(Segment, self).__init__(3, _type='S')

        tokens = self.split(str_val)

        self.name = tokens[1]
        self.seq = None
        self.length = 0
        # '*' marks a dummy sequence
        if tokens[2] != '*':
            if not skip_seq:
                self.seq = tokens[2]
            self.length = len(tokens[2])

        self.optionals = self.read_optionals(tokens)

    def _id(self) -> str:
        return self.name


class Link(BaseHelper):
    """
    A Link record.
    """
    def __init__(self, str_val: str) -> None:
        super(Link, self).__init__(6, _type='L')

        tokens = self.split(str_val)

        self.src = tokens[1]
        self.src_orient = tokens[2]
        self.dest = tokens[3]
        self.dest_orient = tokens[4]

        self.overlap = None
        # overlaps can have dummy entries
        if tokens[5] != '*':
            self.overlap = tokens[5]

        self.optionals = self.read_optionals(tokens)

    def _id(self) -> tuple[str, str, str, str]:
        return self.src, self.src_orient, self.dest, self.dest_orient


class Containment(BaseHelper):
    """
    A Containment record
    """
    def __init__(self, str_val: str) -> None:
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

    def _id(self) -> tuple[str, str, str, str]:
        return self.container, self.container_orient, self.contained, self.contained_orient


class Path(BaseHelper):
    """
    A Path record.
    """
    def __init__(self, str_val: str) -> None:
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

    def _id(self) -> str:
        return self.name


class GFA(object):

    def to_networkx(self, include_seq: bool=False, annotate_paths: bool=False, collections_to_str: bool=True,
                    add_label: bool=False, progress: bool=False) -> nx.MultiGraph:
        """
        Convert the instance of a Networkx DiGraph.
        :param include_seq: Include the segment sequences as node attributes
        :param annotate_paths: Add path membership as node attributes
        :param collections_to_str: Convert attributes which are collections to strings.
        :param add_label: Include segment names also as node attributes (label)
        :param progress: Show progress
        :return: A networkx.MultiGraph
        """
        g = nx.MultiGraph()

        # Add nodes. We do this first to include unconnected nodes
        for si in tqdm.tqdm(self.segments.values(), desc='Creating nodes', disable=not progress):
            attrs = {'length': si.length}
            for k, v in si.optionals.items():
                attrs[k] = v.get()
            if include_seq:
                attrs['seq'] = si.seq
            if add_label:
                attrs['label'] = si.name
            g.add_node(si.name, **attrs)

        # Add edges
        for li in tqdm.tqdm(self.links.values(), desc='Creating edges', disable=not progress):
            attrs = {
                'src_orient': li.src_orient,
                'dest_orient': li.dest_orient,
                'ovlp': li.overlap
            }
            g.add_edge(li.src, li.dest, **attrs)

        if annotate_paths:
            # Collect all the path steps for each node
            for n, p in enumerate(self.paths):
                last = len(p.segment_names)
                for m, (sid, so) in enumerate(p.segment_names):
                    if m == last-1:
                        # mark the last element so it's easy to discern
                        g.nodes[sid].setdefault('paths', []).append('p{}.{}|'.format(n, m))
                    else:
                        g.nodes[sid].setdefault('paths', []).append('p{}.{}'.format(n, m))

        if collections_to_str:
            # Stringify collection attributes.
            # For instance, networkx graphml does not support collections.
            for n, d in g.nodes(data=True):
                for k, v in d.items():
                    if not isinstance(v, str) and isinstance(v, Iterable):
                        d[k] = ' '.join(vv for vv in sorted(v))

        return g

    def to_fasta(self, output: str | io.BufferedWriter) -> None:
        """
        Write all segments to Fasta
        :param output: Output file name or handle
        """
        if isinstance(output, str):
            out_handle = open(output, 'wt')
        else:
            out_handle = output

        try:
            num_seq = 0
            for _id, _s in self.segments.items():
                desc = ' '.join(['{}:{}'.format(_o._id(), _o.get()) for _o in _s.optionals.values()])
                SeqIO.write(SeqRecord(Seq(_s.seq), id=_id, description=desc), out_handle, 'fasta')
                num_seq += 1
            logger.info('Wrote {} segments'.format(num_seq))

        finally:
            if isinstance(output, str):
                out_handle.close()

    def __repr__(self) -> str:
        return 'Segments: {}, Paths: {}, Links: {}, Containments: {}'.format(
            len(self.segments), len(self.paths), len(self.links), len(self.containments))

    def __str__(self) -> str:
        return self.__repr__()

    def write(self, output: str|io.BufferedWriter) -> None:

        if isinstance(output, str):
            out_handle = open(output, 'wt')
        else:
            out_handle = output

        logger.debug('Writing GFA to {}'.format(out_handle.name))

        # write comments
        for ci in self.comments:
            out_handle.write('#{}\n'.format(ci))

        # write header
        logger.debug(self.header)
        if 'VN' in self.header:
            out_handle.write('{}\n'.format('\t'.join(['H', str(self.header['VN'])])))

        def to_str(f:str|int|float|List) -> str:
            if not f:
                return '*'
            else:
                if isinstance(f, str):
                    return f
                elif isinstance(f, list):
                    return ','.join(str(fi) for fi in f)
                else:
                    return str(f)

        # write segments
        for si in self.segments.values():

            fields = ['S', si.name, to_str(si.seq)]
            for oi in si.optionals.values():
                fields.append(str(oi))
            out_handle.write('{}\n'.format('\t'.join(fields)))

        # write links
        for li in self.links.values():

            fields = ['L',
                      li.src, li.src_orient,
                      li.dest, li.dest_orient,
                      to_str(li.overlap)]
            for oi in li.optionals.values():
                fields.append(str(oi))
            out_handle.write('{}\n'.format('\t'.join(fields)))

        # write containments
        for ci in self.containments.values():

            fields = ['C',
                      ci.container, ci.container_orient,
                      ci.contained, ci.contained_orient,
                      to_str(ci.pos), to_str(ci.overlap)]
            for oi in ci.optionals.values():
                fields.append(str(oi))
            out_handle.write('{}\n'.format('\t'.join(fields)))

        # write paths
        for pi in self.paths.values():

            segment_names_str = ','.join('{}{}'.format(sni[0], sni[1]) for sni in pi.segment_names)

            if not pi.overlaps:
                overlap_str = '*'
            else:
                overlap_str = ','.join(ovi for ovi in pi.overlaps)

            fields = ['P', pi.name, segment_names_str, overlap_str]
            for oi in pi.optionals.values():
                fields.append(str(oi))
            out_handle.write('{}\n'.format('\t'.join(fields)))

    def _add_version(self, ver: float=1.0) -> None:
        logger.debug(f'Adding version to header: {ver}')
        o = OptionalField(f'VN:Z:{ver}')
        self.header[o.tag] = o

    def _validate(self) -> None:
        """
        Apply some simple validation rules
        """
        if 'VN' not in self.header:
            logger.warning('No version specified in header. Assuming 1.0')
            self._add_version()
        elif len(self.header) == 0:
            logger.warning('No header found in GFA')
            self._add_version()

    def __init__(self, filename: str, ignore_isolate_paths: bool=False,
                 skip_sequence_data: bool=False, progress: bool=False) -> None:
        """
        Instantiate from a file. Only a single GFA record is expected to be found.
        :param filename: The name of gfa file to read
        :param ignore_isolate_paths: Some assemblers specify non-compliant, single node paths.
        :param progress: Show progress while reading
        """
        self.filename = filename
        self.ignore_isolate_paths = ignore_isolate_paths
        self.skip_sequence_data = skip_sequence_data
        self.comments = []
        self.header = OrderedDict()
        self.segments = OrderedDict()
        self.links = OrderedDict()
        self.containments = OrderedDict()
        self.paths = OrderedDict()
        self.version = '{unset}'

        with open_input(filename) as input_h:

            for line in tqdm.tqdm(input_h, desc='Reading GFA', disable=not progress):

                line = line.strip()
                if not line:
                    break

                if line.startswith('#'):
                    # a comment that be could possibly read in
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

                    # check the version? Why bother
                    if found_ver:
                        if not re.match(r'^1[\\.0-9]+', self.version):
                            raise IOError('unsupported GFA version [{}]'.format(self.version))

                elif line.startswith('S'):
                    obj = Segment(line, self.skip_sequence_data)
                    self.segments[obj.name] = obj

                elif line.startswith('L'):
                    obj = Link(line)
                    self.links[obj] = obj

                elif line.startswith('C'):
                    obj = Containment(line)
                    self.containments[obj] = obj

                elif line.startswith('P'):
                    try:
                        obj = Path(line)
                        self.paths[obj] = obj
                    except MissingFieldsError as e:
                        if not ignore_isolate_paths:
                            raise e

        self._validate()