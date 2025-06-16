import bz2
import gzip
import logging
import re
import typing
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


def update_segments(gfa_file: str, fasta_file: str, output_stream: typing.TextIO) -> None:
    """
    Updates segments in a GFA file with sequences from a FASTA file.

    This function reads a GFA file and updates the segment sequences and their
    lengths based on the corresponding sequences in the provided FASTA file.
    The GFA file is parsed, segments are updated if matches are found in the
    FASTA file, and the updated GFA data is written to the output stream.
    Segments in the FASTA file that are not found in the GFA file result in
    a KeyError. If duplicate segment entries are found in the FASTA file, a
    ValueError is raised. Progress of the updates is displayed, and logs are
    generated for each updated segment and summary.

    Parameters:
        gfa_file: str
            Path to the input GFA file.
        fasta_file: str
            Path to the input FASTA file containing updated sequences.
        output_stream: TextIO
            Stream to write the updated GFA content.

    Raises:
        KeyError:
            If a sequence ID in the FASTA file is not found as a segment in
            the GFA file.
        ValueError:
            If duplicate segment entries are present in the FASTA file.

    Returns:
        None
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

    Reciprocal edges is an additional constraint that impacts the strength of filtering and affects
    the identification of suspected segments. If reciprocal edges are required, the method excludes
    non-reciprocal edges during graph construction.

    Parameters:
    gfa_file: Filename of the GFA file to analyze.
    min_len: Minimum segment length to consider, below which segments are ignored.
    max_len: Maximum segment length to consider, above which segments are ignored.
    require_circular: Boolean indicator for requiring a self-loop as evidence of circularity.
    require_reciprocal: Boolean indicator for requiring reciprocal edges in the graph.

    Returns:
    List of SegmentInfo objects describing the suspected circular isolated segments.
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
    Opens a file based on its compression type for reading as text. Supports
    plain text files, as well as files compressed in bz2 and gzip formats.

    Parameters:
        filename: str
            The path to the file to be opened. File extension determines
            if decompression is needed.

    Returns:
        TextIO
            A file object opened in text mode, allowing for reading operations.

    Raises:
        FileNotFoundError
            If the file specified by the filename does not exist.
        OSError
            If there is an issue opening the file or handling the compression format.
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
    Indicates that a record contains fewer fields than expected by the specification.

    This exception is used to notify situations where the number of
    delimited fields in a record does not meet the minimum required count
    specified in the provided arguments. It encapsulates both the expected
    number of fields and the problematic record for diagnostic purposes.

    Attributes:
        min_fields: int
            The minimum number of fields expected in the record.
        sep: str
            The string delimiter used to split the record into fields.
        str_val: str
            The problematic record that caused the exception.
    """
    def __init__(self, min_fields: int, sep: str, str_val: str) -> None:
        super(MissingFieldsError, self).__init__(
            'less than {} fields delimited by "{}" were found for record: {}'.format(
                min_fields, sep, str_val))
        self.min_fields = min_fields
        self.sep = sep
        self.str_val = str_val


class BaseHelper(ABC):
    """
    An abstract base class for handling common parsing tasks.

    This class provides basic functionality to manage delimited string
    records, ensure a minimum number of required fields, and handle
    optional fields. It serves as a base for implementing more specific
    parsing logic in subclasses.

    Attributes
    ----------
    _type : Optional[str]
        An optional field to specify the type of the record.
    min_fields : int
        The minimum number of required fields in the record.
    sep : str
        The delimiter used to split fields in the record string.
    """
    def __init__(self, min_fields: int, _type: Optional[str]=None, sep: str='\t') -> None:
        self._type = _type
        self.min_fields = min_fields
        self.sep = sep

    @abstractmethod
    def _id(self) -> Hashable:
        """
        Defines an abstract method `_id`. This method must be implemented in any subclass.
        It enforces providing a unique identifier for the implementing class.

        Methods:
            _id: Abstract method that should return a hashable value unique to the implementing class.

        Raises:
            NotImplementedError: If the `_id` method is not implemented in the subclass.

        Returns:
            Hashable: A unique identifier specific to an instance of the implementing class.
        """
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
        Split the delimited string into fields while performing basic validation
        of the string format and content.

        Checks include ensuring the string is not empty or contains only whitespace,
        and verifying that the resulting fields meet the required minimum count.

        :param str_val: A string representation of a record.
        :type str_val: str

        :raises IOError: If the input string is empty or contains only whitespace.
        :raises MissingFieldsError: If the number of fields after splitting is less than
            the minimum required.

        :return: A list of fields obtained by splitting the input string based on the separator.
        :rtype: List[str]
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
        If a record contains additional optional fields beyond the required ones, these fields
        are processed and returned as a dictionary, where each entry is keyed by its unique tag.

        Parameters:
        fields (List[str]): A list of strings representing the fields of a record.

        Returns:
        dict: A dictionary where each key corresponds to the unique tag of an optional field
              and the value is an OptionalField object.

        Raises:
        KeyError: If there are multiple optional fields with the same tag in the given record.
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
    Represents a specialized field object for optional data fields tagged with a specific format.

    The OptionalField class is intended as a utility to parse, handle, and retrieve values associated
    with optional fields in a specified format. It validates the tag format, identifies the type of
    field, and provides the means to extract data in the appropriate type. Multiple field value types
    are supported, including string, integer, float, bytearray, and numeric arrays. Instances of this
    class also allow easy conversion back to the original string representation.

    Attributes:
        tag (str): The identifier for the field, adhering to naming specification rules.
        type (str): The single-character type descriptor of the field.
        val (str): The value part of the field as a raw string.

    Raises:
        IOError: If the field tag fails to meet the specified naming rules.
        ValueError: If an unsupported or unknown field type is encountered during initialization.
        NotImplementedError: If attempting to convert byte arrays or numeric arrays (types H or B)
            back to string representation.
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
    Represents a Segment record (S) in GFA1 and GFA2.

    This class is used to manage and store information about a biological
    sequence segment. It provides functionality to parse a string
    representation of the segment, extracting its name, sequence, and
    metadata, and optionally skipping sequence extraction.

    Attributes:
        name (str): The name of the segment extracted from the input string.
        seq (Optional[str]): The sequence of the segment, which may be None
            if skipping sequence extraction or if the sequence is marked as
            a dummy sequence ('*').
        length (int): The length of the sequence. It is set to 0 if the
            sequence is not provided or is marked as dummy.
        optionals (dict): A dictionary containing optional attributes
            extracted from the input string.
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
    Represents a Link record.

    This class is used to handle and process a GFA Link (L) record, performing
    initialization and providing a mechanism to retrieve an identifier tuple
    that uniquely represents the link. It includes fields for source, destination,
    orientations, overlaps, and optional attributes.
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
    A Containment record (C) in GFA2.

    Represents a containment record in the GFA2 format, encapsulating details
    about the contained element and its relationship with the container element.
    Provides functionality for initialization, parsing, and extracting optional
    fields and identifiers.

    Attributes:
        container: The identifier of the container element.
        container_orient: Orientation of the container element.
        contained: The identifier of the contained element.
        contained_orient: Orientation of the contained element.
        pos: An integer representing the position of the containment.
        overlap: String representing the overlap information unless set to
            a dummy entry ('*').
        optionals: A dictionary storing any additional optional fields
            from the record.
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
    A representation of a Path record (P) in the GFA format.

    This class is designed to parse and represent the information contained
    within a Path record in the GFA (Graphical Fragment Assembly) format.
    A Path record denotes a sequence of segments within a graph, including
    their orientations (forward or reverse) and overlaps between adjacent
    segments.

    The class provides attributes to access the name of the path, the sequence
    of segment names along with their orientations, and the overlaps between
    consecutive segments. It extends the functionality of its parent class
    BaseHelper with specific processing for Path records.
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
    """
    Provides functionality for working with GFA (Graphical Fragment Assembly) data.

    The GFA class represents a graphical fragment assembly structure often used in
    bioinformatics for representing sequences and their relationships in a graph format.
    It facilitates the conversion of GFA entities into other formats, such as NetworkX
    MultiGraphs and FASTA, as well as writing GFA data to files in a standardized way.
    The class supports various operations for handling nodes, edges, paths, and metadata
    within the GFA.

    Attributes:
        comments: List[str]
            A list of comments or metadata that are included at the beginning of the GFA data.
        header: Dict
            Dictionary containing GFA header information, such as version number ('VN').
        segments: Dict
            Represents segment entries, where each key is a segment name and its value is an
            object containing sequence information and optional attributes.
        links: Dict
            Stores link objects describing connections between segments. Each link specifies
            source/destination orientations, overlap information, and can include additional
            attributes.
        containments: Dict
            Defines containment relationships between segments. Each containment specifies
            involved segments, orientations, positions, and optional attributes.
        paths: Dict
            Holds paths through the assembly graph. Each path associates segment names, overlaps,
            and other optional attributes.
    """

    def to_networkx(self,
                    include_seq: bool=False,
                    annotate_paths: bool=False,
                    collections_to_str: bool=True,
                    add_label: bool=False,
                    progress: bool=False) -> nx.MultiGraph:
        """
        Converts the instance into a NetworkX MultiGraph.

        This method generates a MultiGraph representation using the NetworkX library.
        The graph can include node sequences, labeled nodes, annotated paths, and collection
        attributes converted into strings based on the provided options. This representation
        serves as a useful format for graph-based computations or visualizations.

        Parameters:
            include_seq: bool
                If True, node sequences will be included as attributes in the graph.
            annotate_paths: bool
                If True, paths through the graph will be annotated on the nodes.
            collections_to_str: bool
                If True, collection attributes will be converted into strings to ensure
                compatibility with graph serialization formats like GraphML.
            add_label: bool
                If True, labels will be added to the nodes in the graph.
            progress: bool
                If True, progress bars will be displayed during node and edge creation.

        Returns:
            nx.MultiGraph
                A NetworkX MultiGraph representing the instance.
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

    def to_fasta(self, output: str | typing.TextIO) -> None:
        """
        Writes sequence data to a file or file-like object in FASTA format.

        Writes the segments contained in the object to the specified output in the FASTA format,
        including optional descriptive fields for each sequence. This method ensures proper
        handling of file objects and logs the number of sequences written for auditing purposes.

        Arguments:
            output: str | typing.TextIO
                The path to the file or a file-like object where the FASTA-formatted data will be written.
                If a string is provided, it is treated as a file path that will be opened and written to.

        Returns:
            None
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

    def write(self, output: str | typing.TextIO) -> None:
        """
        Writes the GFA (Graphical Fragment Assembly) structure data to a specified output. The output
        can be either a file path or a TextIO stream. This function serializes the different components
        of the GFA structure, such as comments, headers, segments, links, containments, and paths,
        into the GFA format before writing them to the output.

        Arguments:
            output: str | typing.TextIO
                The target output where the GFA data will be written. If a string is provided, it is
                treated as a file path where the data will be saved. If a TextIO stream is provided,
                the data will be written directly to the stream.

        Raises:
            None
        """
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

    def __init__(self,
                 filename: str,
                 ignore_isolate_paths: bool=False,
                 skip_sequence_data: bool=False,
                 progress: bool=False) -> None:
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