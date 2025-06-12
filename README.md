# GFA IO module

Currently this module reads and writes GFA v1

## Motivations:
- Python 3.7 support.
- Provide a simple means of loading GFA files into a standard framework (Networkx).
- Avoid problematic modules

## Example

A small utility `gfa_utils` is included in the module, which included the following subcommands.

### update-segments
Replace sequences in the GFA file with those in a supplied FASTA file. It is expected that the sequence
records match between the GFA and FASTA file.

One use-case is updating a GFA with post-assembly polished sequences.

```bash
gfa_utils update-segments GFA_IN FASTA_IN GFA_OUT
```

### isolates

Analyse a GFA file for segments which are network isolates (unattached) and report in CSV format.

Additional simple criteria can be imposed on the result, such as that a segment must refer to itself (circular), min and maximul segment length.

```bash
gfa_utils isolates GFA_IN ISOLATE_CSV
```

### dump-segments

Extract the segment sequences in FASTA format.

```bash
gfa_utils dump-segments GFA_IN FASTA_OUT
``` 

### convert

The following will read a GFA file and write a corresponding graph using NetworkX. File format can be one of GraphML, GML or an edge list.

```bash
> python gfa_utils convert -f graphml $GFA_IN $GRAPHML_OUT
```



In a piece of code, one can get a networkx.DiGraph in the following manner

```python
import gfa_io
import networkx

# read in the GFA file
gfa = gfa_io.GFA('assem.gfa', ignore_isolate_paths=True)
# create a DiGraph()
g = gfa.to_graph(annotate_paths=True, collections_to_str=True, include_seq=True)
# write to GraphML format
networkx.write_graphml(g, 'assem.graphml')
```
