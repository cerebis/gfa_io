# GFA IO module

Currently this module only reads GFA v1 and does not write either v1 or v2.

## Motivations:
- Python 2.7 support.
- Provide a simple means of loading GFA files into a standard framework (Networkx).
- Avoid problematic modules

## Limitations
- cannot write any GFA format.
- appears to be resource hungry for large assemblies.

## Example

There is a small UX included in the module for testing. This will read a GFAv1 file and write GraphML via Networkx. 
Networkx appears to be a tad hungry when writing this format and there are constraints on what types it is 
capable of serialising.

The following will read a SPADES GFAv1 file and write GraphML, where `--paths` are included as node attributes. Since
SPADES includes empty (non-compliant) paths in its GFA file, the `--ignore` flag is used to skip over them.
```bash
> python gfa_io/gfa_io.py --paths --ignore $GFA_IN $GRAPHML_OUT
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
