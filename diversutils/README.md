# DiversUtils

DiversUtils is a C code base to measure diversity in linguistic data.
Measuring diversity of linguistic resources as well as system predictions
allows for insights in the representativeness of data.

# How to use

## Old interface

(this is the old C interface, which is technically usable but not very
friendly; see the new Python interface below for more convenience)

In the standard use case, you have some data in files (currently, `*.conllu`,
`*.cupt`, and `*.jsonl` are supported) and you wish to assess how diverse they
are.
You may first compile with
```bash
make measurement_standard
```
which should yield a binary under `bin/`.
Once you've defined all the files you with to process in
`measurement_files.txt`, you may run the executable.
It should generate a `measurement_output.tsv`.

## New interface

A Python interface is being built.
You can access it by running
```bash
make measurement_python_interface
```
(this may require changing some paths in the `makefile`, especially with
regards to Python headers)

Once in Python, you may run
```python
import diversutils
graph_index = diversutils.measurement_from_cfg("default.cfg")
```
which will read `default.cfg` and run all enabled diversity functions on your
data.

If you with to create your own data points in a script and measure their
diversity, you may run
```python
import diversutils

# Create an empty graph; graph_index serves as an identifier
graph_index = diversutils.create_empty_graph(0, 0)

# Add nodes, respectively with 1 and 5 in absolute proportions
diversutils.add_node(graph_index, 1)
diversutils.add_node(graph_index, 5)

# Compute relative proportions of nodes
diversutils.compute_relative_proportion(graph_index)

# Compute diversities that you need
entropy, hill_number = diversutils.individual_measure(graph_index, diversutils.ENTROPY_SHANNON_WEAVER)
entropy, hill_number = diversutils.individual_measure(graph_index, diversutils.ENTROPY_RENYI, 2.0)

# Free the underlying graph
diversutils.free_graph(graph_index)
```

NOTE THAT THE PYTHON API IS UNSTABLE.

# Available diversity functions

[...]

# Licensing

```
Copyright (c) 2024  LISN / Université Paris-Saclay / CNRS  Louis Estève (louis.esteve@universite-paris-saclay.fr)
All rights reserved.
```

This work was made possible by the financial support of the "Plan Blanc" from
Université Paris-Saclay.
