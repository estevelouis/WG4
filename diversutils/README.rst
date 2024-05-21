DiversUtils
===========

DiversUtils is a C code base to measure diversity in linguistic data.
Measuring diversity of linguistic resources as well as system
predictions allows for insights in the representativeness of data.

How to use
==========

Both C and Python interfaces are available. The C interface has better
performance and is recommended for large-scale diversity benchmarking,
as is the purpose of WG4’s efforts. The Python interface is better for
local-scale diversity measurement, and is intended to be usable for a
wider community, as well as to benefit of Python’s existing stack.

C interface
-----------

Repository structure
~~~~~~~~~~~~~~~~~~~~

``src`` contains the source code, with headers under ``src/include``.
Most important files include: \* ``src/include/graph.h`` and
``src/include/dfunctions.h`` for diversity functions \*
``src/include/cupt/*`` and ``src/include/jsonl/*`` for parsers \*
``src/include/udpipe_interface.hpp`` and
``src/include/udpipe_interface/*`` for UDPipe interface

Large-scale benchmarking
~~~~~~~~~~~~~~~~~~~~~~~~

In the standard use case, you have some data in files (currently,
``*.conllu``, ``*.cupt``, and ``*.jsonl`` are supported) and you wish to
assess how diverse they are, using a set of diversity functions.

makefile
^^^^^^^^

You may change default parameters in ``makefile``. Most parameters can
be overriden in CLI. For example, if in the makefile you wrote

.. code:: text

   ...
   ENABLE_DISPARITY_FUNCTIONS = 0
   ...

the default behaviour will be to ignore disparity functions, but when
running the program you may add ``--enable_disparity_functions=1`` in
the command line to change the behaviour.

You may compile the program with

.. code:: bash

   make measurement_full

which should yield a binary ``bin/main_measurement``. Once you’ve
defined all the files you with to process in ``measurement_files.txt``,
you may run the executable.

Outputs
^^^^^^^

It should generate a ``measurement_output.tsv``,
``measurement_output_timing.tsv``, and ``measurement_output_memory.tsv``
(or what is set respectively by ``OUTPUT_PATH``, ``OUTPUT_PATH_TIMING``,
and ``OUTPUT_PATH_MEMORY``).

Available macros
~~~~~~~~~~~~~~~~

Performance macros
^^^^^^^^^^^^^^^^^^

-  ``NUM_ROW_THREADS``: integer, number of threads to use if distance
   computation is done by row batches. See
   ``ENABLE_MULTITHREADED_ROW_GENERATION``.
-  ``NUM_MATRIX_THREADS``: integer, number of threads to use if distance
   computation is done directly for the full matrix. See
   ``ENABLE_MULTITHREADED_MATRIX_GENERATION``.
-  ``ROW_GENERATION_BATCH_SIZE``: integer, number of rows to compute at
   once if distance computation is done by row batches.
   ``ROW_GENERATION_BATCH_SIZE <= NUM_ROW_THREADS`` must be true in
   order to have at least one thread per computed row. See
   ``ENABLE_MULTITHREADED_ROW_GENERATION`` and ``NUM_ROW_THREADS``.
-  ``ENABLE_AVX256``: boolean, whether to activate AVX256 instructions.
   This vectorises operations and improves speed. **Setting
   ``ENABLE_AVX256=1`` may break builds on older machines.**

Diversity function macros
^^^^^^^^^^^^^^^^^^^^^^^^^

-  ``ENABLE_DISPARITY_FUNCTIONS``: boolean, whether to activate
   disparity functions. Listed below are the functions that can be
   deactivated.
-  ``ENABLE_STIRLING``: boolean.
-  ``ENABLE_RICOTTA_SZEIDL``: boolean.
-  ``ENABLE_PAIRWISE``: boolean.
-  ``ENABLE_LEXICOGRAPHIC``: boolean.
-  ``ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY``: boolean.
-  ``ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY``:
   boolean.
-  ``ENABLE_LEINSTER_COBBOLD_DIVERSITY``: boolean.
-  ``ENABLE_FUNCTIONAL_EVENNESS``: boolean.
-  ``ENABLE_FUNCTIONAL_DISPERSION``: boolean.
-  ``ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED``: boolean.
-  ``ENABLE_NON_DISPARITY_FUNCTIONS``: boolean, whether to activate
   non-disparity functions. Listed below are the functions that can be
   deactivated.

…

Linking
~~~~~~~

For linking purposes, the ``LD_LIBRARY_PATH`` environment variable must
contain both ``$HOME/.local/lib/diversutils`` and the standard lib path
of the system (usually ``/usr/lib/x86_64-linux-gnu``).

Python interface
----------------

A Python interface is being designed. By running

.. code:: bash

   python3 -m pip install -e .

``setup.py`` will be used to install the ``diversutils`` package.

Once in Python, you may run

.. code:: python

   import diversutils
   graph_index = diversutils.measurement_from_cfg("default.cfg")

which will read ``default.cfg`` and run all enabled diversity functions
on your data.

If you with to create your own data points in a script and measure their
diversity, you may run

.. code:: python

   import diversutils

   # Create an empty graph; graph_index serves as an identifier
   graph_index = diversutils.create_empty_graph(0, 0)

   # Add nodes, respectively with 1 and 5 in absolute proportions
   diversutils.add_node(graph_index, 1)
   diversutils.add_node(graph_index, 5)

   # Compute relative proportions of nodes
   diversutils.compute_relative_proportion(graph_index)

   # Compute diversities that you need
   entropy, hill_number = diversutils.individual_measure(graph_index, diversutils.DF_ENTROPY_SHANNON_WEAVER)
   entropy, hill_number = diversutils.individual_measure(graph_index, diversutils.DF_ENTROPY_RENYI, 2.0)

   # Free the underlying graph
   diversutils.free_graph(graph_index)

You may also compute disparities, but this requires setting some vector
space

.. code:: python

   import diversutils

   # Create an empty graph, specifying the number of dimensions
   graph_index = diversutils.create_empty_graph(0, 100)

   # Create a Word2Vec vector set
   w2v_index = diversutils.load_w2v("some/path/to/a/word2vec/binary.bin")

   # Bind the Word2Vec to the graph
   diversutils.bind_w2v(graph_index, w2v_index)

   # Add nodes, along with their key
   diversutils.add_node(graph_index, 1, "great")
   diversutils.add_node(graph_index, 5, "wonderful")

   # Compute relative proportions of nodes
   diversutils.compute_relative_proportion(graph_index)

   # Compute disparity
   pairwise = diversutils.individual_measure(graph_index, diversutils.DF_PAIRWISE)

   # Free the underlying graph
   diversutils.free_graph(graph_index)

   # Free the Word2Vec set
   diversutils.free_w2v(word2vec_index)

NOTE THAT THE PYTHON API IS UNSTABLE.

Available diversity functions
=============================

[…]

Licensing
=========

::

   Copyright (c) 2024  LISN / Université Paris-Saclay / CNRS  Louis Estève (louis.esteve@universite-paris-saclay.fr)
   All rights reserved.

This work was made possible by the financial support of the “Plan Blanc”
from Université Paris-Saclay.
