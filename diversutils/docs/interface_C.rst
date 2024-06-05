=============
C Interface
=============

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

TODO: group the macros by types of functions?

Use cases
^^^^^^^^^^

For disparity functions, you must first enable disparities generally (ENABLE_DISPARITY_FUNCTIONS or --enable_disparity_functions) and then enable the individual disparities you want (for example with ENABLE_PAIRWISE or --enable_pairwise).

Linking
~~~~~~~

For linking purposes, the ``LD_LIBRARY_PATH`` environment variable must
contain both ``$HOME/.local/lib/diversutils`` and the standard lib path
of the system (usually ``/usr/lib/x86_64-linux-gnu``).
