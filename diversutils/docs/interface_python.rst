=========================
Python Interface
=========================


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