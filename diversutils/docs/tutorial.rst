
Benchmarking Diversity Measures
================================================================


1 Extract texts from jsonl into txt
------------------------------------

To extract texts from JSONL into a TXT file, use the `jsonl_to_word2vec_format` process. This process operates on the .jsonl files located under the `data/` directory, starting with the two letters of the language. It generates a .txt file containing `TARGET_NUMBER_DOCUMENTS` (defined in the makefile, default 1M) documents, where each "document" corresponds to a line in the JSONL file. Subsequently, it trains a Word2Vec model on this data and stores the resulting model under `word2vec/bin/`, with the file name starting with the two letters of the language.

2 Train word2vec for each language
----------------------------------

This step requires the use of the `word2vec` tool available at https://github.com/dav/word2vec. To set up the `word2vec` tool, execute either of the following commands:

- `make setup_word2vec_learner`
- `git clone https://github.com/dav/word2vec.git` followed by `cd word2vec` and `make build`

In this case, the training algorithm adopted is cbow (Continuous Bag of Words), and each word vector will have 100 dimensions.

3 Diversity
-----------

For device-specific parameters, consider the following configurations:

- `ENABLE_AVX256=1`: This setting vectorizes computation on machines capable of AVX256.
- `NUM_ROW_THREADS`, `NUM_MATRIX_THREADS`, and `ROW_GENERATION_BATCH_SIZE`: Adjust these parameters to match the number of threads available on the device for optimal performance.

Regarding diversity functions:

- For the non-disparity run, set `ENABLE_DISPARITY_FUNCTIONS` to 0 and `ENABLE_NON_DISPARITY_FUNCTIONS` to 1. For the disparity run, do the opposite.

Additionally, the process supports taking snapshots at logarithmic intervals instead of linear for languages with massive amounts of data.

Data parameters:

- `DOCUMENT_COUNT_RECOMPUTE_STEP = 1`
- `DOCUMENT_RECOMPUTE_STEP_USE_LOG10 = 1`
- `DOCUMENT_COUNT_RECOMPUTE_STEP_LOG10 = 0.1`

Note: The provided URLs are for reference and may need to be updated based on the actual location of the resources.

