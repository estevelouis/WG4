# DiversUtils

```
Copyright (c) 2024  LISN / Université Paris-Saclay / CNRS  Louis Estève (louis.esteve@universite-paris-saclay.fr)
All rights reserved.
```

This work was made possible by the financial support of the "Plan Blanc" from
Université Paris-Saclay.

# How to use

In the `makefile` you may change macro-parameters; this includes which functions
to test, multithreading, as well as other parameters.

For the use of disparity functions, `W2V_PATH` should be defined, pointing to
a Word2Vec-formatted binary.

You may then run `make measurement_standard`, which will compile to
`bin/main_measurement`.
Go to `bin` and create a file `measurement_files.txt` (or whatever file name
you put as `INPUT_PATH`).
This file must contain the names of all files that are to be processed (one per
line); they can be of `*.cupt` (extension of `*.conllu`) or `*.jsonl` format.
Lines starting with `#` are ignored.

The program will make measurements at multiple occasions of parsing the data,
based on how you defined recomputing options in the `makefile`.
This yields a `measurement_output.tsv` (or whatever file name you put in
`OUTPUT_PATH`).
