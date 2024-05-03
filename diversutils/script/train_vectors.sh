#!/bin/bash
echo "processing lang $1"
mkdir -p data
./bin/jsonl_to_word2vec_format data/$1_raw.txt data/$1*.jsonl
./word2vec/bin/word2vec -train data/$1_raw.txt -output word2vec/bin/$1_vec.bin -cbow 1 -size 100 -window 10 -negative 10 -hs 0 -threads 6 -binary 1 -iter 3 -min-count 1
