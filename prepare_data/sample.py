
from pathlib import Path
import time
import io
import zstandard
import jsonlines
import argparse
import csv
import pandas as pd
import random

NUM_TEXT_1M = 1000000
NUM_TEXT_5M = 5000000
NUM_TEXT_10M = 10000000


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process language data.")
    parser.add_argument('--lang_code', type=str, default='mt',
                        help="Specify the language code (default: mt)")
    parser.add_argument('--sample_rate', type=float, default=0.1)
    parser.add_argument('--output_dir', type=str,
                        default='/scratch/project_462000447/members/shaoxion/WG4/hpltv1.2_samples')
    return parser.parse_args()



def proc_monohplt(data_dir, output_dir, lang_code, parsed_args):
    print(f"Processing {lang_code}...")

    zst_files = [file for file in data_dir.iterdir() if file.suffix == '.zst']
    num_texts = 0
    
    lang_output_dir = Path(output_dir)/lang_code
    if not lang_output_dir.exists():
        lang_output_dir.mkdir()
    print("output dir:", str(lang_output_dir))
    
    for zst_file_path in zst_files:
        if (lang_output_dir/zst_file_path.stem).exists(): 
            continue
        print(zst_file_path)
        texts = []

        with open(zst_file_path, 'rb') as compressed_file:
            decompressor = zstandard.ZstdDecompressor()
            with decompressor.stream_reader(compressed_file) as reader:
                buffered_reader = io.BufferedReader(reader)
                with jsonlines.Reader(buffered_reader) as json_reader:
                    texts = list(json_reader)

        
        sample_size = int(len(texts) * parsed_args.sample_rate)
        texts = random.sample(texts, sample_size)
        with jsonlines.open(lang_output_dir/f'{zst_file_path.stem}', mode='w') as writer:
            writer.write_all(texts)
        num_texts += len(texts)

    return {"lang_code": lang_code,
            "num_texts": num_texts}


    
if __name__ == "__main__":

    args = parse_arguments()

    data_dir = Path("/scratch/project_462000506/source_data/hplt-v1.2/cleaned")
    # folder_names = [item.name for item in data_dir.iterdir() if item.is_dir()]

    csv_file_path = "./stats.csv"
    if Path(csv_file_path).exists():
        languages_done = pd.read_csv(csv_file_path)['lang_code'].to_list()
        print("languages done:", languages_done)
    else:
        with open(csv_file_path, mode='a', newline='', encoding='utf-8') as csv_file:
            writer = csv.DictWriter(
                csv_file, fieldnames=['lang_code', 'num_texts'])
            writer.writeheader()
        languages_done = []

    lang_code = args.lang_code
    print("*"*89)
    print("Processing", lang_code)
    if lang_code in languages_done:
        exit()
    start_time = time.time()
    print(str(data_dir/lang_code))
    stats_dict = proc_monohplt(data_dir=data_dir/lang_code, output_dir=args.output_dir, lang_code=lang_code, parsed_args=args)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Elapsed time: {} seconds".format(elapsed_time))
    if stats_dict != 0 and stats_dict is not None:
        with open(csv_file_path, mode='a', newline='', encoding='utf-8') as csv_file:
            writer = csv.DictWriter(
                csv_file, fieldnames=list(stats_dict.keys()))
            if not Path(csv_file_path).exists():
                writer.writeheader()
            writer.writerow(stats_dict)
