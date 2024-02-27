#
#      DiversUtils - Functions to measure diversity
#
# Copyright (c) 2024  LISN / Université Paris-Saclay / CNRS  Louis Estève (louis.esteve@universite-paris-saclay.fr)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import lingua
import conllu
import numpy as np
import math

PREFIX_PATH = "/home/esteve/Documents/thesis/workrepo/ups_lingdiv_utils/download_data/universal_dependencies/ud-treebanks-v2.13/"

CODE_SWITCHING_TREEBANKS = {
	"UD_Turkish_German-SAGT": {
		"files": {
			"train": "qtd_sagt-ud-train.conllu",
			"dev": "qtd_sagt-ud-dev.conllu",
			"test": "qtd_sagt-ud-test.conllu",
		},
		"languages": [
			lingua.Language.TURKISH,
			lingua.Language.GERMAN,
		]
	},
	"UD_Maghrebi_Arabic_French-Arabizi": {
		"files": {
			"train": "qaf_arabizi-ud-train.conllu",
			"dev": "qaf_arabizi-ud-dev.conllu",
			"test": "qaf_arabizi-ud-test.conllu",
		},
		"languages": [
			lingua.Language.ARABIC,
			lingua.Language.FRENCH,
		]
	}
}

#TARGET_TREEBANK = "UD_Turkish_German-SAGT"
TARGET_TREEBANK = "UD_Maghrebi_Arabic_French-Arabizi"

def code_switching_m_index(dict_abundances):
	k = len(dict_abundances)
	num_tokens = sum(dict_abundances[i] for i in dict_abundances)
	sigma_sum = sum(map(lambda x:(x / num_tokens) ** 2, (dict_abundances[i] for i in dict_abundances)))
	result = (1 - sigma_sum) / ((k - 1) * sigma_sum)
	return result

def code_switching_language_entropy(dict_abundances):
	logarithmic_base = 2.0
	k = len(dict_abundances)
	num_tokens = sum(dict_abundances[i] for i in dict_abundances)
	sum_ = 0.0
	for i in dict_abundances:
		if dict_abundances[i] == 0:
			continue
		p_i = dict_abundances[i] / num_tokens
		sum_ += p_i * (math.log(p_i) / math.log(logarithmic_base))
	result = -sum_
	return result
	

def main():
	if TARGET_TREEBANK not in CODE_SWITCHING_TREEBANKS:
		print("Please select a treebank defined in CODE_SWITCHING_TREEBANKS")
		return 1

	detector = lingua.LanguageDetectorBuilder.from_languages(*(CODE_SWITCHING_TREEBANKS[TARGET_TREEBANK]["languages"])).build()


	treebank_language_counts = {}

	treebank_i_index_upper_count = 0
	treebank_i_index_lower_count = 0
	treebank_i_index_previous_value = None

	treebank_cmi_values = []

	treebank_path = f"{PREFIX_PATH}{TARGET_TREEBANK}/"
	for subtreebank in CODE_SWITCHING_TREEBANKS[TARGET_TREEBANK]["files"]:
		current_path = f"{treebank_path}{CODE_SWITCHING_TREEBANKS[TARGET_TREEBANK]['files'][subtreebank]}"
		print("=" * 32)
		print(f"{subtreebank}:", current_path)
		f = open(current_path, "rt", encoding="utf-8")
		r = f.read()
		f.close()

		subtreebank_language_counts = {}

		subtreebank_i_index_upper_count = 0
		subtreebank_i_index_lower_count = 0
		subtreebank_i_index_previous_value = treebank_i_index_previous_value

		subtreebank_cmi_values = []

		sentences = conllu.parse(r)
		for sent in sentences:
			#print(dir(sent))
			#breakpoint()
			sent_raw = "".join([f"{tok['form']} " if tok["misc"] == None or "SpaceAfter" not in tok["misc"] or tok["misc"]["SpaceAfter"].lower() != "no" else tok["form"] for tok in sent])
			sent_lang = detector.detect_language_of(sent_raw)

			sent_language_counts = {}

			#print(f"sent_lang: {sent_lang}")
			for tok in sent:
				tok_lang = detector.detect_language_of(tok["form"])
				#print(f"tok_lang: {tok_lang}")
				if tok_lang != None:
					if tok_lang.iso_code_639_3.name not in subtreebank_language_counts:
						subtreebank_language_counts[tok_lang.iso_code_639_3.name] = 1
					else:
						subtreebank_language_counts[tok_lang.iso_code_639_3.name] += 1

					# CMI part
					if tok_lang.iso_code_639_3.name not in sent_language_counts:
						sent_language_counts[tok_lang.iso_code_639_3.name] = 1
					else:
						sent_language_counts[tok_lang.iso_code_639_3.name] += 1

					# I-index part
					if subtreebank_i_index_previous_value != None:
						if tok_lang.iso_code_639_3.name != subtreebank_i_index_previous_value:
							subtreebank_i_index_upper_count += 1
						subtreebank_i_index_lower_count += 1
					subtreebank_i_index_previous_value = tok_lang.iso_code_639_3.name
					if treebank_i_index_previous_value != None:
						if tok_lang.iso_code_639_3.name != treebank_i_index_previous_value:
							treebank_i_index_upper_count += 1
						treebank_i_index_lower_count += 1
					treebank_i_index_previous_value = tok_lang.iso_code_639_3.name

			local_cmi = (sum(sent_language_counts[i] for i in sent_language_counts) - max(sent_language_counts[i] for i in sent_language_counts)) / sum(sent_language_counts[i] for i in sent_language_counts) # is this really ok simply on the basis that tokens without language are ignored?
			subtreebank_cmi_values.append(local_cmi)

			#break
		for key in subtreebank_language_counts:
			print(f"{key}: {subtreebank_language_counts[key]}")
			if key not in treebank_language_counts:
				treebank_language_counts[key] = subtreebank_language_counts[key]
			else:
				treebank_language_counts[key] += subtreebank_language_counts[key]

		treebank_cmi_values += subtreebank_cmi_values

		local_m_index = code_switching_m_index(subtreebank_language_counts)
		print(f"m_index: {local_m_index}")
		local_i_index = subtreebank_i_index_upper_count / subtreebank_i_index_lower_count
		print(f"i_index: {local_i_index}")
		local_cmi_avg = np.mean(subtreebank_cmi_values)
		local_cmi_std = np.std(subtreebank_cmi_values)
		print(f"cmi, avg: {local_cmi_avg}; std: {local_cmi_std}")
		local_language_entropy = code_switching_language_entropy(subtreebank_language_counts)
		print(f"language_entropy: {local_language_entropy}")
		#break
	print("=" * 32)
	print("global:")
	for key in treebank_language_counts:
		print(f"{key}: {treebank_language_counts[key]}")
	local_m_index = code_switching_m_index(treebank_language_counts)
	print(f"m_index: {local_m_index}")
	local_i_index = treebank_i_index_upper_count / treebank_i_index_lower_count
	print(f"i_index: {local_i_index}")
	local_cmi_avg = np.mean(treebank_cmi_values)
	local_cmi_std = np.std(treebank_cmi_values)
	print(f"cmi, avg: {local_cmi_avg}; std: {local_cmi_std}")
	local_language_entropy = code_switching_language_entropy(treebank_language_counts)
	print(f"language_entropy: {local_language_entropy}")
	return 0

if __name__ == "__main__":
	exit(main())
