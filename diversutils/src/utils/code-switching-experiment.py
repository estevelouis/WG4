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

ISO_639_3_TO_ISO_639_2_MAP = {
	"FRA": "FRE",
	"ARA": "ARA",
	"TUR": "TUR",
	"DEU": "GER",
}

class PhylNode():
	def __init__(self, key="", forth_degree=0, neighbours=[], parent=None):
		self.key = key
		self.forth_degree = forth_degree
		self.neighbours = neighbours[:] # copy matters for reasons of pointers
		self.parent = parent

path_counter = 0
#num_children_per_node = {}

def compute_path_recursively_phyl_node(pn, all_paths, current_depth, current_path, reset_path_counter):
	global path_counter
	#global num_children_per_node

	if reset_path_counter:
		path_counter = 0

	current_path.append(pn.key)

	all_paths.append({"key": pn.key, "keys": current_path[:]}) # copy matters for reasons of pointers, as after we pop
	print(all_paths[-1])

	path_counter += 1
	for i in range(len(pn.neighbours)):
		compute_path_recursively_phyl_node(pn.neighbours[i], all_paths, current_depth + 1, current_path, 0)

	current_path.pop()


def parse_nex_file(path):
	num_created_nodes = 0

	pn_root = PhylNode(key="root")

	num_created_nodes += 1

	f = open(path, "rt", encoding="utf-8")
	
	readlines = f.readlines()

	reached_begin_tree = 0

	for line in readlines:
		if "BEGIN TREES;" in line:
			reached_begin_tree = 1
			break

	if not reached_begin_tree:
		f.close()
		raise Exception("no 'BEGIN TREES;'")

	depth = 0
	max_depth = 0

	for line in readlines:
		if "END;" in line:
			break

		index_tree = line.find("tree ")
		if index_tree == -1:
			continue
		index_tree_name = index_tree + 5
		index_equal = line.find("=", index_tree_name)
		index_tree_name_end = index_equal

		pn_current = pn_root

		index_char = index_equal

		while not (index_char == len(line) or line[index_char] == "("):
			index_char += 1

		while True:
			while line[index_char] != ";":
				index_end = -1
				index_next_parenthesis_closing = -1
				index_next_parenthesis_opening = -1
				index_next_comma = -1
				index_next_semi_colon = -1
				if line[index_char] == '(':
					local_node = PhylNode()
					local_node.parent = pn_current
					pn_current.neighbours.append(local_node)
					pn_current.forth_degree += 1
					pn_current = local_node
					depth += 1
					if depth > max_depth:
						max_depth = depth
					index_char += 1
				elif line[index_char] == ')':
					index_char += 1

					index_end = -1
					index_next_parenthesis_closing = line.find(")", index_char)
					index_next_parenthesis_opening = line.find("(", index_char)
					index_next_comma = line.find(",", index_char)
					index_next_semi_colon = line.find(";", index_char)

					if index_next_parenthesis_closing > -1 and (index_end == -1 or index_next_parenthesis_closing < index_end):
						index_end = index_next_parenthesis_closing
					if index_next_parenthesis_opening > -1 and (index_end == -1 or index_next_parenthesis_opening < index_end):
						index_end = index_next_parenthesis_opening
					if index_next_comma > -1 and (index_end == -1 or index_next_comma < index_end):
						index_end = index_next_comma
					if index_next_semi_colon > -1 and (index_end == -1 or index_next_semi_colon < index_end):
						index_end = index_next_semi_colon

					
					key_len = 0
					if index_end == -1:
						key_len = len(line) - index_char
					else:
						key_len = index_end - index_char
					

					pn_current.key = line[index_char:index_char+key_len]

					if index_end == -1:
						index_char += key_len
					else:
						index_char = index_end

					pn_current = pn_current.parent

					depth -= 1
				elif line[index_char] == ',':
					index_char += 1
				elif line[index_char] == ';':
					pass
				else:
					index_end = -1
					index_next_parenthesis_closing = line.find(")", index_char)
					index_next_parenthesis_opening = line.find("(", index_char)
					index_next_comma = line.find(",", index_char)
					index_next_semi_colon = line.find(";", index_char)

					if index_next_parenthesis_closing > -1 and (index_end == -1 or index_next_parenthesis_closing < index_end):
						index_end = index_next_parenthesis_closing
					if index_next_parenthesis_opening > -1 and (index_end == -1 or index_next_parenthesis_opening < index_end):
						index_end = index_next_parenthesis_opening
					if index_next_comma > -1 and (index_end == -1 or index_next_comma < index_end):
						index_end = index_next_comma
					if index_next_semi_colon > -1 and (index_end == -1 or index_next_semi_colon < index_end):
						index_end = index_next_semi_colon

					
					key_len = 0
					if index_end == -1:
						key_len = len(line) - index_char
					else:
						key_len = index_end - index_char

					local_node = PhylNode()

					num_created_nodes += 1

					local_node.key = line[index_char:index_char+key_len]

					pn_current.neighbours.append(local_node)
					pn_current.forth_degree += 1

					if index_end == -1:
						index_char += key_len
					else:
						index_char = index_end
			if index_char >= len(line):
				break
			if line[index_char] == ";":
				break

	f.close()

	current_path = []

	all_paths = []

	compute_path_recursively_phyl_node(pn_root, all_paths, 0, current_path, 1)

	dict_phyl_distances = {}
	for i in all_paths:
		dict_phyl_distances[i["key"]] = {}
		for j in all_paths:
			k = 0
			k_limit = min(len(i["keys"]), len(j["keys"]))
			while k < k_limit:
				if i["keys"][k] == j["keys"][k]:
					k += 1
				else:
					break
			dict_phyl_distances[i["key"]][j["key"]] = max(len(i["keys"]), len(j["keys"])) - k
			#print(i["key"], j["key"], dict_phyl_distances[i["key"]][j["key"]])


	num_children_per_node = {}
	for i in all_paths:
		for subpath_length in range(1, len(i["keys"]) + 1):
			subpath = i["keys"][:subpath_length]
			#if tuple(*subpath) not in num_children_per_node:
			if tuple(subpath) not in num_children_per_node:
				num_children_per_node[*subpath] = 1
			else:
				num_children_per_node[*subpath] += 1

	return ({i["key"]: i["keys"] for i in all_paths}, dict_phyl_distances, num_children_per_node)

def scheiner_phylogenetic_diversity(n_i, l_i, q=2.0):
	sum_standard = 0.0
	for i in range(len(n_i)):
		sum_standard += n_i[i] * l_i[i]
	sum_global = 0.0
	for i in range(len(n_i)):
		sum_global += ((n_i[i] * l_i[i]) / sum_standard) ** q
	return sum_global ** (1 / (1 - q))

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

	# --------

	key_to_phyl_path, dict_phyl_distances, num_children_per_node = parse_nex_file("/home/esteve/Documents/thesis/workrepo/publishing/survey_diversity/cldf-datasets-wals-878ea47/cldf/genealogy.nex")
	for key, value in num_children_per_node.items():
		print(key, value)

	# --------

	detector = lingua.LanguageDetectorBuilder.from_languages(*(CODE_SWITCHING_TREEBANKS[TARGET_TREEBANK]["languages"])).build()


	treebank_language_counts = {}

	treebank_i_index_upper_count = 0
	treebank_i_index_lower_count = 0
	treebank_i_index_previous_value = None

	treebank_cmi_values = []

	treebank_iso_639_2_counts = {}

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

		subtreebank_iso_639_2_counts = {}

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

					# Scheiner part
					"""
					if tok_lang.iso_code_639_2.name not in subtreebank_iso_639_2_counts:
						subtreebank_iso_639_2_counts[tok_lang.iso_code_639_2.name] = 1
					else:
						subtreebank_iso_639_2_counts[tok_lang.iso_code_639_2.name] += 1
					"""
					if tok_lang.iso_code_639_3.name in ISO_639_3_TO_ISO_639_2_MAP:
						iso_code_639_2 = ISO_639_3_TO_ISO_639_2_MAP[tok_lang.iso_code_639_3.name]
						if iso_code_639_2 not in subtreebank_iso_639_2_counts:
							subtreebank_iso_639_2_counts[iso_code_639_2] = 1
						else:
							subtreebank_iso_639_2_counts[iso_code_639_2] += 1

			local_cmi = (sum(sent_language_counts[i] for i in sent_language_counts) - max(sent_language_counts[i] for i in sent_language_counts)) / sum(sent_language_counts[i] for i in sent_language_counts) # is this really ok simply on the basis that tokens without language are ignored?
			subtreebank_cmi_values.append(local_cmi)

			#break
		for key in subtreebank_language_counts:
			print(f"{key}: {subtreebank_language_counts[key]}")
			if key not in treebank_language_counts:
				treebank_language_counts[key] = subtreebank_language_counts[key]
			else:
				treebank_language_counts[key] += subtreebank_language_counts[key]

		for key in subtreebank_iso_639_2_counts:
			if key not in treebank_iso_639_2_counts:
				treebank_iso_639_2_counts[key] = subtreebank_iso_639_2_counts[key]
			else:
				treebank_iso_639_2_counts[key] += subtreebank_iso_639_2_counts[key]

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

	n_i = []
	l_i = []
	for key in treebank_iso_639_2_counts:
		print(f"{key}: {treebank_iso_639_2_counts[key]}; {key_to_phyl_path[key.lower()]}")
		n_i.append(treebank_iso_639_2_counts[key])
		l_i.append(sum(1.0/num_children_per_node[*(key_to_phyl_path[key.lower()][:i])] for i in range(1, len(key_to_phyl_path[key.lower()]) + 1)))
	scheiner_res = scheiner_phylogenetic_diversity(n_i, l_i, q=2.0)
	print(f"scheiner_phylogenetic_diversity: {scheiner_res}")

	return 0

if __name__ == "__main__":
	exit(main())
