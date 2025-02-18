#!/usr/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

files = sys.argv[1:]

range_ = 3
column = f"protein_coding_range_{range_}_NamOfGen"

protein_coding = pd.read_csv("additional/proteinCodingGenes.csv", usecols=["gene_name"])["gene_name"]

result_file = pd.DataFrame(index=files, columns=protein_coding, data = 0)

for file in files:
	gene_set = set()
	genes = pd.read_csv(file, usecols=[column])
	genes.dropna(inplace=True)
	for index, gene in genes.iterrows():
		gene_set = gene_set | set(*[j.split("|") for j in gene])
	gene_list = list(gene_set)
	result_file.loc[file, gene_list] = 1

result_file.to_csv("gene_similarity_matrix.csv")
