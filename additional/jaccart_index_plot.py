#!/usr/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

filename =  sys.argv[1]
data = pd.read_csv(filename)
columns = data.columns 
columns = columns[[i.endswith("NamOfGen") for i in columns]]
ranges = []
for i in columns:
	ranges.extend([int(j) for j in i.split("_") if j.isnumeric()])


gene_lists = []
for i, val_i in data.iterrows():
	gene_list = []
	for j in columns:
		try:
			gene_list.append(set(val_i[j].split("|")))
		except:
			gene_list.append(None)
	gene_lists.append(gene_list)

genes = pd.DataFrame(columns=ranges, data = gene_lists)
print(ranges)

matrix = []
for index_i, mutation_i in genes.iterrows():
	for index_j, mutation_j in genes[index_i + 1:].iterrows():
		jaccard_indexis = []
		for i in range(len(ranges)):
			if genes.loc[index_i, ranges[i]] == None or\
					genes.loc[index_j, ranges[i]] == None:
				jaccard_indexis.append(0)
				continue
			intersection = set.intersection(\
				genes.loc[index_i, ranges[i]], genes.loc[index_j, ranges[i]])
			union = set.union(\
				genes.loc[index_i, ranges[i]], genes.loc[index_j, ranges[i]])
			jaccard_index = round(len(intersection)/len(union), 3)
			jaccard_indexis.append(jaccard_index)
			#print(f"{index_i}-{index_j}: {jaccard_index}")
		#print(jaccard_indexis)
		if sum(jaccard_indexis) == 0:
			continue
		matrix.append(jaccard_indexis) #plt.errorbar(ranges, jaccard_indexis)


matrix = pd.DataFrame(matrix)

plt.boxplot(matrix, tick_labels=ranges)
plt.xlabel("Ranges")
plt.ylabel("Jaccard index")
plt.title("Jaccard Index of Mutation Pairs based on Genes for Different Ranges ")
#plt.xticks(ranges)
plt.show()

