#!/usr/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1]

data = pd.read_csv(filename)

columns = data.columns 

columns = columns[[i.endswith("NamOfGen") for i in columns]]

ranges = []
for i in columns:
	ranges.extend([int(j) for j in i.split("_") if j.isnumeric()])

print(ranges)

matrix = []

for index_i, mutation_i in data.iterrows():
	gene_list_i = dict()
	for column in columns:
		try:
			gene_list_i[column] = set(data.loc[index_i, column].split("|"))
		except:
			gene_list_i[column] = None
	for index_j, mutation_j in data[index_i + 1:].iterrows():
		jaccard_indexis = []
		gene_list_j = dict()
		for column in columns:
			try:
				gene_list_j[column] = set(data.loc[index_j, column].split("|"))
			except:
				gene_list_j[column] = None
		for i in range(len(columns)):
			if gene_list_i[columns[i]] == None or gene_list_j[columns[i]] == None:
				jaccard_indexis.append(0)
				continue
			intersection = set.intersection(\
				gene_list_i[columns[i]], gene_list_j[columns[i]])
			union = set.union(\
				gene_list_i[columns[i]], gene_list_j[columns[i]])
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

