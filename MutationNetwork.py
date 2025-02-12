#!/usr/bin/python3

from multiprocessing import Pool
import networkx as nx
import pandas as pd
import numpy as np
import argparse
import pickle
import time
import sys
import os


def read_bedpe_file(filename):
	# read bedpe file
	bedpe_file = pd.read_csv(filename, delimiter="\t",\
			names = ["col1_chr", "start1", "end1", "col2_chr", "start2", "end2", "p2"],\
			header=None, dtype=\
			{0: object, 1: int, 2: int, 3: object, 4: int, 5: int, 6: int})
	
	bedpe_file["col1_chr"] = bedpe_file["col1_chr"].str.replace("chr", "")
	bedpe_file["col2_chr"] = bedpe_file["col2_chr"].str.replace("chr", "")
	return bedpe_file


# check this part to optimize !!!
def filter_bedpe_file(bedpe_file, chromosome):
	# filter based on chromosome
	selected_intervals = bedpe_file[ bedpe_file["col1_chr"] == chromosome ].copy()
	selected_intervals.drop(columns=["col1_chr", "col2_chr"], inplace=True)
	selected_intervals.reset_index(drop=True, inplace=True)
	indexStop = selected_intervals.index.stop
	
	scores = selected_intervals["p2"]
	
	# add positive index to intervals on left
	selected_intervals.insert(2, "p1", range(1, 1 + indexStop ))
	
	# add negative index to intervals on right
	selected_intervals["p2"] = range(-1, -1 -1 * indexStop, -1 )
	
	# separate them into different arrays
	int1 = selected_intervals.iloc[:, :3].values
	int2 = selected_intervals.iloc[:, 3:].values
	del selected_intervals
	
	# merge two array into one array, 0th index will be (0,0,0) array
	# In each sub array, the rightmost value will be its index (including minus index)
	selected_intervals = np.concatenate(( int1, int2[::-1] ))
	return selected_intervals, scores


def sort_intervals(selected_intervals):
	sorted_intervals = selected_intervals[selected_intervals[:, 0].argsort()]
	
	return sorted_intervals


def intervals_to_array(intervals):
	lt = len( intervals )
	
	array = np.full( lt + 1, None, dtype=set)
	for i in range( lt +1):
		array[i] = set()
	
	i = -1
	for inv in intervals:
		i += 1
		j = i + 1
		while j < lt and inv[1] > intervals[j][0] :
			index1 = inv[2]
			index2 = intervals[j][2]
			
			array[ index1 ].add(index2)
			array[ index2 ].add(index1)
			
			j += 1
	
	return array


def write_file(bedpe_file_name):
	
	bedpe_file = read_bedpe_file(bedpe_file_name)
	base_file_name = os.path.basename(bedpe_file_name).split(".")[0]
	chromosomes = bedpe_file.loc[:, "col1_chr"].unique()
	
	# add chromosome list as file
	for chromosome in chromosomes:
		
		filtered_intervals, scores = filter_bedpe_file(bedpe_file, chromosome)
		sorted_intervals = sort_intervals(filtered_intervals)
		array = intervals_to_array(sorted_intervals)
		
		with open(f".pickles/{base_file_name}_{chromosome}_intervals.pickle", "wb") as f:
			if verbose:
				print(f"{base_file_name}_{chromosome}_intervals.pickle has been written")
			pickle.dump(sorted_intervals, f)
		with open(f".pickles/{base_file_name}_{chromosome}_array.pickle", "wb") as f:
			if verbose:
				print(f"{base_file_name}_{chromosome}_array.pickle has been written")
			pickle.dump(array, f)
		with open(f".pickles/{base_file_name}_{chromosome}_scores.pickle", "wb") as f:
			if verbose:
				print(f"{base_file_name}_{chromosome}_scores.pickle has been written")
			pickle.dump(scores, f)
	np.save(f".pickles/{base_file_name}_chromosomes", chromosomes)
	
	return 0


def check_pickle_file(bedpe_file_name):
	base_file_name = os.path.basename(bedpe_file_name).split(".")[0]
	if os.path.isfile(f".pickles/{base_file_name}_chromosomes.npy"):
		if verbose:
			print(f"pickle file for {base_file_name} exists")
		file = np.load(f".pickles/{base_file_name}_chromosomes.npy", allow_pickle=True)
	else:
		if verbose:
			print(f"pickle file for {base_file_name} doesn't exist")
		write_file(bedpe_file_name)
		file = np.load(f".pickles/{base_file_name}_chromosomes.npy", allow_pickle=True)
	
	return file


def read_bed(bed_filename):
	bed_file = pd.read_csv(bed_filename, delimiter=",")
	for i in ["chr", "start", "end"]:
		if i not in bed_file.columns:
			print(f"{i} column is not in mutation file")
			print(bed_file.columns)
			sys.exit(0)
	cols = bed_file.columns.tolist()
	cols.remove("chr")
	cols.remove("start")
	cols.remove("end")
	cols = ["chr", "start", "end"] + cols
	bed_file = bed_file.astype({"chr": str})
	return bed_file[cols]


def initial_intervals(sorted_intervals, mutation):
	mutants = set()
	
	for i in sorted_intervals[(( sorted_intervals[:, 1] > mutation["start"]) \
			& ( sorted_intervals[:, 0] <= mutation["start"]) \
			| (( sorted_intervals[:, 1] > mutation["end"]) \
			& ( sorted_intervals[:, 0] <= mutation["end"]) ))]:
		mutants.add(i[2])
	
	return mutants


def read_driver_genes(filename):
	
	if filename == None:
		return None
	
	genes = pd.read_csv(filename)
	genes.Chromosome = genes.Chromosome.str.replace("chr", "")
	
	return genes


def find_driver_overlaps(genes, intervals):
	
	gene_interval = {}
	for geneType in genes.gene_type.unique():
		gene_interval[geneType] = {}
	
	for geneType in genes.gene_type.unique():
		for index, row in genes.loc[genes.gene_type == geneType, :].iterrows():
			for interval in intervals[\
					((intervals[:, 0] > row["Start"]) & \
						(intervals[:, 0] < row["End"])) |\
					((intervals[:, 1] > row["Start"]) &\
						(intervals[:, 1] < row["End"])) ]:
				
				if interval[2] in gene_interval[geneType]:
					gene_interval[geneType][interval[2]].append( row["gene_name"] )
				else:
					gene_interval[geneType][interval[2]] = [ row["gene_name"] ]
					
	return gene_interval


def counter(array, initials, scores, gene_interval, genes, mutation):
	# if initials set is empty, there is no overlap at all
	if len(initials) == 0:
		return None
	
	# assing initial interval to visited intervals
	visited_intervals = set( initials )
	indices = set( [ -1 * i for i in visited_intervals ] )
	
	G= nx.Graph()
	for i in initials:
		G.add_edge( i, -i )
		G.add_edge( 0, i)
	
	score_counter = 0
	
	visited_intervals = visited_intervals | indices
	overlaps = set()
	
	while True:		
		# working in TRUBA with tuple(indices) instead of *indices
		overlaps = set.union( *array[[ *indices ]] )
		for i in indices:
			for j in array[i]:
				G.add_edge(i, j)
			G.add_edge(i, -i)
		
		indices = set( [-i for i in overlaps] )
		indices = indices.difference(visited_intervals)
		
		if not indices:
			break
		
		visited_intervals = visited_intervals | indices | overlaps
	
	interaction_counter = (G.number_of_nodes() - 1)//2
	overlap_counter = G.number_of_edges() - interaction_counter - len(initials)
	cycle_counter = G.number_of_edges() - G.number_of_nodes() + 1
	
	for i in G.nodes:
		if i > 0:
			score_counter += scores[i-1]
	
	result = [ G.number_of_nodes() - 1, interaction_counter, overlap_counter, \
				cycle_counter, round(score_counter/interaction_counter,2)]
	
	if gene_interval == None:
		return 	result
	
	for ind, mut in genes.loc[ ( (genes.Start <= mutation["start"]) &\
			(genes.End >= mutation["start"]) ) |\
			( (genes.Start <= mutation["end"]) &\
			(genes.End >= mutation["end"]) ),:].iterrows():
		#gene_interval[ mut["gene_type"] ][0] = mut["gene_name"]
		if 0 in gene_interval[ mut["gene_type"] ]:
			gene_interval[ mut["gene_type"] ][0].append( mut["gene_name"] )
		else:
			gene_interval[ mut["gene_type"] ][0] = [ mut["gene_name"] ]
	
	paths = nx.single_source_dijkstra_path_length(G, 0, cutoff= int(ranges[-1]) )
	for geneType in gene_interval:
		inv_ranges = len(ranges)*[0]
		inv_ranges[-1] = len(paths.keys() & gene_interval[geneType].keys())
		set_ranges = len(ranges) * [None]
		
		for i in range(len(ranges)):
			set_ranges[i] = set()
		
		for  i in paths.keys() & gene_interval[geneType].keys():
			for j in range(len(ranges)-1):
				if paths[i] <= int(ranges[j]):
					inv_ranges[j] += 1
					set_ranges[j] = set_ranges[j] | set(gene_interval[geneType][i])
			set_ranges[-1] = set_ranges[-1] | set(gene_interval[geneType][i])
		gene_ranges = len(ranges) * [0]
		
		for i in range(len(ranges)):
			gene_ranges[i] = len(set_ranges[i])
		
		for i in range(len(ranges)):
			result.append(inv_ranges[i])
			result.append(gene_ranges[i])
			result.append("|".join(list(set_ranges[i])))
		
		#if 0 in gene_interval[geneType]:
		gene_interval[geneType].pop(0, None)
	return result


def worker(bedpe_filename, bed_filenames):
	base_bedpe_name = os.path.basename(bedpe_filename).split(".")[0]
	bedpe_chromosomes = check_pickle_file(bedpe_filename)
	if only_write:
		return 0
	for bed_filename in bed_filenames:
		t1 = time.time()
		base_bed_name = os.path.basename(bed_filename).split(".")[0]
		bed_file = read_bed(bed_filename)
		total_mutation = bed_file.index.stop
		n = 0
		bed_chromosomes = bed_file.loc[:, "chr"].values
		
		result_file = bed_file.copy()
		columns =["intervals", "interactions", "overlaps", \
				"cycle_rank", "score"]
		
		#ranges = ["5", "10", "20", "30"]
		if genes is not None:
			for i in genes["gene_type"].unique():
				for j in ranges:
					columns.append(f"{i}_range_{j}_NumOfInv")
					columns.append(f"{i}_range_{j}_NumOfGen")
					columns.append(f"{i}_range_{j}_NamOfGen")
		
		# adding columns to result file
		for i in columns:
			if "score" == i:
				result_file[i] = 0.0
			elif i.endswith("NamOfGen"):
				result_file[i] = ""
			else:
				result_file[i] = 0
		for common_chromosome in np.intersect1d(bed_chromosomes, bedpe_chromosomes):
			mutations = bed_file.loc[ bed_file.chr == common_chromosome , ["start", "end"]]
			with open(f".pickles/{base_bedpe_name}_{common_chromosome}_intervals.pickle", "rb") as f:
				intervals = pickle.load(f)
			with open(f".pickles/{base_bedpe_name}_{common_chromosome}_array.pickle", "rb") as f:
				array = pickle.load(f)
			with open(f".pickles/{base_bedpe_name}_{common_chromosome}_scores.pickle", "rb") as f:
				scores = pickle.load(f)
			
			if genes is None:
				f_genes = None
				gene_interval = None
			else:
				f_genes = genes.loc[genes.Chromosome == common_chromosome, :]
				gene_interval = \
						find_driver_overlaps(f_genes, intervals)
			
			for index, mutation in mutations.iterrows():
				initials = initial_intervals(intervals, mutation)
				count_values = counter(array, initials, scores,\
						gene_interval, f_genes, mutation)
				if count_values is not None:
					result_file.loc[ (result_file.chr == common_chromosome) & \
							(result_file.start == mutation.start), columns] = count_values
				
				if verbose:
					t2 = time.time()
					print(f"for {base_bedpe_name}_{base_bed_name}_result.csv:", end=" ")
					sys.stdout.write(f"{n}/{total_mutation}, time: {round(t2-t1, 2)} seconds\r")
					sys.stdout.flush()
					n += 1
		
		t2 = time.time()
		biosample = metadata.loc[base_bedpe_name, "Biosample term name"]
		target = metadata.loc[base_bedpe_name, "Experiment target"]
		result_file.to_csv(f"{output_dir}/{base_bedpe_name}_{base_bed_name}.csv",\
				index=False)
		if verbose:
			print(f"{output_dir}/{base_bedpe_name}_{base_bed_name}.csv has been written")
			print(f"{output_dir}/{base_bedpe_name}_{base_bed_name}.csv took {t2-t1} second(s)")
			
	return 0


def main():
	parser = argparse.ArgumentParser(description="InvMutMapper.py script")
	
	parser.add_argument('--bed_files',\
			required = True, nargs='+', help="input mutation file(s)")
	parser.add_argument('-o', "--output",\
			nargs='?', default="result", help="output directory")
	parser.add_argument('--bedpe_files', \
			required = True, nargs='+', help="input bedpe file(s)")
	parser.add_argument('--debug', type=bool,\
			nargs='?', default = False, help="Debug mode (default false)")
	parser.add_argument("-ow", "--only_write", action="store_true",\
			help="If True, .pickle files are writen and stop")
	parser.add_argument("-v", "--verbose", action="store_true")
	parser.add_argument("-dg", "--drivergenes", help="oncogenes file",\
			nargs="?", default=None)
	parser.add_argument("-md", "--metadata", help="metadata for bedpe files",\
			required = True)
	parser.add_argument("-r", "--remove_pickles", action="store_true",\
			help="If True, .pickle files will be removed after calulation finished")
	parser.add_argument("-s", "--serial", action="store_true",\
			help="If True, code will be run in serial")
	parser.add_argument("--ranges", \
			help="custom ranges (shortest path from mutation) should be greater than 0 integers",\
			nargs="+", default = ["0", "5", "10"])
	
	args = parser.parse_args()
	
	print("output directory:", end=" ")
	print(args.output)
	print("debug mode:", end=" ")
	print(args.debug)
	print("only-write mode:", end=" ")
	print(args.only_write)
	print("mode:", end=" ")
	if args.serial:
		print("serial")
	else:
		print("parallel")
	global output_dir, only_write, verbose, genes, metadata, ranges
	genes = read_driver_genes(args.drivergenes)
	metadata = pd.read_csv(args.metadata, sep="\t")
	metadata = metadata.loc[:, ["File accession", "Biosample term name", "Experiment target"]]
	metadata["Experiment target"] = metadata["Experiment target"].str.replace("-human", "")
	metadata["Biosample term name"] = metadata["Biosample term name"].str.replace(" ", "_")
	metadata.set_index("File accession", inplace=True)
	
	only_write = args.only_write
	verbose = args.verbose
	
	for i in args.ranges:
		if (not i.isnumeric()):
			print("not numeric or less than 0 or not an integer")
			sys.exit()
	ranges = sorted(list(set(map(int, args.ranges))))
	ranges = list(map(str, ranges))
	
	# create output directory
	current_dir = os.path.abspath(os.getcwd())
	output_dir = os.path.join(current_dir, args.output)
	tmp_dir = os.path.join(current_dir, ".pickles")
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	if not os.path.isdir(tmp_dir):
		os.makedirs(tmp_dir)
	
	if args.serial:
		for bedpe_file_name in args.bedpe_files:
			worker(bedpe_file_name, args.bed_files)
	else:
		pool = Pool(processes=os.cpu_count())
		jobs = []
		for bedpe_file_name in args.bedpe_files:
			jobs.append(pool.apply_async(worker, \
					args = (bedpe_file_name, args.bed_files)))
		pool.close()
		pool.join()
	
	if args.remove_pickles:
		for i in os.listdir(tmp_dir):
			os.remove(os.path.join(tmp_dir, i))
			if verbose:
				print(f"{os.path.join(tmp_dir, i)} removed")
	
	print("Finished")

if __name__ == "__main__":
	main()
