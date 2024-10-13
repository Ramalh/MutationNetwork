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
	bed_file = pd.read_csv(bed_filename, delimiter=",",\
						names = ["chr", "start", "end", "start_hg19", "driver"],\
						header = None, dtype = \
						{0: object, 1: int, 2: int, 3: object},\
						skiprows=1)
	
	return bed_file


def initial_intervals(sorted_intervals, mutation):
	mutants = set()
	
	for i in sorted_intervals[( sorted_intervals[:, 1] > mutation) \
			& ( sorted_intervals[:, 0] <= mutation)]:
		mutants.add(i[2])
	
	return mutants


def read_driver_genes(filename):
	
	if filename == None:
		return None
	
	genes = pd.read_csv(filename,
		dtype = {0: object, 1: int, 2: object, 3: object, 4: int, 5: int})
	#genes.chromosome = "chr" + genes.chromosome
	genes.chromosome = genes.chromosome.str.replace("chr", "")
	genes_array = genes.iloc[:,[0, 2, 3, 4, 5]].values
	
	return genes_array


def find_driver_overlaps(genes, intervals):
	
	onco_interval_id = dict()
	ts_interval_id = dict()
	for i in genes[genes[:,0] == "oncogene", :]:
		for j in intervals[ ((i[3] < intervals[:, 0]) & (intervals[:, 0] < i[4])) \
				| ((i[3] < intervals[:, 1]) & (intervals[:, 1] < i[4])) ]:
			onco_interval_id[j[2]] = i[1]
	
	for i in genes[genes[:,0] == "tumorSuppressor", :]:
		for j in intervals[ ((i[3] < intervals[:, 0]) & (intervals[:, 0] < i[4])) \
				| ((i[3] < intervals[:, 1]) & (intervals[:, 1] < i[4])) ]:
			ts_interval_id[j[2]] = i[1]
	
	return onco_interval_id, ts_interval_id


def counter(array, initials, scores, onco_interval_id, ts_interval_id):
	# if initials set is empty, there is no overlap at all
	if len(initials) == 0:
		if onco_interval_id == None:
			return [0, 0, 0, 0, 0.0]
		else:
			return [0, 0, 0, 0, 0.0, 0, 0, 0, 0, 0, 0, \
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	
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
	
	if onco_interval_id == None:
		return [ G.number_of_nodes() - 1, interaction_counter, overlap_counter, \
				cycle_counter, round(score_counter/interaction_counter,2)]
	else:		
		paths = nx.single_source_dijkstra_path_length(G, 0, cutoff=30)
		onco_inv_5 = 0
		onco_inv_10 = 0
		onco_inv_20 = 0
		onco_inv_30 = len(paths.keys() & onco_interval_id.keys())
		onco_5_set = set()
		onco_10_set = set()
		onco_20_set = set()
		onco_30_set = set()
		for i in paths.keys() & onco_interval_id.keys():
			if paths[i] <= 5:
				onco_inv_5 += 1
				onco_5_set.add(onco_interval_id[i])
			if paths[i] <= 10:
				onco_inv_10 += 1
				onco_10_set.add(onco_interval_id[i])
			if paths[i] <= 20:
				onco_inv_20 += 1
				onco_20_set.add(onco_interval_id[i])
			onco_30_set.add(onco_interval_id[i])
		onco_gene_5 = len(onco_5_set)
		onco_gene_10 = len(onco_10_set)
		onco_gene_20 = len(onco_20_set)
		onco_gene_30 = len(onco_30_set)
		
		ts_inv_5 = 0
		ts_inv_10 = 0
		ts_inv_20 = 0
		ts_inv_30 = len(paths.keys() & onco_interval_id.keys())
		ts_5_set = set()
		ts_10_set = set()
		ts_20_set = set()
		ts_30_set = set()
		for i in paths.keys() & ts_interval_id.keys():
			if paths[i] <= 5:
				ts_inv_5 += 1
				ts_5_set.add(ts_interval_id[i])
			if paths[i] <= 10:
				ts_inv_10 += 1
				ts_10_set.add(ts_interval_id[i])
			if paths[i] <= 20:
				ts_inv_20 += 1
				ts_20_set.add(ts_interval_id[i])
			ts_30_set.add(ts_interval_id[i])
		ts_gene_5 = len(ts_5_set)
		ts_gene_10 = len(ts_10_set)
		ts_gene_20 = len(ts_20_set)
		ts_gene_30 = len(ts_30_set)
		
		return [ G.number_of_nodes() - 1, interaction_counter, overlap_counter, \
				cycle_counter, round(score_counter/interaction_counter,2), \
				onco_inv_5, onco_gene_5, "|".join(list(onco_5_set)), onco_inv_10, \
				onco_gene_10, "|".join(list(onco_10_set)), onco_inv_20, onco_gene_20,\
				"|".join(list(onco_20_set)), onco_inv_30, onco_gene_30, \
				"|".join(list(onco_30_set)), ts_inv_5, ts_gene_5, \
				"|".join(list(ts_5_set)), ts_inv_10, ts_gene_10, \
				"|".join(list(ts_10_set)), ts_inv_20, ts_gene_20, \
				"|".join(list(ts_20_set)), ts_inv_30, ts_gene_30, \
				"|".join(list(ts_30_set))]

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
		
		result_file = bed_file.loc[:, ["chr", "start", "end", "start_hg19", "driver"]].copy()
		columns =["intervals", "interactions", "overlaps", \
				"cycles", "score", "onco_range_5", "onco_range_5_gene",\
				"onco_range_5_list", "onco_range_10", "onco_range_10_gene",\
				"onco_range_10_list", "onco_range_20", "onco_range_20_gene",\
				"onco_range_20_list", "onco_range_30", "onco_range_30_gene",\
				"onco_range_30_list", "ts_range_5", "ts_range_5_gene", \
				"ts_range_5_list", "ts_range_10", "ts_range_10_gene", \
				"ts_range_10_list", "ts_range_20", "ts_range_20_gene",\
				"ts_range_20_list", "ts_range_30", "ts_range_30_gene", \
				"ts_range_30_list"]
		
		if genes is None:
			columns = columns[:5]
		
		# adding columns to result file
		for i in columns:
			if "score" == i:
				result_file[i] = 0.0
			elif i.endswith("list"):
				result_file[i] = ""
			else:
				result_file[i] = 0
	
		for common_chromosome in np.intersect1d(bed_chromosomes, bedpe_chromosomes):
			mutations = bed_file.loc[ bed_file.chr == common_chromosome , "start"]
			with open(f".pickles/{base_bedpe_name}_{common_chromosome}_intervals.pickle", "rb") as f:
				intervals = pickle.load(f)
			with open(f".pickles/{base_bedpe_name}_{common_chromosome}_array.pickle", "rb") as f:
				array = pickle.load(f)
			with open(f".pickles/{base_bedpe_name}_{common_chromosome}_scores.pickle", "rb") as f:
				scores = pickle.load(f)
			
			if genes is None:
				onco_interval_id, ts_interval_id = None, None
			else:
				filtered_genes = genes[genes[:, 2] == common_chromosome, :]
				onco_interval_id, ts_interval_id = find_driver_overlaps(filtered_genes, intervals)
			
			for mutation in mutations:
				initials = initial_intervals(intervals, mutation)
				count_values = counter(array, initials, scores, \
						onco_interval_id, ts_interval_id)
				
				result_file.loc[ (result_file.chr == common_chromosome) & \
						(result_file.start == mutation), columns] = count_values
				
				if verbose:
					t2 = time.time()
					print(f"for {base_bedpe_name}_{base_bed_name}_result.csv:", end=" ")
					sys.stdout.write(f"{n}/{total_mutation}, time: {round(t2-t1, 2)} seconds\r")
					sys.stdout.flush()
					n += 1
		
		t2 = time.time()
		biosample = metadata.loc[base_bedpe_name, "Biosample term name"]
		target = metadata.loc[base_bedpe_name, "Experiment target"]
		result_file.to_csv(f"{output_dir}/{base_bedpe_name}_{biosample}_{target}.csv",\
				index=False)
		if verbose:
			print(f"{output_dir}/{base_bedpe_name}_{biosample}_{target}.csv has been written")
			print(f"{output_dir}/{base_bedpe_name}_{biosample}_{target}.csv took {t2-t1} second(s)")
			
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
	
	args = parser.parse_args()
	
	print("output directory:", end=" ")
	print(args.output)
	print("debug mode:", end=" ")
	print(args.debug)
	print("only-write mode:", end=" ")
	print(args.only_write)
	
	global output_dir, only_write, verbose, genes, metadata
	genes = read_driver_genes(args.drivergenes)
	metadata = pd.read_csv(args.metadata, sep="\t")
	metadata = metadata.loc[:, ["File accession", "Biosample term name", "Experiment target"]]
	metadata["Experiment target"] = metadata["Experiment target"].str.replace("-human", "")
	metadata["Biosample term name"] = metadata["Biosample term name"].str.replace(" ", "_")
	metadata.set_index("File accession", inplace=True)
	
	only_write = args.only_write
	verbose = args.verbose
	
	# create output directory
	current_dir = os.path.abspath(os.getcwd())
	output_dir = os.path.join(current_dir, args.output)
	tmp_dir = os.path.join(current_dir, ".pickles")
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	if not os.path.isdir(tmp_dir):
		os.makedirs(tmp_dir)
	
	for bedpe_file_name in args.bedpe_files:
		worker(bedpe_file_name, args.bed_files)
	
	#pool = Pool(processes=os.cpu_count())
	#jobs = []
	#for bedpe_file_name in args.bedpe_files:
	#	jobs.append(pool.apply_async(worker, \
	#			args = (bedpe_file_name, args.bed_files)))
	#pool.close()
	#pool.join()

	if args.remove_pickles:
		for i in os.listdir(tmp_dir):
			os.remove(os.path.join(tmp_dir, i))
			if verbose:
				print(f"{os.path.join(tmp_dir, i)} removed")
	
	print("Finished")

if __name__ == "__main__":
	main()
