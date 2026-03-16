#!/usr/bin/python3

from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Pool
import networkx as nx
import pyranges as pr
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


def read_vcf(vcf_filename):
	vcf_file = pd.read_csv(vcf_filename, comment = "#", sep = "\t", usecols=[0, 1],\
			names = ["chr", "pos"], dtype = {"chr": str, "pos": int})
	vcf_file["chr"] = vcf_file["chr"].str.replace("chr", "")
	return vcf_file
	for i in ["chr", "start", "end"]: #change here
		if i not in bed_file.columns:
			print(f"{i} column is not in mutation file")
			print(bed_file.columns)
	cols = bed_file.columns.tolist()
	cols.remove("chr")
	cols.remove("start")
	cols.remove("end")
	cols = ["chr", "start", "end"] + cols
	bed_file = bed_file.astype({"chr": str})
	return bed_file[cols]


def initial_intervals(sorted_intervals, mutation):
	mutants = set()
	
	for i in sorted_intervals[ (mutation > sorted_intervals[:, 0]) \
			& (mutation < sorted_intervals[:, 1] )]:
		mutants.add( i[2] )
	
	return mutants


def read_genes(gtf_filename):
	genes_gtf = pr.read_gtf(gtf_filename)
	
	genes = genes_gtf[["Chromosome", "Start", "End", "Strand", "gene_name"]][(genes_gtf.Feature=="gene")].df
	
	genes.drop_duplicates(subset=["gene_name"], ignore_index=True, inplace=True)
	genes.Chromosome = genes.Chromosome.str.replace("chr", "")
	
	genes["gene_number"] = genes.index
	
	return genes


def find_driver_overlaps(genes, intervals):
	
	gene_interval = dict()
	
	np_genes = genes[["Start", "End", "gene_number"]].copy()
	np_genes["typ"] = 1
	np_genes = np_genes.values
	intervals = np.c_[intervals, np.zeros( intervals.shape[0] ) ]
	array = np.concatenate( (intervals, np_genes), axis = 0)
	array = array[array[:, 0].argsort()]
	array = array.astype(int)
	
	lt = len( array )
	i = -1
	for inv in array:
		i += 1
		j = i + 1
		while j < lt and inv[1] > array[j][0] :
			
			if array[j, 3] and not inv[3]:
				if inv[2] in gene_interval:
					gene_interval[inv[2]].add( array[j, 2] )
				else:
					gene_interval[inv[2]] = { array[j, 2] }
			if inv[3] and not array[j, 3]:	
				if array[j, 2] in gene_interval:
					gene_interval[array[j, 2]].add( inv[2] )
				else:
					gene_interval[array[j, 2]] = { inv[2] }
			j += 1
		
	gene_interval[0] = set()
	
	return gene_interval


def counter(array, initials, scores, gene_interval, genes_array):

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
	
	c = 0
	max_range = max(ranges)
	while c <= max_range:
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
	
	if gene_interval == None:
		return result
	
	paths = nx.single_source_dijkstra_path_length(G, 0, cutoff= max_range )
	
	set_ranges = [set()] * (max_range + 1)
	
	for  i in paths.keys() & gene_interval.keys():
		set_ranges[ paths[i] ] = set_ranges[ paths[i] ] | gene_interval[i]
	
	for i in range(len(ranges)):
		for j in range(ranges[i]+1 ):
			genes_array[ i , list(set_ranges[j]) ] = 1


def workerParallelVCF(bedpe_filenames, vcf_filename):
	
	base_vcf_name = os.path.basename(vcf_filename).split(".")[0]
	vcf_file = read_vcf(vcf_filename)
	vcf_chromosomes = vcf_file.loc[:, "chr"].values
	
	if only_write:
		return 0
	
	for bedpe_filename in bedpe_filenames:
		t11 = time.time()
		
		genes_array = np.zeros(( len(ranges) , len(genes.gene_name) ), dtype = int)
		base_bedpe_name = os.path.basename(bedpe_filename).split(".")[0]
		bedpe_chromosomes = check_pickle_file(bedpe_filename)
		
		for common_chromosome in np.intersect1d(vcf_chromosomes, bedpe_chromosomes):
			t10 = time.time()
			mutations = vcf_file.loc[ vcf_file.chr == common_chromosome , "pos"].reset_index(drop=True)
			with open(f".pickles/{base_bedpe_name}_{common_chromosome}_intervals.pickle", "rb") as f:
				intervals = pickle.load(f)
			with open(f".pickles/{base_bedpe_name}_{common_chromosome}_array.pickle", "rb") as f:
				array = pickle.load(f)
			with open(f".pickles/{base_bedpe_name}_{common_chromosome}_scores.pickle", "rb") as f:
				scores = pickle.load(f)
					
			f_genes = genes.loc[genes.Chromosome == common_chromosome, :]
			gene_interval = find_driver_overlaps(f_genes, intervals)
			range_0_genes = set()
			
			for index, mutation in mutations.items():
				for ind, mut in f_genes.loc[ ( ( mutation > f_genes.Start) &\
					( mutation < f_genes.End ) ),:].iterrows():
					range_0_genes.add(mut["gene_number"])
		
			genes_array[:, list(range_0_genes)] = 1
			for index, mutation in mutations.items():
				
				initials = initial_intervals(intervals, mutation)
				if len(initials) == 0:
					continue
				counter(array, initials, scores,\
						gene_interval, genes_array)
		
		t12 = time.time()
		print(f"{base_vcf_name} took {t12 - t11: .2f} sec")
		for i in range(len(ranges)):
			output_file = f"{output_dir}/gene_similarity_range_{ranges[i]}_all.csv"
			
			# Build the full row as a string
			row_data = ",".join(map(str, genes_array[i, :]))
			line = f"{base_vcf_name},,{row_data}\n"
			
			with open(output_file, "a") as f:
			    f.write(line)



def main():
	t1 = time.time()
	
	parser = argparse.ArgumentParser(description="InvMutMapper.py script")
	group = parser.add_mutually_exclusive_group()
	
	parser.add_argument('--vcf_files',\
			required = True, nargs='+', help="input mutation file(s)")
	parser.add_argument('-o', "--output",\
			nargs='?', default="result", help="output directory")
	parser.add_argument('--bedpe_files', \
			required = True, nargs='+', help="input bedpe file(s)")
	parser.add_argument("-ow", "--only_write", action="store_true",\
			help="If True, .pickle files are writen and stop")
	parser.add_argument("-v", "--verbose", action="store_true")
	parser.add_argument("--genes", help="gene annotation file (gtf)",\
			nargs="?", default=None)
	parser.add_argument("-r", "--remove_pickles", action="store_true",\
			help="If True, pickle files will be removed after calulation finished")
	parser.add_argument("--ranges", \
			help="custom ranges (shortest path from mutation) should be greater than 0 integers",\
			nargs="?", default = "range(11)")
	group.add_argument("-pb", "--parallelBEDPE", help="set parallel mode \
			for bedpe files (defaul is parallelVCF)", action="store_true")
	group.add_argument("-pv", "--parallelVCF", help="set parallel mode \
			for vcf files (defaul is parallelVCF)", action="store_true")
	group.add_argument("-sb", "--serialBEDPE", help="set serial mode \
			for bedpe files (defaul is parallelVCF)", action="store_true")
	group.add_argument("-sv", "--serialVCF", help="set serial mode \
			for vcf files (defaul is parallelVCF)", action="store_true")
	
	args = parser.parse_args()
	print("output directory:", end=" ")
	print(args.output)
	print("only-write mode:", end=" ")
	print(args.only_write)
	print("mode:", end=" ")
	
	if not (args.parallelBEDPE or args.parallelVCF or args.serialBEDPE or args.serialVCF):
		args.parallelVCF = True
	if args.parallelBEDPE:
		print("parallelVCF")
	elif args.parallelVCF:
		print("parallelVCF")
	elif args.serialBEDPE:
		print("serialBEDPE")
	else:
		print("serialVCF")
	
	global output_dir, only_write, verbose, genes, ranges
	genes = read_genes(args.genes)
	
	only_write = args.only_write
	verbose = args.verbose
	
	ranges = eval(args.ranges)
	ranges = sorted(list(set(map(int, ranges))))
	
	# create output directory
	current_dir = os.path.abspath(os.getcwd())
	output_dir = os.path.join(current_dir, args.output)
	tmp_dir = os.path.join(current_dir, ".pickles")
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	if not os.path.isdir(tmp_dir):
		os.makedirs(tmp_dir)
	
	headers = np.array( ["ID", "subtype", *genes.gene_name] )
	for i in range(len(ranges)):
		np.savetxt( f"{output_dir}/gene_similarity_range_{ranges[i]}_all.csv", \
				[ headers ], delimiter = ',', fmt='%s') 
	
	args.vcf_files = sorted(args.vcf_files, key = lambda x: os.path.getsize(x), reverse=True )	
	
	if args.serialBEDPE:
		for bedpe_file_name in args.bedpe_files:
				workerParallelBedpe(bedpe_file_name, args.vcf_files)
	
	elif args.parallelBEDPE:
		with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
			features = [ executor.submit(workerParallelBedpe, bedpe_file_name, \
					args.vcf_files) for bedpe_file_name in args.bedpe_files]
			for feature in as_completed(features):
				feature.result()
			
	elif args.serialVCF:
		for vcf_file in args.vcf_files:
			workerParallelVCF(args.bedpe_files, vcf_file)
			
	elif args.parallelVCF:
		with ProcessPoolExecutor(max_workers= os.cpu_count() ) as executor:
			features = [ executor.submit(workerParallelVCF, args.bedpe_files, \
					vcf_file) for vcf_file in args.vcf_files]
			for feature in as_completed(features):
				feature.result()		
	
	if args.remove_pickles:
		for i in os.listdir(tmp_dir):
			os.remove(os.path.join(tmp_dir, i))
			if verbose:
				print(f"{os.path.join(tmp_dir, i)} removed")
	
	t2 = time.time()
	print(f"Finished: {t2-t1}")



if __name__ == "__main__":
	main()

