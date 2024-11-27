
import pyranges as pr
import pandas as pd
import argparse
import sys


def main():
	# start of argument parser
	parser = argparse.ArgumentParser(description="InvMutMapper.py script")
	parser.add_argument("-n", "--names", required = True,\
			help="input hgnc names")
	parser.add_argument("-t", "--type", required = True,\
			help="gene types")
	parser.add_argument("-g", "--gencode", required = True,\
			help="gencode file")
	parser.add_argument("-o", "--output", nargs='?',\
			default="output", help="output filename")
	
	args = parser.parse_args()
	gene_names_filename = args.names
	gene_names_type = args.type
	genes_filename = args.gencode
	output = args.output
	# end of argument parser
	
	gene_names = pd.read_csv(gene_names_filename, sep="\t").loc[:, "Hugo Symbol"]
	genes = pr.read_gtf(genes_filename)
	
	gene = genes[["Source", "Chromosome", "Start", "End", "gene_id", "gene_name", "gene_type", "hgnc_id"]][genes.Feature=="gene"]
	
	selected_gene =  gene[gene.gene_name.isin(gene_names)].df
	selected_gene.loc[:, "gene_type"] = gene_names_type
	selected_gene.hgnc_id =  selected_gene.hgnc_id.str[5:]
	selected_gene.to_csv(f"{output}.csv", index=False)


if __name__ == "__main__":
	main()
