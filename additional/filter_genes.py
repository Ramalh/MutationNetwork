import pyranges as pr
import pandas as pd
import sys

gene_types = sys.argv[1:]
print(gene_types)

data = pr.read_gtf("gencode.v47.basic.annotation.gtf.gz")

data1 = data[["Source", "Chromosome", "Start", "End", "gene_id", "gene_name", "gene_type", "hgnc_id"]][(data.Feature=="gene") & (data.gene_type==gene_types[0])].df

if len(gene_types) == 1:
	data1.to_csv(f"genes_{gene_types[0]}.csv", index=False)
	sys.exit(0)

for gene_type in gene_types[1:]:
	data2 = data[["Source", "Chromosome", "Start", "End", "gene_id", "gene_name", "gene_type", "hgnc_id"]][(data.Feature=="gene") & (data.gene_type==gene_type)].df
	data1 = pd.concat([data1, data2])

data1.to_csv(f"genes_{'_'.join(gene_types)}.csv", index=False)
