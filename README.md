# MutationNetwork.py

## Installation

Script can be install via github directly with example input data.

### conda env

Before creating conda env, some channels should be added. ***It is going to change your config file!***

``` 
conda config --add channels conda-forge

conda config --add channels bioconda 
```

Now, the code above will create a conda env and install necessary packages.

``` 
conda create --name mutation_network python=3.10

conda install --name mutation_network --file requirements.txt

conda activate mutation_network
```

## Usage

usage: MutationNetwork.py [-h] --vcf\_files VCF\_FILES [VCF\_FILES ...] [-o [OUTPUT]] --bedpe\_files BEDPE\_FILES [BEDPE\_FILES ...] [-ow] [-v] [--genes [GENES]] [-r] [--ranges [RANGES]] 
[--output\_format {count,binary} [{count,binary} ...]] [-pb | -pv | -sb | -sv]

example: 

```
python MutationNetwork.py --vcf_files mutations.vcf -o output --bedpe_files ENCFF597SQA.bedpe.gz --genes gencode.v47.basic.annotation.gtf.gz --ranges "range(11)" -v -sv --output_format binary

python MutationNetwork.py --vcf_files mutations.vcf -o output --bedpe_files ENCFF597SQA.bedpe.gz --genes gencode.v47.basic.annotation.gtf.gz --ranges "[0, 1, 2, 3, 10]" -pb --output_format binary count

```

## Input Files

gencode (\*.annotation.gtf.gz) -- gencode file that has been downloaded from <https://www.gencodegenes.org/human/> (content: "Basic gene annotation", Regions: "CHR", Dowload: "GTF").

loop file (.bedpe) -- loop files store long range chromosomal interactions and can be retrived from <https://www.encodeproject.org/>

## Result file: binary matrix for each range

Each row represents a mutation sample file and each column represents genes given in the annotation file.
If correspoding sample file and gene is 1, it represents for the given range, at least one mutation in the sample reaches the gene via long range interactions.

## Method

### Step 1: processing bedpe

![image](img/method_explanation.png)

### Step 2: flowchart

<p align="center">

<img src="img/flowchart.png"/>

</p>

### Others

![image](img/how_tree_to_array_works.png)

![image](img/multiple_interaction.drawio.png)

![image](img/mutation_effect.png)


