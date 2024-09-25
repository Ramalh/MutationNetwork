# invMutMapper.py -- Interval Mutation Mapper

## Installation

Script can be install via github directly with example input data.

### conda env

conda env create -f invMutMapper.yml

## Usage

usage: python invMutMapper.py [-h] --bed\_files BED\_FILES [BED\_FILES ...] [-o [OUTPUT]] --bedpe\_files BEDPE\_FILES [BEDPE\_FILES ...] [--debug [DEBUG]] [-ow] [-v] [-r] [-dg [DRIVERGENES]] -md METADATA

Example: python InvMutMapper.py --bedpe\_files bedpe/\*.bedpe.gz --bed\_files mutations.csv -dg driverGenes.csv -md metadata.tsv -v

## Input Files

bedpe\_files -- 

bed\_file --

driver\_genes --

metadata --
