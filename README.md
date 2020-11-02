# Sigmap
Sigmap is a streaming method for mapping raw nanopore signal to reference genomes. 

## Installation
First get the repo with submodules (make sure you use `--recursive`):
```
git clone --recursive git@github.com:haowenz/sigmap.git
```
Then just run:
```
cd sigmap && make
```

## Usage
First, an index needs to be built for the reference using ONT pore models (as a submodule in the extern folder). We use yeast genome as an example:
```
./sigmap -i -r yeast.fasta -p extern/kmer_models/r9.4_180mv_450bps_6mer/template_median68pA.model -o yeast_index
```
It will generate two index files *yeast_index.pt* and *yeast_index.si*. After index construction, yeast raw signals in fast5 format can be mapped using
```
./sigmap -m -r yeast.fasta -p extern/kmer_models/r9.4_180mv_450bps_6mer/template_median68pA.model -x yeast_index -s /path/to/yeast/fast5/dir -o yeast_mapping.paf -t 4
```
This command map all the fast5 reads in `/path/to/yeast/fast5/dir` and its subdirectories to the yeast genome using 4 threads. The output will be saved to `yeast_mapping.paf` in a modified PAF format used by [Uncalled](https://github.com/skovaka/UNCALLED). 

Many other parameters can be found in the help information:
```
./sigmap -h
```
