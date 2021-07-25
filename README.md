# Sigmap
Sigmap is a streaming method for mapping raw nanopore signal to reference genomes. 

## Installation
First get the repo with submodules (make sure you use `--recursive`):
```
git clone --recursive git@github.com:haowenz/sigmap.git
```
Make sure you have GCC version > 8. Then just run:
```
cd sigmap && make
```

## Usage
First, an index needs to be built for the reference using ONT pore models (as a submodule in the extern folder). We use yeast genome as an example:
```
./sigmap -i -r yeast.fasta -p extern/kmer_models/r9.4_180mv_450bps_6mer/template_median68pA.model -o yeast_index
```
It will generate the index file *yeast_index.si*. Note that a genome point cloud file *yeast_index.pt* will also be saved. But it can also be generated very quickly on the fly every time before mapping. After index construction, yeast raw signals in fast5 format can be mapped using
```
./sigmap -m -r yeast.fasta -p extern/kmer_models/r9.4_180mv_450bps_6mer/template_median68pA.model -x yeast_index -s /path/to/yeast/fast5/dir -o yeast_mapping.paf -t 4
```
This command map all the fast5 reads in `/path/to/yeast/fast5/dir` and its subdirectories to the yeast genome using 4 threads. The output will be saved to `yeast_mapping.paf` in a modified PAF format used by [Uncalled](https://github.com/skovaka/UNCALLED). 

Many other parameters can be found in the help information:
```
./sigmap -h
```

It is possible that your reads are compressed with the [VBZ compression](https://github.com/nanoporetech/vbz_compression) from Nanopore. Then you have to download the proper HDF5 plugin from [here](https://github.com/nanoporetech/vbz_compression/releases) and make sure it can be found by your HDF5 library:
```
export HDF5_PLUGIN_PATH=/path/to/hdf5/plugins/lib
```

## Getting help
If you encounter bugs or have further questions or requests, you can raise an issue at the [issue page](https://github.com/haowenz/sigmap/issues) or contact me by email hwzhang@gatech.edu.

## Citing Sigmap
If you use Sigmap, please cite:

Zhang, H., Li, H., Jain, C., Cheng, H., Au, K. F., Li, H., & Aluru, S. (2021). Real-time mapping of nanopore raw signals. Bioinformatics, 37(Supplement_1), i477-i483. https://doi.org/10.1093/bioinformatics/btab264
