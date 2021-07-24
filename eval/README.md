# Evaluation
We provide the evaluation script to reproduce the results in the paper.

## Usage
First, run minimap2 to get a groud truth for real sequencing data. Then run Uncalled and Sigmap with parameters specified in the paper and get the PAF outputs. Run [pafstats](https://github.com/skovaka/UNCALLED/blob/master/uncalled/pafstats.py) on the alignment outputs and get the read classification results in PAF. Finally, run our evaluation script and get the accuracy statistics. 

Note that pafstats only counts the mapping time of mapped reads and by default extend the mapping intervals by 1.5x. While in our evalution, we count the mapping time of all the reads in the groud truth and change pafstats to not extend the mapping intervals.
