import os
import sys
import argparse
import subprocess
import fileinput

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

from statistics import median
from statistics import mean

if (len(sys.argv) < 4):
  print("usage: eval.py uncalled_ann.paf sigmap_ann.paf output_prefix")
  sys.exit(1)

fps = []
for tool in [1, 2]:
  fps.append(open(sys.argv[tool]))

uncalled_tp = 0
uncalled_fp = 0
uncalled_fn = 0
uncalled_tn = 0

uncalled_chunk = []
uncalled_time_per_read = []
for line in fps[0]:
  cols = line.rstrip().split()
  mt = float(cols[14].split(":")[2])
  if (cols[15].split(":")[2] != 'na'):
    uncalled_time_per_read.append(mt)
  if (cols[15].split(":")[2] == 'tp'):
    uncalled_tp += 1
  if (cols[15].split(":")[2] == 'fp'):
    uncalled_fp += 1
  if (cols[15].split(":")[2] == 'fn'):
    uncalled_fn += 1
  if (cols[15].split(":")[2] == 'tn'):
    uncalled_tn += 1
print("Uncalled TP: " + str(uncalled_tp))
print("Uncalled FP: " + str(uncalled_fp))
print("Uncalled FN: " + str(uncalled_fn))
print("Uncalled TN: " + str(uncalled_tn))
uncalled_preciosion = uncalled_tp / (uncalled_tp + uncalled_fp)
print("Uncalled precision: " + str(uncalled_preciosion))
uncalled_recall = uncalled_tp / (uncalled_tp + uncalled_fn)
print("Uncalled recall: " + str(uncalled_recall))
print("Uncalled F-1 score: " + str(2 * uncalled_preciosion * uncalled_recall / (uncalled_preciosion + uncalled_recall)))
print("Mean time per read : " + str(mean(uncalled_time_per_read)))
print("Median time per read : " + str(median(uncalled_time_per_read)))
print("#Done with uncalled\n")

sigmap_tp = 0
sigmap_fp = 0
sigmap_fn = 0
sigmap_tn = 0

sigmap_time_per_chunk = []
sigmap_time_per_read = []
for line in fps[1]:
  cols = line.rstrip().split()
  if (len(cols) == 24):
    mt = float(cols[12].split(":")[2])
    if (cols[23].split(":")[2] != 'na'):
      sigmap_time_per_read.append(mt)
    chunk = int(cols[13].split(":")[2])
    cm = int(cols[15].split(":")[2])
    nc = int(cols[16].split(":")[2])
    s1 = float(cols[17].split(":")[2])
    s2 = float(cols[18].split(":")[2])
    sm = float(cols[19].split(":")[2])
    ad = float(cols[20].split(":")[2])
    at = float(cols[21].split(":")[2])
    aq = float(cols[22].split(":")[2])
    if (cols[23].split(":")[2] == 'tp'):
      sigmap_tp += 1
      sigmap_time_per_chunk.append(mt / chunk)
    if (cols[23].split(":")[2] == 'fp'):
      sigmap_fp += 1
      sigmap_time_per_chunk.append(mt / chunk)
    if (cols[23].split(":")[2] == 'fn'):
      sigmap_fn += 1
      sigmap_time_per_chunk.append(mt / chunk)
    if (cols[23].split(":")[2] == 'tn'):
      sigmap_tn += 1
      sigmap_time_per_chunk.append(mt / chunk)
  if (len(cols) == 15):
    mt = float(cols[12].split(":")[2])
    if (cols[14].split(":")[2] != 'na'):
      sigmap_time_per_read.append(mt)
    if (cols[14].split(":")[2] == 'fn'):
      sigmap_fn += 1
    if (cols[14].split(":")[2] == 'tn'):
      sigmap_tn += 1
print("Sigmap TP: " + str(sigmap_tp))
print("Sigmap FP: " + str(sigmap_fp))
print("Sigmap FN: " + str(sigmap_fn))
print("Sigmap TN: " + str(sigmap_tn))
sigmap_preciosion = sigmap_tp / (sigmap_tp + sigmap_fp)
print("Sigmap precision: " + str(sigmap_preciosion))
sigmap_recall = sigmap_tp / (sigmap_tp + sigmap_fn)
print("Sigmap recall: " + str(sigmap_recall))
print("Sigmap F-1 score: " + str(2 * sigmap_preciosion * sigmap_recall / (sigmap_preciosion + sigmap_recall)))
print("Mean time per chunk : " + str(mean(sigmap_time_per_chunk)))
print("Median time per chunk : " + str(median(sigmap_time_per_chunk)))
print("Mean time per read : " + str(mean(sigmap_time_per_read)))
print("Median time per read : " + str(median(sigmap_time_per_read)))
print("#Done with sigmap")

for fp in fps:
  fp.close()
