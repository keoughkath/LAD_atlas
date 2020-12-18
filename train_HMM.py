# !/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Train an HMM using bedGraph files as input
Kathleen Keough 2020 Pollard Lab, UCSF

Usage: 
    train_HMM.py <train_dir> <n_states> <n_jobs> <out> [--modelname=<model_name>]
    train_HMM.py -h

Arguments:
  train_dir             Directory with bedgraphs containing the data to train the HMM
  n_states              # of states desired in the HMM
  n_jobs                # of jobs (processes) you want to use to speed up the computation
  out                   Where to saved the trained model (should end in .json)
  model_name            What to call the model [default: HMM]

"""

__version__ = 0.0

from pomegranate import *
import pandas as pd
import json
from functools import reduce
from docopt import docopt
import os

def main(args):

  print(args)

  canonical_chroms = [ 'chr' + str(x) for x in range(1,23) ] + ['chrY', 'chrX']

  # separate sequences for model training

  counter = -1
  dfs = []
  cols = []
  for f in os.listdir(args['<train_dir>']):
    if f.endswith('.bedgraph'):
      counter += 1
      indf = pd.read_csv(os.path.join(args['<train_dir>'], f),
                         sep='\t', header=None, names = ['chrom',
                                                         'start',
                                                         'stop',
                                                         f'score{counter}']).query('chrom in @canonical_chroms')
      cols.append(f'score{counter}')
      dfs.append(indf)
              
  # merge these, replacing places where value may not exist in one + dfs with 0.0
  merge_df_out = reduce(lambda  left,right: pd.merge(left,right,on=['chrom','start','stop'],
                                      how='outer'), dfs).fillna(0.0)

  # set up initial parameters for model, initialize and train

  model = HiddenMarkovModel(name=args['--modelname'])
  model = model.from_samples(NormalDistribution, n_components=int(args['<n_states>']), X=merge_df_out[cols].to_numpy(),
                                 algorithm='baum-welch', verbose=True, n_jobs=int(args['<n_jobs>']))


  # save the trained model as a json

  trained_model_json = model.to_json()

  with open(args['<out>'], 'w') as outfile:
      json.dump(trained_model_json, outfile)


if __name__ == "__main__":
    arguments = docopt(__doc__, version=__version__)
    main(arguments)
