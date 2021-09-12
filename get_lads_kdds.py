# train_ctspecific_HMMs.py

from pomegranate import *
import pandas as pd
import json
import os
import sys
from lad_utils import *

# paths for LADs
bedgraphs_dir = 'bedgraphs_lb1'
outdir = 'LADs_out'
out_beds = 'LAD_BEDs

## paths for KDDs
#bedgraphs_dir = 'bedgraphs_h3k9me2'
#outdir = 'KDDs_out'
#out_beds = 'KDD_BEDs

org_df = pd.read_excel('organizational_df.xlsx')

# states dict for LADs
states_dict = {
        2:['nonLAD','LAD'],
        3:['nonLAD','T2-LAD','T1-LAD'],
        4:['nonLAD','T3-LAD','T2-LAD','T1-LAD'],
        5:['nonLAD','T4-LAD','T3-LAD','T2-LAD','T1-LAD']
    }
    
# states dict for KDDs
states_dict = {
        2:['nonKDD','KDD'],
        3:['nonKDD','T2-KDD','T1-KDD'],
        4:['nonKDD','T3-KDD','T2-KDD','T1-KDD'],
        5:['nonKDD','T4-KDD','T3-KDD','T2-KDD','T1-KDD']
    }

def parse_bedgraph(f, binsize, blacklistbed):
    # parses bedgraph file for use in the HMM training

    canonical_chroms = [ 'chr' + str(x) for x in range(1,23) ] + ['chrY', 'chrX']
    
    indf_pre = pd.read_csv(f, sep='\t', low_memory=False, header=None, 
        names=['chrom','start','stop','score']).dropna().query('chrom in @canonical_chroms')
    indf_pre['id'] = indf_pre.index
    
    inbed = BedTool.from_dataframe(indf_pre[['chrom','start','stop','id','score']]).sort()
    inbed_filt = inbed.subtract(blacklistbed).sort()
    
    indf = inbed_filt.to_dataframe()
    indf.columns = ['chrom','start','stop','id','score']
    indf = indf[['chrom','start','stop','score']].copy()

    indf['size'] = indf['stop'] - indf['start']

    indf_split = indf.query('size > @binsize').copy()
    if len(indf_split) > 0:
        to_add = []
        for ix, row in indf_split.iterrows():
            n_intervals = row['size'] / binsize
            start = row['start']
            chrom = row['chrom']
            score = row[f'score']
            for val in range(int(n_intervals)):
                to_add.append(pd.DataFrame({
                    'chrom':[row['chrom']],
                    'start':[start],
                    'stop':[start + binsize],
                    'score':[score]
                }))
                start = start + binsize

        indf = pd.concat([indf.query('size == @binsize')[['chrom','start','stop',f'score']],
                          pd.concat(to_add)]).query('chrom in @canonical_chroms')
    else:
        indf = indf.query('size == @binsize')[['chrom','start','stop',f'score']].query('chrom in @canonical_chroms')
    return(indf)

for ix, row in org_df.query('chip == "LB1"').iterrows():
    for n_states in [2, 3, 4, 5]:
        ct = row['cell_type']
        dat_rep1 = parse_bedgraph(f'{bedgraphs_dir}/{ct}_rep1.sorted.bedgraph',20000, encode_blacklist)
        dat_rep2 = parse_bedgraph(f'{bedgraphs_dir}/{ct}_rep2.sorted.bedgraph',20000, encode_blacklist)
        model = HiddenMarkovModel()
        dat_model = model.from_samples(NormalDistribution, n_components=n_states, X=[dat_rep1['score'].tolist(),dat_rep2['score'].tolist()], algorithm='baum-welch', verbose=True, n_jobs=20)
        dat = dat_rep1.merge(dat_rep2, on=['chrom','start','stop'], suffixes=('_rep1', '_rep2'))
        dat['hmm_pred'] = dat_model.predict(dat[['score_rep1','score_rep2']].to_numpy(),
                                        algorithm='viterbi')[1:]
        
        trained_model_json = dat_model.to_json()
        if not os.path.exists(f'{outdir}/{n_states}states'):
            os.makedirs(f'{outdir}/hmms_{n_states}states')
            os.makedirs(f'{outdir}/hmm_calls_{n_states}states')
            os.makedirs(f'{outdir}/BED_files_{n_states}states')
        with open(f'{outdir}/hmms_{n_states}states/{ct}.json', 'w') as outfile:
            json.dump(trained_model_json, outfile)
        dat.to_csv(f'{outdir}/hmm_calls_{n_states}states/{ct}.tsv',
              sep='\t', index=False)
        
        dat = pd.read_table(f'{outdir}/hmm_calls_{n_states}states/{ct}.tsv', names=['chrom','start','stop','score0','score1','hmm_pred'],
                           header=0)
        dat_w_cats = assign_categories(dat, states_dict[n_states])
        for cat in dat_w_cats['category'].unique():
            bed = BedTool.from_dataframe(dat_w_cats.query('category == @cat')).sort().merge(d=1)
            bed.to_dataframe().to_csv(f'{out_beds}/BED_files_{n_states}states/{ct}_{cat}s.bed',
                             sep='\t', header=False, index=False)
