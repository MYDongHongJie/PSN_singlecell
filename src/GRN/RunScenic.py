#!/usr/bin/env python
# pyscenic command line wrapper 

import os
import sys
import numpy as np 
import pandas as pd
import glob
import argparse
import pprint
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names
from pyscenic.utils import modules_from_adjacencies
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.prune import  prune2df, df2regulons
from pyscenic.aucell import aucell



if __name__ == "__main__":

    # parse the command line parameters and store them into a dict
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,
                        help="[REQUIRED]the tab seperated count matrix of the scRNA-seq data.")
    parser.add_argument("-j", "--threads", type=int, default = 10,
                        help="the threads used to run the pyscenic program.")
    parser.add_argument("-r", "--rankdb", type=str,
                        help="[REQUIRED]the directory of databases ranking the whole genome of your \
                            species of interest based on regulatory features (i.e. transcription factors).\
                            Ranking databases are typically in the feather format(*.feather) and can \
                            be downloaded from cisTargetDBs.")
    parser.add_argument("-m", "--motif", type=str,
                        help= "[REQUIRED]the table of motif annotation database providing the link between \
                             an enriched motif and the transcription factor that binds this motif.")
    parser.add_argument( "-s", "--species", type=str,
                         help="the organism prefix for the rankdb.For example: \
                              the prefix for human rankdb hg19-*.feather is hg19,\
                              mm9 for mm9-*.feather.")
    parser.add_argument("-l", "--TF", type=str,
                        help="[REUIRED]the transcript factor list for the species from other source for this run.")
    parser.add_argument("-o", "--outdir", type=str, help="the output directory of this run.")
    args = parser.parse_args()

    if ( not os.path.exists(args.outdir) ):
        os.mkdir(args.outdir) 

    # loading ranking database
    db_fnames = glob.glob(os.path.join(args.rankdb, args.species+"*.feather"))
    dbs = [RankingDatabase(fname=fname , name=os.path.splitext(os.path.basename(fname))[0]) for fname in db_fnames]
    # read in the transcrition factor list
    tf_names = load_tf_names(args.TF)

    # read in the count matrix
    ex_matrix = pd.read_csv(args.input, sep='\t', header=0, index_col=0).T

    # construct the adjacency matrix using the count matrix
    adjacencies = grnboost2(expression_data=ex_matrix, tf_names=tf_names, verbose=True)

    adjacencies.to_csv(os.path.join(args.outdir,"adjacencies.tsv"), index=False, sep='\t')
    # find co-expression gene modules
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
    df = prune2df(dbs, modules, args.motif )
    df.to_csv(os.path.join(args.outdir, "motifs.tsv"),index=False, sep='\t')
    regulons = df2regulons(df)
    auc_mtx = aucell(ex_matrix, regulons, num_workers= args.threads)
    auc_mtx.to_csv(os.path.join(args.outdir,"auc.tsv"), index=False, sep='\t')

    export2loom(ex_mtx=ex_matrix, auc_mtx=auc_mtx,
                regulons=[r.rename(r.name.replace('(+)', ' (' + str(len(r)) + 'g)')) for r in regulons],
                out_fname=os.path.join(args.outdir, "regulons.loom"))

