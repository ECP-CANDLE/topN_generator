#!/usr/bin/env python

import pandas as pd
import argparse
from pathlib import Path

# input files
response_path = Path('./data/combined_single_response_agg')
cell_cancer_types_map_path = Path('./data/combined_cancer_types')
drug_list_path = Path('./data/drugs_1800')
drug_descriptors_path = Path('./data/combined_dragon7_descriptors')
cell_rnaseq_path = Path('./data/combined_rnaseq_data_lincs1000_combat')


def parse_arguments(model_name=''):
    parser = argparse.ArgumentParser()
    parser.add_argument('--top_n', type=int, default=6,
                        help='Number of cancer types to be included')
    parser.add_argument('--drug_descriptor', type=str, default='dragon7',
                        choices=['dragon7'],
                        help='Drug descriptors')
    parser.add_argument('--cell_feature', default='rnaseq',
                        choices=['rnaseq'],
                        help='Cell line features')
    parser.add_argument('--format', default='hdf5',
                        choices=['csv', 'tsv', 'parquet', 'hdf5', 'feather'],
                        help='Dataframe file format')

    args, unparsed = parser.parse_known_args()
    return args, unparsed


def build_dataframe(args):

    # Identify Top N cancer types
    df_response = pd.read_csv(response_path, sep='\t', engine='c', low_memory=False)

    df_uniq_cl_drugs = df_response[['CELL', 'DRUG']].drop_duplicates().reset_index(drop=True)

    df_cl_cancer_map = pd.read_csv(cell_cancer_types_map_path, sep='\t', header=None, names=['CELL', 'CANCER_TYPE'])
    df_cl_cancer_map.set_index('CELL')

    df_cl_cancer_drug = df_cl_cancer_map.merge(df_uniq_cl_drugs, on='CELL', how='left', sort='true')
    df_cl_cancer_drug['CELL_DRUG'] = df_cl_cancer_drug.CELL.astype(str) + '.' + df_cl_cancer_drug.DRUG.astype(str)

    top_n = df_cl_cancer_drug.groupby(['CANCER_TYPE']).count().sort_values('CELL_DRUG', ascending=False).head(args.top_n)
    top_n_cancer_types = top_n.index.to_list()

    print("Identified {} cancer types: ".format(args.top_n), top_n_cancer_types)

    # Indentify cell lines associated with the target cancer types
    df_cl = df_cl_cancer_drug[df_cl_cancer_drug['CANCER_TYPE'].isin(top_n_cancer_types)][['CELL']].drop_duplicates().reset_index(drop=True)

    # Identify drugs associated with the target cancer type & filtered by drug_list

    df_drugs = df_cl_cancer_drug[df_cl_cancer_drug['CANCER_TYPE'].isin(top_n_cancer_types)][['DRUG']].drop_duplicates().reset_index(drop=True)

    drug_list = pd.read_csv(drug_list_path)['DRUG'].to_list()
    df_drugs = df_drugs[df_drugs['DRUG'].isin(drug_list)].reset_index(drop=True)

    # Filter response by cell lines (4882) and drugs (1779)
    cl_filter = df_cl.CELL.to_list()
    dr_filter = df_drugs.DRUG.to_list()

    df_response = df_response[df_response.CELL.isin(cl_filter) & df_response.DRUG.isin(dr_filter)][['CELL', 'DRUG', 'AUC']].drop_duplicates().reset_index(drop=True)

    # Join response data with Drug descriptor & RNASeq
    df_rnaseq = pd.read_csv(cell_rnaseq_path, sep='\t', low_memory=False)
    df_rnaseq = df_rnaseq[df_rnaseq['Sample'].isin(cl_filter)].reset_index(drop=True)

    df_rnaseq.rename(columns={'Sample': 'CELL'}, inplace=True)
    df_rnaseq = df_rnaseq.set_index(['CELL'])

    df_descriptor = pd.read_csv(drug_descriptors_path, sep='\t', low_memory=False, na_values='na')
    df_descriptor = df_descriptor[df_descriptor.DRUG.isin(dr_filter)].set_index(['DRUG']).fillna(0)

    df = df_response.merge(df_rnaseq, on='CELL', how='left', sort='true')
    df.set_index(['DRUG'])

    df_final = df.merge(df_descriptor, on='DRUG', how='left', sort='true')
    df_final.drop(columns=['CELL', 'DRUG'], inplace=True)
    print("Dataframe is built with total {} rows.".format(len(df_final)))

    save_filename = 'top_{}.{}'.format(args.top_n, args.format)
    print("Saving to {}".format(save_filename))

    if args.format == 'feather':
        df_final.to_feather(save_filename)
    elif args.format == 'csv':
        df_final.to_csv(save_filename, index=False)
    elif args.format == 'tsv':
        df_final.to_csv(save_filename, sep='\t', index=False)
    elif args.format == 'parquet':
        df_final.to_parquet(save_filename, index=False)
    elif args.format == 'hdf5':
        df_final.to_hdf(save_filename, key='df', mode='w', complib='blosc:snappy', complevel=9)


if __name__ == '__main__':
    FLAGS, unparsed = parse_arguments()
    build_dataframe(FLAGS)
