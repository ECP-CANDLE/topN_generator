# topN_generator

```
$ python build.py -h
usage: build.py [-h] [--top_n TOP_N] [--drug_descriptor {dragon7,mordred}]
                [--cell_feature {rnaseq,snps}]
                [--cell_feature_subset {lincs1000,oncogenes,all}]
                [--format {csv,tsv,parquet,hdf5,feather}]
                [--response_type {reg,bin}] [--labels LABELS]

optional arguments:
  -h, --help            show this help message and exit
  --top_n TOP_N         Number of cancer types to be included. Default 6
  --drug_descriptor {dragon7,mordred}
                        Drug descriptors
  --cell_feature {rnaseq,snps}
                        Cell line features
  --cell_feature_subset {lincs1000,oncogenes,all}
                        Subset of cell line features. Default lincs1000
  --format {csv,tsv,parquet,hdf5,feather}
                        Dataframe file format. Default hdf5
  --response_type {reg,bin}
                        Response type. Regression(reg) or Binary
                        Classification(bin). Default reg
  --labels LABELS       Contains Cell and Drug label. Default False
```
