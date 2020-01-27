# topN_generator

```
$ python build.py -h
usage: build.py [-h] [--top_n TOP_N] [--drug_descriptor {dragon7,mordred}]
                [--cell_feature {rnaseq,snps}]
                [--cell_feature_subset {lincs1000,oncogenes,all}]
                [--format {csv,tsv,parquet,hdf5,feather}]
                [--response_type {reg,bin}] [--labels]
                [--target {AUC,IC50,EC50,EC50se,R2fit,Einf,HS,AAC1,AUC1,DSS1}]
                [--scaled]

optional arguments:
  -h, --help            show this help message and exit
  --top_n TOP_N         Number of cancer types to be included. Default 6
  --drug_descriptor {dragon7,mordred}
                        Drug descriptors. Default dragon7
  --cell_feature {rnaseq,snps}
                        Cell line features. Default rnaseq
  --cell_feature_subset {lincs1000,oncogenes,all}
                        Subset of cell line features. Default lincs1000
  --format {csv,tsv,parquet,hdf5,feather}
                        Dataframe file format. Default hdf5
  --response_type {reg,bin}
                        Response type. Regression(reg) or Binary
                        Classification(bin). Default reg
  --labels              Contains Cell and Drug label. Default False
  --target {AUC,IC50,EC50,EC50se,R2fit,Einf,HS,AAC1,AUC1,DSS1}
                        Response label value. Default AUC
  --scaled              Apply scaling. Default False
```
