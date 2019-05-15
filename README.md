# topN_generator

```
$ python build.py -h
usage: build.py [-h] [--top_n TOP_N] [--drug_descriptor {dragon7}]
                [--cell_feature {rnaseq}]
                [--format {csv,tsv,parquet,hdf5,feather}]

optional arguments:
  -h, --help            show this help message and exit
  --top_n TOP_N         Number of cancer types to be included. Default 6
  --drug_descriptor {dragon7}
                        Drug descriptors
  --cell_feature {rnaseq}
                        Cell line features
  --format {csv,tsv,parquet,hdf5,feather}
                        Dataframe file format. Default hdf5
  ```
