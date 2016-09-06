"""
helper script to combine a set of insert site tables into one big table
the files should all have extensions of .table and should be in a
heirarchy where the directory name they are in indicates the date from the
run
"""
import pandas as pd
import os
import fnmatch

def read_dataframe(fn):
    date = os.path.basename(os.path.dirname(fn))
    sample = os.path.basename(fn).split("_")[0]
    df = pd.DataFrame.from_csv(fn, sep="\t")
    if df.shape[0] == 0:
        df.loc[1] = "NA"
    df["sample"] = sample
    df["date"] = date
    return df

def read_and_combine_dataframes(fns):
    df = read_dataframe(fns[0])
    for fn in fns[1:]:
        new_df = read_dataframe(fn)
        df = df.append(new_df)
    return df

def find_tables(base_dir):
    matches = []
    for root, dirnames, filenames in os.walk(base_dir):
      for filename in fnmatch.filter(filenames, '*.deduped.table'):
          matches.append(os.path.join(root, filename))
    return matches

tables = find_tables("/Users/rory/cache/li_hiv/align")
df = read_and_combine_dataframes(tables)
df.to_csv("/Users/rory/cache/li_hiv/align/chimeric.tables", sep="\t")
