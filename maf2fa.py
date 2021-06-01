##stitch together the individual alignment fragments from a maf produced by hal2maf. Care is taken to infer and fill in gaps not reported in the maf.
#python 3
#use as: python maf2fa.py $NUM
import pandas as pd
import sys
NUM = int(sys.argv[1])
df = pd.DataFrame()
for x in range(0, NUM):
 if df.empty:
  name = "maf_part"+str(x+1)+".txt"
  a = pd.read_csv(name,sep="\t",header=None,index_col=1)
  df = a[7].loc[~a.index.duplicated(keep='first')] # take fasta string and if split over 2 chr, take first instance
  df = df[pd.to_numeric(df, errors='coerce').isnull()] #if not fasta, remove
  df = pd.concat([df], axis=1, sort=False)
 elif not df.empty:
  name = "maf_part"+str(x+1)+".txt"
  a = pd.read_csv(name,sep="\t",header=None,index_col=1)
  b = a[7].loc[~a.index.duplicated(keep='first')]
  b = b[pd.to_numeric(b, errors='coerce').isnull()]
  df = pd.concat([df,b], axis=1, sort=False)

#make col names
namer = []
for x in range(0, NUM):
 inner = "a"+str(x+1)
 namer.append(inner)

df.columns = namer

##fill missing fasta seqments with the correct amount of missing spaces
for x in namer:
 filler = "-"*int(df[x].str.len().mean())
 df[x].fillna(filler, inplace=True)

##concatenate and output
df = df.apply(lambda x: ''.join(x), axis = 1)
df.to_csv("prep.txt",sep="\t",header=False, index=True)

