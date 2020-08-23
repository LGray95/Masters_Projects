import pandas as pd

df1 = pd.read_csv("/Users/lachlan/Desktop/comparedata.txt", sep="\t")
df2 = pd.read_csv("/Users/lachlan/Desktop/circRNA_list.txt")

outdf = pd.merge(df1, df2, how='inner', on='circRNA_ID')

outdf.to_csv("/Users/lachlan/Desktop/DE_circRNA_merge.txt", sep="\t")