import pandas as pd
import sys

args = sys.argv

file1 = args[1]
outfile = args[2]

df1 = pd.read_csv(file1, delimiter="\t")
#df2= pd.read_csv(file2, delimiter="\t")

concat= pd.concat([df1],axis=0)

#concat.drop(columns=['ID'], inplace=True)
concat.replace(to_replace='.', value='-1', inplace=True)
concat.to_csv(outfile, index=False)
