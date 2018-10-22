#!/usr/bin/env python

#Set python environment in the first line as written above

''' This script is written to do the following:
a) Read a .rnt file and save the Location, Strand and Product columns
to another file called filename_LocStrPro.txt
b) Split the Location into two parts and save the Strand and Location
into another file filename_StrLoc1Loc2.txt
c) Based on the strand either increase or decrease the Loc1 and Loc2
numbers and save the changed file to filename_Upstream_StrLoc1Loc2.txt
d) Open the filename.fna file and get the genome presented in both of
the above locations i.e. at original locations and changed locations.
'''

import pandas as pd
import numpy as np
import sys
import collections 

# Read the file. Modify the below to check for the file and ask the
# user to specify the filename after the script. Also need to store
# the filename (without the extension) and use it save other files.

# try:
#   inp_file = sys.argv[1]
# except ValueError:
#   print("No valid integer in line.")

inp_file = sys.argv[1]
#inp_file = 'NC_004061.rnt'

filename=inp_file.split('.')[0]
print(inp_file, filename)

# Read the file as a data frame but Skip the first two rows
df1=pd.read_csv(inp_file,skiprows=2,sep='\t')
print(df1.head())
print(df1.columns)
print(df1.dtypes)

# Store those first two lines (in case they are required)
inp=open(inp_file)
firstline=inp.readline()
secondline=inp.readline()
inp.close()

# Read the gnome file and save it as a string (to slice later)
genomefile=filename+".fna"
with open(genomefile, 'r') as gf:
  trash=gf.readline()
  genome=gf.read().replace('\n', '')
#print(genome)
#print(trash)

# Save only "Location", "Strand", "Product" columns to a separate
# dataframe
df2=df1[["Location", "Strand", "Product"]].copy()

# Split the Location column into startLoc and endLoc. Finally, add
# these columns to the df2 and remove the Location column from the df2
dfTemp = df2["Location"].str.split("\.\.", n = 1, expand=True)
print(dfTemp.head())

df2["startLoc"] = dfTemp[0]
df2["endLoc"] = dfTemp[1]
print(df2.head())
#df2.drop(columns=['Location'], inplace=True)
df2.drop(['Location'], axis=1, inplace=True)

#Change the type of startLoc and endLoc
df2[["startLoc", "endLoc"]] = df2[["startLoc", "endLoc"]].astype(int)

#Change the order of columns
df2 = df2[["Strand", "startLoc", "endLoc", "Product"]]

df_length=len(df2.index)
#Create a numpy array to store product, AT_geneome_percent,
#AT_upstream_percent
product=[]
at_array =np.zeros((df_length,2))

# Print the Genome from start to end location 
genefna='gene_'+filename+'.fna'
gfna=open(genefna,'w')

count=0
for index, row in df2.iterrows():
   this_gene=genome[row['startLoc']:row['endLoc']+1]
   atgc_count = collections.Counter(this_gene)
   at_count = atgc_count['A'] + atgc_count['T']
   gene_length = len(this_gene)
   at_percent = 100*(at_count/gene_length)

   product.append(row['Product'])
   at_array[count,0] = at_percent
   
   gfna.write("%a %a %i %i %a %6.3f\n" % (str('>').strip(''),
row['Strand'], row['startLoc'], row['endLoc'], row['Product'],
at_percent))
   gfna.write("%a \n" % (genome[row['startLoc']:row['endLoc']+1]))
   count = count + 1
   
print(this_gene)
print(at_percent)
print(gene_length)

# Write df2 to separate file with space as the separator.
LSPfile="LSP_gene_"+filename+".txt"
df2.to_csv(LSPfile, sep="\t", index=False)

print(df2.head())
print(df2.columns)
print(df2.dtypes)

# If Strand is -, change startLoc to endLoc and endLoc to endLoc + 50
# If Strand is +, change endLoc to startLoc and startLoc to endLoc - 50
df2.loc[df2.Strand == "-", ['startLoc']] = df2.loc[:,'endLoc'] 
df2.loc[df2.Strand == "-", ['endLoc']] =  df2.loc[:,'endLoc']+50 

df2.loc[df2.Strand == "+", ['endLoc']] =  df2.loc[:,'startLoc']
df2.loc[df2.Strand == "+", ['startLoc']] = df2.loc[:,'startLoc']-50
print(df2.head())

# Write df2 to separate file with space as the separator.
LSPfile="LSP_upstream_"+filename+".txt"
df2.to_csv(LSPfile, sep="\t", index=False)

# Print the upstream Genome from start to end location 
upfna='upstream_'+filename+'.fna'
ufna=open(upfna,'w')

count=0
for index, row in df2.iterrows():
   this_gene=genome[row['startLoc']:row['endLoc']+1]
   atgc_count = collections.Counter(this_gene)
   at_count = atgc_count['A'] + atgc_count['T']
   gene_length = len(this_gene)
   at_percent = 100*(at_count/gene_length)

   at_array[count,1] = at_percent
   
   ufna.write("%a %a %i %i %a %6.3f\n" % (str('>').strip(''),
row['Strand'], row['startLoc'], row['endLoc'], row['Product'],
at_percent))
   ufna.write("%a \n" % (genome[row['startLoc']:row['endLoc']+1]))
   count = count + 1

# Print the both the genome's AT percentage
atPercent='atPercent_'+filename+'.txt'
atp=open(atPercent,'w')
atp.write("%a %a\t %a\n" % ('genAT%', 'upstreamAT%', 'Product'))
for i in range(0,len(product)):
   atp.write("%6.3f\t %6.3f\t %a \n" % (at_array[i,0],
at_array[i,1],str(product[i]).strip("\'")))
