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
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import collections 
from matplotlib.backends.backend_pdf import PdfPages

# Read the file. Modify the below to check for the file and ask the
# user to specify the filename after the script. Also need to store
# the filename (without the extension) and use it save other files.


inp_file1 = sys.argv[1]
inp_file2 = sys.argv[2]
inp_file3 = sys.argv[3]
inp_file4 = sys.argv[4]

filename1=inp_file1.split('.')[0]
print(inp_file1, filename1)

filename2=inp_file2.split('.')[0]
print(inp_file2, filename2)

filename3=inp_file3.split('.')[0]
print(inp_file3, filename3)

filename4=inp_file4.split('.')[0]
print(inp_file4, filename4)

# Read the file as a data frame but Skip the first two rows
dff1=pd.read_csv(inp_file1,skiprows=2,sep='\t')
dff1['Label'] = filename1 + dff1.Product
# print(dff1.head())
# print(dff1.columns)
# print(dff1.dtypes)

# Read the gnome file and save it as a string (to slice later)
genomefile1=filename1+".fna"
with open(genomefile1, 'r') as gf1:
  trash=gf1.readline()
  genome1=gf1.read().replace('\n', '')
  


# Read the file as a data frame but Skip the first two rows
dff2=pd.read_csv(inp_file2,skiprows=2,sep='\t')
dff2['Label'] =filename2+ dff2.Product
# print(dff2.head())
# print(dff2.columns)
# print(dff2.dtypes)


# Read the gnome file and save it as a string (to slice later)
genomefile2=filename2+".fna"
with open(genomefile2, 'r') as gf2:
  trash=gf2.readline()
  genome2=gf2.read().replace('\n', '')


# Read the file as a data frame but Skip the first two rows
dff3=pd.read_csv(inp_file3,skiprows=2,sep='\t')
dff3['Label'] =filename3+ dff3.Product
# print(dff3.head())
# print(dff3.columns)
# print(dff3.dtypes)


# Read the gnome file and save it as a string (to slice later)
genomefile3=filename3+".fna"
with open(genomefile3, 'r') as gf3:
  trash=gf3.readline()
  genome3=gf3.read().replace('\n', '')




# Read the file as a data frame but Skip the first two rows
dff4=pd.read_csv(inp_file4,skiprows=2,sep='\t')
dff4['Label'] = filename4 + dff4.Product
# print(dff4.head())
# print(dff4.columns)
# print(dff4.dtypes)


# Read the gnome file and save it as a string (to slice later)
genomefile4=filename4+".fna"
with open(genomefile4, 'r') as gf4:
  trash=gf4.readline()
  genome4=gf4.read().replace('\n', '')


dfList1 = dff1['Product'].tolist()
dfList2 = dff2['Product'].tolist()
dfList3 = dff3['Product'].tolist()
dfList4 = dff4['Product'].tolist()
dfList5 = dfList1+dfList2+dfList3+dfList4

uniqproduct= list(set(dfList5))

dflabList1 = dff1['Label'].tolist()
dflabList2 = dff2['Label'].tolist()
dflabList3 = dff3['Label'].tolist()
dflabList4 = dff4['Label'].tolist()
labelist = dflabList1+dflabList2+dflabList3+dflabList4

uniqlabel= list(set(labelist))



print ("uniq list",  uniqproduct, len(uniqproduct), len(uniqlabel))



#Combining all the 4 dataframes and getting one dataframe dfcomb
dftemp= dff2.append(dff1)
dftemp2 = dff3.append(dftemp)
dfcomb = dff4.append(dftemp2)
dfcomb.reset_index(inplace=True, drop=True) 
#print(dfcomb)


filename = "comb"
# Write dfcomb to separate file with space as the separator.
LSPfile="LSP_gene_"+filename+".txt"
dfcomb.to_csv(LSPfile, sep="\t", index=False)


# Save only "Location", "Strand", "Product" columns to a separate
# dataframe
df2=dfcomb[["Location", "Strand", "Product", "Label", "Length"]].copy()

# Split the Location column into startLoc and endLoc. Finally, add
# these columns to the df2 and remove the Location column from the df2
dfTemp = df2["Location"].str.split("\.\.", n = 1, expand=True)
#print(dfTemp.head())

df2["startLoc"] = dfTemp[0]
df2["endLoc"] = dfTemp[1]
df2.drop(['Location'], axis=1, inplace=True)

#Change the type of startLoc and endLoc
df2[["startLoc", "endLoc"]] = df2[["startLoc", "endLoc"]].astype(int)

#Change the order of columns
df2 = df2[["Strand", "startLoc", "endLoc", "Product", "Length", "Label"]]
print(df2)

df_length=len(df2.index)

# Write df2 to separate file with space as the separator.
LSPfile="LSP_gene_"+filename+".txt"
df2.to_csv(LSPfile, sep="\t", index=False)

#print(df2.head(),"gene")
#print(df2.columns)
#print(df2.dtypes)


# If Strand is -, change startLoc to endLoc and endLoc to endLoc + 50
# If Strand is +, change endLoc to startLoc and startLoc to endLoc - 50
df2.loc[df2.Strand == "-", ['startLoc']] = df2.loc[:,'endLoc'] 
df2.loc[df2.Strand == "-", ['endLoc']] =  df2.loc[:,'endLoc']+50 

df2.loc[df2.Strand == "+", ['endLoc']] =  df2.loc[:,'startLoc']
df2.loc[df2.Strand == "+", ['startLoc']] = df2.loc[:,'startLoc']-50


#print(df2.head(),"upstream")

# Write df2 to separate file with space as the separator.
LSPfile="LSP_upstream"+filename+".txt"
df2.to_csv(LSPfile, sep="\t", index=False)

print("condition--------------------------------------")
print(df2.loc[(df2.Strand == "-") & (df2.Length <60)])

# If Strand is -, change startLoc to endLoc and endLoc to endLoc + 50
# If Strand is +, change endLoc to startLoc and startLoc to endLoc - 50
df2.loc[(df2.Strand == "-") & (df2.Length <60), ['startLoc']] = df2.loc[:,'startLoc']-35
df2.loc[(df2.Strand == "-") & (df2.Length <60), ['endLoc']] =  df2.loc[:,'endLoc']

df2.loc[(df2.Strand == "+") & (df2.Length <60), ['endLoc']] =  df2.loc[:,'endLoc']+35
df2.loc[(df2.Strand == "+") & (df2.Length <60), ['startLoc']] = df2.loc[:,'startLoc']

#df2['window']=(df2['endLoc']- df2['startLoc'])-5


# If Strand is -, change startLoc to endLoc and endLoc to endLoc + 50
# If Strand is +, change endLoc to startLoc and startLoc to endLoc - 50
df2.loc[(df2.Strand == "-") & (df2.Length >60), ['startLoc']] = df2.loc[:,'startLoc']-60
df2.loc[(df2.Strand == "-") & (df2.Length >60), ['endLoc']] =  df2.loc[:,'endLoc']

df2.loc[(df2.Strand == "+") & (df2.Length >60), ['endLoc']] =  df2.loc[:,'endLoc']+60
df2.loc[(df2.Strand == "+") & (df2.Length >60), ['startLoc']] = df2.loc[:,'startLoc']

df2['window']=(df2['endLoc']- df2['startLoc'])-5

#print(df2.head(),"upstreamgene")

# Write df2 to separate file with space as the separator.
LSPfile="LSP_upstreamgene"+filename+".txt"
df2.to_csv(LSPfile, sep="\t", index=False)


# Slidingwindow defintion to know the window(ex:(0,7,'ATCG'),(1,8,'TCGTA'))
def slidingWindow(sequence,winSize=7,step=1):

        #Pre-compute number of chunks to emit
         numOfChunks = ((len(sequence)-winSize)/step)+1
         print(numOfChunks)
         print(numOfChunks*step)
         # Do the work
         count = 0
         for i in range(0,int(numOfChunks*step),step):
            count+=1
            yield count,i,i+winSize,sequence[i:i+winSize]
             #print (i, i+winSize, sequence[i:i+winSize])
         #print (count)

df_max = 105
print (df_max ,"dfmaxvalue---------------")
#Create a numpy array to store productin rows and windows in columns.
product=[]
at_array =np.zeros((df_length,df_max))
#print(at_array)

#Print the gen+upstream and At%  from start to end location into fasta file
geneupfna='upstreamgene'+filename+'.fna'
geneupfna=open(geneupfna,'w')


count=0
for index, row in df2.iterrows():
 print(row.Label)
 #print(filename4)
 #print(df2.Label)
 
 if filename4 in row.Label:
        print("in file4")
        this_gene=genome4[row['startLoc']:row['endLoc']+1]
        print(this_gene)
        chunks = slidingWindow(this_gene,7,1)
        #product.append(row['Product'])
        for chunk in chunks :
          print(chunk[0],chunk[1],chunk[2],chunk[3])
          atgc_count = collections.Counter(chunk[3])
          at_count = atgc_count['A'] + atgc_count['T']
          gene_length = len(chunk[3])
          #gene_length = len(this_gene)
          #product.append(row['Product'])
          at_percent = 100*(at_count/gene_length)
     
        geneupfna.write("%a %a %i %i %a %6.3f\n" % (str('>').strip(''),
row['Strand'], row['startLoc'], row['endLoc'], row['Label'],
at_percent))
        geneupfna.write("%a \n" % (genome4[row['startLoc']:row['endLoc']+1]))
 if filename1 in row.Label:
        print("in file1")
        this_gene=genome1[row['startLoc']:row['endLoc']+1]
        print(this_gene)
        chunks = slidingWindow(this_gene,7,1)
        #product.append(row['Product'])
        for chunk in chunks :
          # print(chunk[0],chunk[1],chunk[2],chunk[3])
          atgc_count = collections.Counter(chunk[3])
          at_count = atgc_count['A'] + atgc_count['T']
          gene_length = len(chunk[3])
          #product.append(row['Product'])
          at_percent = 100*(at_count/gene_length)
     
        geneupfna.write("%a %a %i %i %a %6.3f\n" % (str('>').strip(''),
row['Strand'], row['startLoc'], row['endLoc'], row['Label'],
at_percent))
        geneupfna.write("%a \n" % (genome1[row['startLoc']:row['endLoc']+1]))
 if filename3 in row.Label:
      print("in file3")
      this_gene=genome3[row['startLoc']:row['endLoc']+1]
      chunks = slidingWindow(this_gene,7,1)
      #product.append(row['Product'])
      for chunk in chunks :
         #print(chunk[0],chunk[1],chunk[2],chunk[3])
         atgc_count = collections.Counter(chunk[3])
         at_count = atgc_count['A'] + atgc_count['T']
         gene_length = len(chunk[3])
         #product.append(row['Product'])
         at_percent = 100*(at_count/gene_length)
     
      geneupfna.write("%a %a %i %i %a %6.3f\n" % (str('>').strip(''),
row['Strand'], row['startLoc'], row['endLoc'], row['Label'],
at_percent))
      geneupfna.write("%a \n" % (genome3[row['startLoc']:row['endLoc']+1]))
 if  filename2 in row.Label:
      print("in file2")
      this_gene=genome2[row['startLoc']:row['endLoc']+1]
      chunks = slidingWindow(this_gene,7,1)
      #product.append(row['Product'])
      for chunk in chunks :
         # print(chunk[0],chunk[1],chunk[2],chunk[3])
         atgc_count = collections.Counter(chunk[3])
         at_count = atgc_count['A'] + atgc_count['T']
         gene_length = len(chunk[3])
         #product.append(row['Product'])
         at_percent = 100*(at_count/gene_length)
     
      geneupfna.write("%a %a %i %i %a %6.3f\n" % (str('>').strip(''),
row['Strand'], row['startLoc'], row['endLoc'], row['Label'],
at_percent))
      geneupfna.write("%a \n" % (genome2[row['startLoc']:row['endLoc']+1]))
count = count + 1

## Print the genenames along with windows and  AT percentage in each window
window='window'+filename+'.txt'
win=open(window,'w')
win.write("%a\t %a\t %a\t %a\t %a\t %a\n" % ('wcount','Wstart','Wend','wseq','AT%', 'Label'))
for index, row in df2.iterrows():
 product.append(row['Label'])
 print(row.Label)
 print(filename4)
 print(df2.Label)
 
 if filename4 in row.Label:
       print("in file4")
       #product.append(row['Label'])
       # print (index, "indexxx+====")
       this_gene=genome4[row['startLoc']:row['endLoc']+1]
       chunks = slidingWindow(this_gene,7,1)
       count =0
       for chunk in chunks :
     # count =0
     # print(chunk[0],chunk[1],chunk[2],chunk[3])
        atgc_count = collections.Counter(chunk[3])
        at_count = atgc_count['A'] + atgc_count['T']
        gene_length = len(chunk[3])
        at_percent = 100*(at_count/gene_length)
        at_array[index,count] = at_percent
     #win.write("%a\t %a\t %a\t %a\n" % ('Wstart','Wend','Product','AT%'))
        win.write("%i\t %i\t %i\t %a\t %6.3f\t  %a\n" % (chunk[0],chunk[1],chunk[2],chunk[3],at_percent,row['Label']))
     #geneupfna.write("%a \n" % (genome[row['startLoc']:row['endLoc']+1]))
        count = count + 1

 if filename1 in row.Label:
       print("in file1")
       #product.append(row['Label'])
       # print (index, "indexxx+====")
       this_gene=genome1[row['startLoc']:row['endLoc']+1]
       chunks = slidingWindow(this_gene,7,1)
       count =0
       for chunk in chunks :
     # count =0
     # print(chunk[0],chunk[1],chunk[2],chunk[3])
        atgc_count = collections.Counter(chunk[3])
        at_count = atgc_count['A'] + atgc_count['T']
        gene_length = len(chunk[3])
        at_percent = 100*(at_count/gene_length)
        at_array[index,count] = at_percent
     #win.write("%a\t %a\t %a\t %a\n" % ('Wstart','Wend','Product','AT%'))
        win.write("%i\t %i\t %i\t %a\t %6.3f\t  %a\n" % (chunk[0],chunk[1],chunk[2],chunk[3],at_percent,row['Label']))
     #geneupfna.write("%a \n" % (genome[row['startLoc']:row['endLoc']+1]))
        count = count + 1

 if filename3 in row.Label:
      print("in file3")
      #product.append(row['Label'])
      # print (index, "indexxx+====")
      this_gene=genome3[row['startLoc']:row['endLoc']+1]
      chunks = slidingWindow(this_gene,7,1)
      count =0
      for chunk in chunks :
     # count =0
     # print(chunk[0],chunk[1],chunk[2],chunk[3])
        atgc_count = collections.Counter(chunk[3])
        at_count = atgc_count['A'] + atgc_count['T']
        gene_length = len(chunk[3])
     #product.append(row['Product'])
        at_percent = 100*(at_count/gene_length)
        at_array[index,count] = at_percent
     #win.write("%a\t %a\t %a\t %a\n" % ('Wstart','Wend','Product','AT%'))
        win.write("%i\t %i\t %i\t %a\t %6.3f\t  %a\n" % (chunk[0],chunk[1],chunk[2],chunk[3],at_percent,row['Label']))
     #geneupfna.write("%a \n" % (genome[row['startLoc']:row['endLoc']+1]))
        count = count + 1

 if  filename2 in row.Label:
      print("in file2")
      #product.append(row['Label'])
      # print (index, "indexxx+====")
      this_gene=genome2[row['startLoc']:row['endLoc']+1]
      chunks = slidingWindow(this_gene,7,1)
      count =0
      for chunk in chunks :
     # count =0
     # print(chunk[0],chunk[1],chunk[2],chunk[3])
        atgc_count = collections.Counter(chunk[3])
        at_count = atgc_count['A'] + atgc_count['T']
        gene_length = len(chunk[3])
     #product.append(row['Product'])
        at_percent = 100*(at_count/gene_length)
        at_array[index,count] = at_percent
     #win.write("%a\t %a\t %a\t %a\n" % ('Wstart','Wend','Product','AT%'))
        win.write("%i\t %i\t %i\t %a\t %6.3f\t  %a\n" % (chunk[0],chunk[1],chunk[2],chunk[3],at_percent,row['Label']))
     #geneupfna.write("%a \n" % (genome[row['startLoc']:row['endLoc']+1]))
        count = count + 1




# Read the file as a data frame but Skip the first two rows
dfwindow=pd.read_csv("windowcomb.txt",sep='\t')
print(dfwindow.head())
print(dfwindow.columns)
print(dfwindow.dtypes)


# Save only "wcount", "AT%", "Label" columns to a separate
# dataframe
dfdrop=dfwindow.copy()
dffinal=dfdrop.drop(dfdrop.columns[[1, 2, 3]], axis=1) 
print(dffinal.head())
print(dffinal.columns)
print(dffinal.dtypes)


#writing the SlidingwindowAnalaysis to excel  
print(at_array)
print(np.shape(at_array))

df5= pd.DataFrame(at_array)
print(len(df5.columns),"---------df4columns")
df5.columns = pd.RangeIndex(1, len(df5.columns)+1) 
df5['Label']=df2['Label']
print (df5)


excel= 'final_'+filename+'.xlsx'
writer = pd.ExcelWriter(excel)
df5.to_excel(writer,'Sheet1')
#df10.to_excel(writer,'Sheet6')
writer.save()


def colour(df):
    for i in list(df.columns): 
       if filename1 in  i:
          return'darkgreen'
       if filename2 in  i:
          return 'navy'
       if filename3 in  i:
          return 'black'
       if filename4 in  i:
          return 'darkorange'
          
def color(df):
    my_list = []
    for i in list(df.columns): 
       if i.startswith(filename1):
          return my_list.append('darkgreen')
       if i.startswith(filename2):
          return my_list.append('navy')
       if i.startswith(filename3):
          return my_list.append('black')
       if i.startswith(filename4):
          return my_list.append('darkorange')
    return my_list



def mycolour(df):
    my_list = []
    for i in list(df.columns): 
       if filename1 in  i:
          return my_list.append('darkgreen')
       if filename2 in  i:
          return my_list.append('navy')
       if filename3 in  i:
          return my_list.append('black')
       if filename4 in  i:
          return my_list.append('darkorange')
    print("mylist-----------------"+ my_list)
    return my_list  
    

with PdfPages('multipage_pdf.pdf') as pdf:
    for prod in uniqproduct:
      dfname = "df"+prod
      print(dfname)
      dfname =df5[df5['Label'].str.contains(prod)]
      dfname.reset_index(inplace=True, drop=True)
      #dfdrop = dfname.drop('Label', 1)
      print(dfname)
      dftran= dfname.T
      dftran.columns =dftran.iloc[105]
      #dftran.reindex(dftran.index.drop(80))
      dftran.drop(dftran.tail(1).index,inplace=True)
      print("dftran index++++++++++++")
      print(dftran.index)
      #dftran['windows']=dftran.index
      #dftran=dftran.drop(dftran.columns[[0]], axis=1)
      print(dftran)
      #plt.plot(dftran.windows,dftran[dftran.columns[0]])
      plt.figure(figsize=(10,10))
      my_list = []
      #plt.plot(dftran.index,dftran[[i for i in list(dftran.columns)]],c=color(dftran))
      #plt.plot(dftran.index,dftran[[i for i in list(dftran.columns)]])
      for i in list(dftran.columns): 
           if "NC_004061" in i:
             plt.plot(dftran.index,dftran[i],color='dodgerblue')
             my_list.append(i)
           if "NC_011833" in i:
              plt.plot(dftran.index,dftran[i],color='seagreen')
              my_list.append(i)
           if "NC_017256" in i:
              plt.plot(dftran.index,dftran[i],color='firebrick')
              my_list.append(i)
           if "NC_017259" in i:
              plt.plot(dftran.index,dftran[i],color='black')
              my_list.append(i)
           
      plt.xlim([0, 120])
      plt.ylim([0, 100])
      #print(i)  
      plt.title(prod)
      plt.xlabel("Windows")
      plt.ylabel("AT%")
      #plt.legend([i for i in list(dftran.columns) if i != 'windows'],loc='lower right',ncol=2,borderaxespad=1, )
      plt.legend([i for i in list(dftran.columns) ],loc='upper center', bbox_to_anchor=(0.5, -0.05),
fancybox=True, shadow=True, ncol=3)
      print(i for i in list(dftran.columns) if i != 'windows')
      pdf.savefig()
      plt.close()
      
      













