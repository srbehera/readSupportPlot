import sys
from pysam import VariantFile
import matplotlib.pyplot as plt
import pandas as pd

f = open("plot.csv", "w")
f0 = open("plotlt10.csv", "w")
#f.write(" ,FP,TP\n")
#f0.write(" ,FP,TP\n")
vcf_in = VariantFile(sys.argv[1])  # auto-detect input format
#vcf_out = VariantFile('-', 'w', header=bcf_in.header)
total = 0
arr_TP=[]
arr_FP=[]
arr_TP = [0 for i in range(30)] 
arr_FP = [0 for i in range(30)]
for rec in vcf_in.fetch():
    record = str(rec)
    info = rec.info
    if rec.info.keys()[-1] == "NIB_READ_SUPPORTS":
        if record.find("HG2count=0;") != -1:
            val = int(rec.info.values()[-1][0])
            if val < 10:
                arr_FP[val] = arr_FP[val] + 1
            elif val < 100:
                index = 9 + int(val/10)
                arr_FP[index] = arr_FP[index] + 1
            elif val < 1000:
                index = 18 + int(val/100)
                arr_FP[index] = arr_FP[index] + 1
            else:
                arr_FP[29] = arr_FP[29] + 1
        else:
            #print(rec)
            #exit()
            val = int(rec.info.values()[-1][0])
            if val < 10:
                arr_TP[val] = arr_TP[val] + 1
            elif val < 100:
                index = 9 + int(val/10)
                arr_TP[index] = arr_TP[index] + 1
            elif val < 1000:
                index = 18 + int(val/100)
                arr_TP[index] = arr_TP[index] + 1
            else:
                arr_TP[29] = arr_TP[29] + 1
#print(arr)
for i in range(1, 30):
    if i < 10:
        f.write(str(i)+","+str(arr_FP[i])+","+str(arr_TP[i])+"\n")
        f0.write(str(i)+","+str(arr_FP[i])+","+str(arr_TP[i])+"\n")
    elif i < 19:
        index = i - 9
        #print(i, index)
        #print(str(index)+"0-"+str(index)+"9,") 
        f.write(str(index)+"0-"+str(index)+"9,"+str(arr_FP[i])+","+str(arr_TP[i])+"\n")
    elif i < 28:
        index = i - 18
        f.write("<"+str(index+1)+"00,"+str(arr_FP[i])+","+str(arr_TP[i])+"\n")
        #f.write(str(index)+"00-"+str(index)+"99,"+str(arr[i])+"\n")
    else:
        print(arr_FP[i], arr_TP[i])
        #f.write(">1000,"+arr[i]+"\n")
f.close()
f0.close()
data = pd.read_csv('plot.csv', sep=',',header=None, index_col =0)
data.columns=['FP', 'TP']
data.plot(kind='line')
plt.ylabel('# SVs')
plt.xlabel('# Read supports')
plt.title('Read supports vs SV counts')
plt.savefig('plot.png')
plt.show()

data = pd.read_csv('plotlt10.csv', sep=',',header=None, index_col =0)
data.columns=['FP', 'TP']
data.plot(kind='line')
plt.ylabel('# SVs')
plt.xlabel('# Read supports')
plt.title('Read supports vs SV counts')
plt.savefig('plotlt10.png')
plt.show()

data = pd.read_csv('plot.csv', sep=',',header=None, index_col =0)
data.columns=['FP', 'TP']
data.plot(kind='bar')
plt.ylabel('# SVs')
plt.xlabel('# Read supports')
plt.title('Read supports vs SV counts')
plt.savefig('plot_bar.png')
plt.show()

data = pd.read_csv('plotlt10.csv', sep=',',header=None, index_col =0)
data.columns=['FP', 'TP']
data.plot(kind='bar')
plt.ylabel('# SVs')
plt.xlabel('# Read supports')
plt.title('Read supports vs SV counts')
plt.savefig('plotlt10_bar.png')
plt.show()
