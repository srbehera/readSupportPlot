import sys
from pysam import VariantFile
import matplotlib.pyplot as plt
import pandas as pd

f = open("plot.csv", "w")
f0 = open("plotlt10.csv", "w")
vcf_in = VariantFile(sys.argv[1])  # auto-detect input format
#vcf_out = VariantFile('-', 'w', header=bcf_in.header)
total = 0
arr=[]
arr = [0 for i in range(30)] 
for rec in vcf_in.fetch():
    info = rec.info
    #info = info.split(";")[-1]
    if rec.info.keys()[-1] == "NIB_READ_SUPPORTS":
        val = int(rec.info.values()[-1][0])
        if val < 10:
            arr[val] = arr[val] + 1
        elif val < 100:
            index = 9 + int(val/10)
            arr[index] = arr[index] + 1
        elif val < 1000:
            index = 18 + int(val/100)
            arr[index] = arr[index] + 1
        else:
            arr[29] = arr[29] + 1
    total = total + 1
print(arr)
for i in range(1, 30):
    if i < 10:
        f.write(str(i)+","+str(arr[i])+"\n")
        f0.write(str(i)+","+str(arr[i])+"\n")
    elif i < 19:
        index = i - 9
        #print(i, index)
        #print(str(index)+"0-"+str(index)+"9,") 
        f.write(str(index)+"0-"+str(index)+"9,"+str(arr[i])+"\n")
    elif i < 28:
        index = i - 18
        f.write("<"+str(index+1)+"00,"+str(arr[i])+"\n")
        #f.write(str(index)+"00-"+str(index)+"99,"+str(arr[i])+"\n")
    else:
        print(arr[i])
        #f.write(">1000,"+arr[i]+"\n")
f.close()
f0.close()
data = pd.read_csv('plot.csv', sep=',',header=None, index_col =0)
data.columns=["called SVs"]
data.plot(kind='line')
plt.ylabel('# SVs')
plt.xlabel('# Read supports')
plt.title('Read supports vs SV counts')
plt.savefig('plot.png')
plt.show()

data = pd.read_csv('plotlt10.csv', sep=',',header=None, index_col =0)
data.columns=["called SVs"]
data.plot(kind='line')
plt.ylabel('# SVs')
plt.xlabel('# Read supports')
plt.title('Read supports vs SV counts')
plt.savefig('plotlt10.png')
plt.show()
