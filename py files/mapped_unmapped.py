import pandas as pd
import csv
import math

mapped=[]
unmapped=[]
def calculate_sym(item1,item2):
	val=0
	for i in range(0,len(item1)):
		val=((item1[i]-item2[i])**2)+val
	return round(math.sqrt(val),7)

def find_mapped_unmapped_index():
	list_=[]
	input_file='yeast_genes.txt'
	dataset = pd.read_csv(input_file,sep=" ",header= None)
	csv_list = dataset.iloc[0:,0:].values
	print csv_list.shape
	input_file='unmapped.txt'
	dataset1 = pd.read_csv(input_file,sep=" ",header = None)
	csv_list1 = dataset1.iloc[0,0:].values
	print csv_list1.shape
	for i in range(0,237):
		list_.append(csv_list1[i])
	for i in range(0,2884):
		if (csv_list[i,0] in list_):
			unmapped.append(i)
		else:
			mapped.append(i)
find_mapped_unmapped_index()
lencol=0
lenrow=0
input_file='microarray_2884.txt'
dataset = pd.read_csv(input_file,sep=" ",header=None)
csv_list = dataset.iloc[0:,0:].values
output_file='microarray_2647.csv'
with open(output_file,'w') as outputcsv_file:
	spamwriter = csv.writer(outputcsv_file,delimiter=' ')
	for i in range(0,2884):
		lencol=0
		if i not in unmapped:
			lenrow=lenrow+1
			col1=[]
			for j in range (0,2884):
				if j not in unmapped:
					lencol=lencol+1
					item1=csv_list[0:,i]
					item2=csv_list[0:,j]
					result=calculate_sym(item1,item2)
					res_str=str(result)
					col1.append(res_str)
			if len(col1)==2647:
				print "ok"
			else:
				print "exception"
			spamwriter.writerow(col1)
if len(lenrow)!=2647:
	print "exception in row"
		

