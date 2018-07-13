import pandas as pd
import csv
import math
mapped=[]
unmapped=[]
input_file='norm_yst_microarray_a.txt'
dataset = pd.read_csv(input_file,sep=" ",header=None)
csv_list = dataset.iloc[0:,0:2884].values
# print csv_list.shape
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
	# return list_
	for i in range(0,2884):
		if (csv_list[i,0] in list_):
			unmapped.append(i)
		else:
			mapped.append(i)
find_mapped_unmapped_index()
# print mapped
# print len(unmapped)
rowlen=0
output_file='mapped_microarray_2647_17_new.csv'
with open(output_file,'w') as outputcsv_file:
	spamwriter = csv.writer(outputcsv_file,delimiter=' ')
	for i in range(0,2884):
		col1=[]
		if i in mapped:
			rowlen=rowlen+1
			print i
			item1=csv_list[0:,i]
			for j in range(0,len(item1)):
				col1.append(str(item1[j]))
			# print len(item1)
			spamwriter.writerow(item1)
print rowlen
