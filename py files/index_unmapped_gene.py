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
	input_file='yeast.genes.txt'
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
print len(mapped)
print len(unmapped)
input_file='norm_yst_microarray_a.txt'
dataset = pd.read_csv(input_file,sep=" ")
csv_list = dataset.iloc[0:,0:].values
output_file='view_dist_matix.csv'
with open(output_file,'w') as outputcsv_file:
	spamwriter = csv.writer(outputcsv_file,delimiter=' ')
	for i in range(0,2884):
		if i not in unmapped:
			col1=[]
			for j in range (0,2884):
				if j not in unmapped:
					item1=csv_list[0:,i]
					item2=csv_list[0:,j]
					result=calculate_sym(item1,item2)
					res_str=str(result)
					col1.append(res_str)
			print len(col1)
			spamwriter.writerow(col1)
