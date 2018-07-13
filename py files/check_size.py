import pandas as pd
import csv
import math
import numpy



def convert_into_percentage(csv_list,input_list):
	with open(csv_list,'w') as outputcsv_file1:
		spamwriter = csv.writer(outputcsv_file1,delimiter=' ')
		for i in range(0,2545):
			col1=[]
			for j in range(0,2545):
				val=input_list[i,j]
				val=float(val*100);
				col1.append(round(val,7))
			spamwriter.writerow(col1)


def subset(out_list,input_list):
	with open(out_list,'w') as outputcsv_file1:
		spamwriter = csv.writer(outputcsv_file1,delimiter=' ')
		for i in range(0,500):
			col1=[]
			for j in range(0,500):
				val=input_list[i,j]
				col1.append(round(val,8))
			print len(col1)
			spamwriter.writerow(col1)
def diagonal_zero(out_list,input_list):
	with open(out_list,'w') as outputcsv_file1:
		spamwriter = csv.writer(outputcsv_file1,delimiter=' ')
		for i in range(0,2545):
			col1=[]
			for j in range(0,2545):
				if i==j:
					val=0
				else:
					val=input_list[i,j]
				col1.append(round(val,8))
			print len(col1)
			spamwriter.writerow(col1)

# convert_into_percentage('view_dissimilarity_2545_per.txt',csv_list1)
# subset('view_dissim_500.txt',csv_list1)
# subset('view_eucledian_500.txt',csv_list2)

input_file='microarray_2647_17.txt'
dataset = pd.read_csv(input_file,sep=" ",header= None)
csv_list1 = dataset.iloc[0:,0:].values
print csv_list1.shape

# input_file="view_eucledian_dist_2545.txt"
# dataset = pd.read_csv(input_file,sep=" ",header= None)
# csv_list2 = dataset.iloc[0:,0:].values
# print csv_list2.shape
# print csv_list1[5,10]
# print csv_list1[10,5]
# diagonal_zero('new_view_dissimilarity_2545_per.txt',csv_list1)


