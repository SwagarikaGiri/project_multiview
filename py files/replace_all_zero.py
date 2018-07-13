import pandas as pd
import csv
import math
import numpy
#initially we need take the similarity matrix 
input_file='r12dst.txt'
dataset = pd.read_csv(input_file,sep="\t",header= None)
csv_list = dataset.iloc[0:,0:2647].values
# once we get the dissimilarity matix then we do the following operation take it as input and find zero rows
input_file_dissimilarity='dissimilarity_view.txt'
dataset_dissimilarity = pd.read_csv(input_file_dissimilarity,sep=" ",header= None)
csv_list_dissimilarity = dataset_dissimilarity.iloc[0:,0:].values

#view euclian distance matrix
input_file_eucledian='view1_new.txt'
dataset_eucledian = pd.read_csv(input_file_eucledian,sep=" ",header= None)
csv_list_eucledian = dataset_eucledian.iloc[0:,0:].values

index_with_all_zero=[]
# view dissimilarity 
output_dissim='view_dissimiraity_matrix_2545.csv'
output_eucledian='view_eucledian_dist_2545.csv'

def convert_into_dissimilarity(output_file,csv_list):
	with open(output_file,'w') as outputcsv_file:
		spamwriter = csv.writer(outputcsv_file,delimiter=" ")
		for i in range(0,2647):
			col1=[]
			for j in range(0,2647):
				value=csv_list[i,j]
				new_value=float(1-value)
				col1.append(round(new_value,7))
			spamwriter.writerow(col1)
			print len(col1)


# output_file='dissimilarity_view_t_0.csv'
# convert_into_dissimilarity(output_file,csv_list)

def indentify_index_of_zero_row(csv_list_dissimilarity):
	for i in range(0,2647):
		num_zero=0
		for j in range(0,2647):
			if csv_list_dissimilarity[i,j]==0:
				num_zero=num_zero+1
		if num_zero==2647:
			index_with_all_zero.append(i)

# here we identify the zero row 
# total 102 rows are there which are zero this will result in elimination of those rows from the column too as they are zero column
indentify_index_of_zero_row(csv_list_dissimilarity)
# print index_with_all_zero
# print len(index_with_all_zero)
def store_zero_and_nonzero(zero_output,non_zero,index_with_all_zero):
	with open(zero_output,'w') as outputcsv_file1:
		spamwriter_zero = csv.writer(outputcsv_file1,delimiter=' ')
		with open(non_zero,'w') as outputcsv_file2:
			spamwriter_nonzero = csv.writer(outputcsv_file2,delimiter=' ')
			for i in range(0,2647):
				col1=[]
				if i in index_with_all_zero:
					col1.append(i)
					spamwriter_zero.writerow(col1)
				else:
					col1.append(i)
					spamwriter_nonzero.writerow(col1)

def create_file_eliminate(csv_list,index_with_all_zero,output_file):
	with open(output_file,'w') as outputcsv_file:
		spamwriter = csv.writer(outputcsv_file,delimiter=' ')
		for i in range(0,2647):
			if i not in index_with_all_zero:
				col1=[]
				for j in range(0,2647):
					if j not in index_with_all_zero:
						col1.append(csv_list[i,j])
				print len(col1)
				spamwriter.writerow(col1)
# csv_list_dissimilarity is the dissimilarity matrix we have
# store_zero_and_nonzero('zero.txt','non_zero.txt',index_with_all_zero)
# create_file_eliminate(csv_list_dissimilarity,index_with_all_zero,output_dissim)
# create_file_eliminate(csv_list_eucledian,index_with_all_zero,output_eucledian)
microarray_input='microarray_2647_17.txt'
dataset_microarray= pd.read_csv(microarray_input,sep=" ",header= None)
csv_list_microarray = dataset_microarray.iloc[0:,0:].values
def microarray_ele_row(csv_list_microarray,output_file,index_with_all_zero):
	with open(output_file,'w') as outputcsv_file1:
		spamwriter = csv.writer(outputcsv_file1,delimiter=' ')
		for i in range(0,2647):
			if i not in index_with_all_zero:
				print i
				col1=[]
				for j in range(0,17):
					col1.append(csv_list_microarray[i,j])
				print len(col1)
				spamwriter.writerow(col1)
microarray_ele_row(csv_list_microarray,'microarray_2545_17.txt',index_with_all_zero)






