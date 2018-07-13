import pandas as pd
import csv
import math
def calculate_sym(item1,item2):
	val=0
	for i in range(0,len(item1)):
		val=((item1[i]-item2[i])**2)+val
	return round(math.sqrt(val),7)

input_file='norm_yst_microarray_a.txt'
dataset = pd.read_csv(input_file,sep=" ",header=None)
csv_list = dataset.iloc[0:,0:2884].values
output_file='microarray_2884.csv'
with open(output_file,'w') as outputcsv_file:
	spamwriter = csv.writer(outputcsv_file,delimiter=' ')
	for i in range(0,2884):
		col1=[]
		for j in range (0,2884):
			item1=csv_list[0:,i]
			item2=csv_list[0:,j]
			result=calculate_sym(item1,item2)
			res_str=str(result)
			col1.append(res_str)
		print col1
		spamwriter.writerow(col1)
		


