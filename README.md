To run the code u need gcc compiler 
type
g++ code_phase2.c -o outputfile_name
./outputfile_name view1 view2 name_of_output_file dimension_of_the_similarity_matix
ex:
./output view1.txt view2.txt result12 3000
where 3000 is the dimension of sililarity matrix



#############################################################################################
file_name                   specification
#############################################################################################
check_size.py               it was used to check if the array read is of proper size or not
                            also if u need a subset to test ur result
                            diagonal '0' as dissimilarity matrix had non-zero diagonal value
                            converted into percentage 
_____________________________________________________________________________________________

code_phase2.c             it is the actual amosa code priority given to view2 in case 
                          of conflict
______________________________________________________________________________________________
compute_silhotte_index.py  here the cluster result thus formed we are checking the validity 
                           of it from the actual dataset here actual dataset is have 17      
                           dimensional vectors
_____________________________________________________________________________________________
dissimilarity_view.txt    it is the dissimilarity between the gene 2647 genes are present

_____________________________________________________________________________________________
distance.cpp and distance.h they are the file used along with silhoutte for distance calculation

_______________________________________________________________________________________________
file1.txt and file2.txt were initial 30*30 file that were used for trial

_______________________________________________________________________________________________

index_unmapped_gene.py         there were 2884 gene from which 237 were unmapped hence     
                               needed to be removed  
_______________________________________________________________________________________________

interchange_rows_col.py      if we need to interchange rows with column and column with rows


____________________________________________________________________________________________
mapped_unmapped.py       same seperates mapped gene from unmapped

______________________________________________________________________________________________
microarray_2545_17.txt    a matrix as a result of gene annotation data has 2545 rows 17 columns obtained as a result of eliminating non mapped and zero row and column

______________________________________________________________________________________________
microarray_2647.txt   gene gene eucledian distance from gene annotation data

_____________________________________________________________________________________________
micrroarray_2647_17.txt   it has result of mapped gene and its 17 dimension

_____________________________________________________________________________________________
microarray_2884.txt       it is 2884* 2884 gene gene eucludian distance obtained from gene annotation data
____________________________________________________________________________________________
new_view_dissimilarity_2545_per.txt   here as the gene gene dissimilarity was in percentage converted it into value multipying with 100

____________________________________________________________________________________________
non_zero.txt         here the gene gene dissimilarity had almost 102 zero columns that caused problem when we had compute membership function hence here store the index of those non-zero for further use

_____________________________________________________________________________________________
norm_yst_microarray_a.txt   it has 2884*17 the actual microarray data the corresponding gene are present in yeast_genes.txt

____________________________________________________________________________________________
r12              it is wrong data 

________________________________________________________________________________________________
r12dst.txt   it is gene gene similarity data but here doesnot have diagonal as 1 so changes had to be made

_______________________________________________________________________________________________
replace_all_zero.py   here we meke changes in the file where the  we needed to remove the all zero rows
we make the following changes here store zero and non zero index
analyse the zero row and non zero row
also eliminate those rows columns where there were all zero conditions

________________________________________________________________________________________________

row_col.py      is used to calculate the eucledian distance from necessary row row combination

_________________________________________________________________________________________________
silhoutte.cpp the code for silhoutte

_________________________________________________________________________________________________
unmapped.txt & unmapped237genesforyeast.txt  these are the genes which were not mapped

_________________________________________________________________________________________________
view_dissim_500.txt  here it is the subset of 500 genes for trial

____________________________________________________________________________
view_dissim_diag_0.txt  dissimilarity matrix with   0 in the diagonal as gene gene dissimilarity should be 0

_____________________________________________________________________________
view_dissimilarity.txt it is the view of dissimilarity matrix 
____________________________________________________________________
  

__________________________________________________________________________


