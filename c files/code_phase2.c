 /*g++ Amosa_new_new.cpp -o amosa -I/home/iitp/ANN/include -L/home/iitp/ANN/lib -lANN*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<limits.h>
#include<string.h>
#include<time.h>
#define DIMENSION 300 //may get changed
#define randm(num) 	(rand()%(num))
#define DUMMY 99999.99
#define DUMMY1 99999
#define MIN_LEN 2
#define MAX_LEN 60
#define MAX_GEN  20//100
#define MAX 5
#define MINFS 25000.0
#define hardl 10
#define softl 20
#define PMUT 0.8
#define NO_OF_ELM 1000000
#define MUTATE_NORMAL 0.75
#define MUTATE_DELETE 0.15
#define MUTATE_INSERT 0.1
#define irand() rand()%1000
#define frand() irand()/1000.0
#define Max(a,b) a>b?a:b
#define Min(a,b) a<b?a:b
#define size 400


struct pattern
{
double x[DIMENSION];
};

int knear=2;
double	eps= 0.00;			// error bound
int maxPts= 3000;
double epsilon=0.00;
double threshold=0.8;

double ***distance_matrix;

int track_count=0.0;
FILE *fp_track;
char func_name[100],mutat[100];

int n,d,r,K,no_testpoint,total_pt;
int arcsize;
double Temperature, TMin = 0.00001, TMax = 100, alpha = 0.75;
int  MaxLen, MinLen;
double **area;
int POP_SIZE;
char file[30];
double *Obj;
double *Obj_current;
int total_index;
int *FitnessIndex;
double **points, **testpoints;
double *min,*max;
FILE *fpo,*testp;
int *fclass, *testclass;
struct archive_element
{
	int len;
	struct pattern z[MAX_LEN];
	int index1[MAX_LEN];
	double ***Membership;
	int **Membership_int;
};
struct archive_element *archive;
struct archive_element current;
struct archive_element new_pool1;
struct archive_element pool[softl];
struct archive_element new_pool;
//struct archive_element pur[pursize];
int current_in_arc,pos,generations;
long seed;


double *euc_dist;
double *range_min, *range_max;
double Dunn_index,*Dunn_pur;
int *views,no_view,clas;

double minkwski(char *);
void main_process();
void mutation();
void WriteChromosome();
//void ComputeMembership_sym(struct archive_element *Chrom);
void ComputeMembership(struct archive_element *Chrom,int);
//void ComputeMembership(struct archive_element *Chrom);
void CombineCenter(struct archive_element *Chrom);
double ComputeFitnessPBM(struct archive_element*,int);
double ComputeFitnessDAI(struct archive_element*);
double ComputeFitnessXB(struct archive_element*,int);
double find_dom(double *,double *);
void burn_in_period();
void clustering();
double FindDistance(double *, double *,int);
int flip(double);
int see_similar(struct archive_element);
int similarity_in_points();
void UpdateCenter(struct archive_element *,int);
//void UpdateCenter(struct archive_element *);
void form_nondominated_archive();
void menu(void) ;
int compute_length();
double find_dom1(double *f1,double *f2);
void printarchive();
void process_current_dominates_new(void);
void process_new_dominates_current(void);
void process_new_current_nondominating(void);
void print_archive1();
void read_file(char *,char *);
void writechromosome(struct archive_element);
void InitializePop();
void printclustering(FILE *c_ptr,struct archive_element *y);
void calculate_len(struct archive_element *,int);

int doc_count;

main(int argc, char *argv[]) // argv[1]= distance using tf-idf argv[2]= distance using scp argv[3]= output file name argv[4] = no of documents
{

	int i,m,j,k,count,dn,*flag2,ii;
	int dim,p;
	FILE *c_ptr, *plt;
	char c_file_name[600],string[2][500];
	double maxnearest_dist,mindist,dist,ob,minkow;
	POP_SIZE=softl;
	flag2=(int *)malloc(POP_SIZE*sizeof(int));
	strcpy(string[0],argv[1]);
	strcpy(string[1],argv[2]);
	doc_count=atoi(argv[4]);
	n=doc_count;
	printf("\n**** doc count= %d ****\n",n);

	read_file(string[0],string[1]);



	dim=d;
	strcpy(func_name,"name");
	//strcpy(file,argv[3]);
	//euc_dist=(double *)malloc(n*sizeof(double));


	/*for(i=0;i<n;i++)
	{
			printf("\n%d-  ",i+1);
		   for(j=0;j<d;j++)
		    {
		      printf("%lf\t",points[i][j]);

		    }
	}
	printf("\n");

	for(j=0;j<d;j++)
	{
	    printf("%d\t",views[j]);

	}*/

	/*kdTree = new ANNkd_tree(					// build search structure
						dataPts,		// the data points
						nPts,			// number of points
						dim);			// dimension of space



	maxnearest_dist=0;
	printf("\n\n here\n\n");
	for(i=0;i<n;i++)
	{

	       mindist=999999;
	       for(j=0;j<n;j++)
	        {
			   if(i!=j)
			    {
				     dist=0;
				     for(p=0;p<d;p++)
				      	dist=dist+pow((points[i][p]-points[j][p]),2);

				     dist=sqrt(dist);

				     if(mindist>dist)
				      	mindist=dist;
			    }
		 	}

		    if(maxnearest_dist<mindist)
		          maxnearest_dist=mindist;
	}*/
	//  threshold=maxnearest_dist;
	threshold=1.2;
	//printf("\n Max nearest distance=%lf\n\n",maxnearest_dist);
	printf("\n Enter maximum number of clusters ");
	//scanf("%d",&MaxLen);
	MaxLen=sqrt(n);
	//printf("%d\n",MaxLen);
	//printf("\n Enter minimum number of clusters ");
	//scanf("%d",&MinLen);
	MinLen=2;
	//printf("%d\n",MinLen);

	fp_track=fopen("FindDistance","w");

	//doc_count=atoi(argv[4]);

	archive=(struct archive_element *)malloc((softl)*sizeof(struct
							archive_element));

	for(i=0;i<(softl);i++)
	{
		if(no_view>1)
		{
			if ( (archive[i].Membership = (double ***) malloc( (no_view+1)* sizeof (double **))) == NULL)
					exit(0);
			if ( (pool[i].Membership = (double ***) malloc( (no_view+1) * sizeof (double **))) == NULL)
				exit(0);
		}
		else
		{
			if ( (archive[i].Membership = (double ***) malloc(sizeof (double **))) == NULL)
					exit(0);
			if ( (pool[i].Membership = (double ***) malloc(sizeof (double **))) == NULL)
				exit(0);
		}

		for(k=0;k<=no_view;k++)
		{
			if ( (archive[i].Membership[k] = (double **) malloc( MaxLen* sizeof (double *))) == NULL)
				exit(0);
			if ( (pool[i].Membership[k] = (double **) malloc( MaxLen* sizeof (double *))) == NULL)
				exit(0);
			for (j=0;j<MaxLen; j++)
			{

				//printf("\nEnd allocating Membership %d of size %d: ",i,n);
				if ( (archive[i].Membership[k][j] = (double *) malloc (n* sizeof(double))) == NULL)
					exit(0);
				if ( (pool[i].Membership[k][j] = (double *) malloc (n* sizeof(double))) == NULL)
					exit(0);
			}
		}
	}


	//printf("\nEnter how many indices u want to optimize(value of total index=3 Enter 3):");
	//scanf("%d",&total_index);
	if(no_view>1)
		total_index=3;
	else
		total_index=2;
	FitnessIndex=(int *)malloc(total_index*sizeof(int));
	range_min=(double *)malloc(total_index*sizeof(double));
	range_max=(double *)malloc(total_index*sizeof(double));

	//menu();
	if(no_view>1)
		new_pool.Membership=(double ***)malloc((no_view+1)*sizeof(double **));
	else
		new_pool.Membership=(double ***)malloc(sizeof(double **));

	for(j=0;j<=no_view;j++)
	{
		new_pool.Membership[j]=(double **)malloc(MaxLen*sizeof(double *));
		for(i=0;i<MaxLen;i++)
				new_pool.Membership[j][i]=(double *)malloc(n*sizeof(double));
	}
	/*for(i=0;i<total_index;i++)
	{
		printf("\nEnter %dth index: ",i);
		scanf("%d",&FitnessIndex[i]);
	}*/
	FitnessIndex[0] = 1;
	FitnessIndex[1] = 2;
	FitnessIndex[2] = 3;
	//printf("\nstart");
	//printf("\nEnter seed : ");
	//scanf("%lf",&seed);

	srand((unsigned)time(NULL));

	InitializePop();

	//printf("\n After This Initialization phase");
	//printf("running 5");
	//WriteChromosome();
	for(m=0;m<POP_SIZE;m++)
	{
		//printf("\n m=%d",m);
		if(no_view>1)
		{
			for(i=0;i<30;i++)
			{
				track_count=i;
				ComputeMembership(&(pool[m]),1);

				//printf("Membership computed for view 1\n\n");
				UpdateCenter(&(pool[m]),1);
printf("\n4\n");


			}

			//getchar();
			for(i=0;i<30;i++)
			{
				    track_count=i;
				ComputeMembership(&(pool[m]),2);
				//printf("Membership computed for view 2\n\n");
				UpdateCenter(&(pool[m]),2);
printf("\n04\n");

			}

			track_count=0;
			CombineCenter(&(pool[m]));

		}
		else
		{
			for(i=0;i<30;i++)
			{
				track_count=i;
				ComputeMembership(&(pool[m]),0);


				UpdateCenter(&(pool[m]),0);


			}
		}

	}
	//printf("\n After Updation");


	form_nondominated_archive();
	printf("\nout form_nondominated_archive\n");
	//getchar();

	POP_SIZE=hardl;
	//printf("\n After initialization phase");
	//printf("running 11");
	printarchive();
	//printf("\n arcsize=%d",arcsize);
	Temperature = TMax;
	//printf("\n tmax=%lf",TMax);
	//alpha=exp((log(TMin)-log(TMax))/MAX_GEN);
	//printf("\n alpha=%lf",alpha);
	//printf("running 12");
	/*printf("\n %lf",Temperature);*/
	//printf("\n arcsize=%d",arcsize);
	if(arcsize==1)
	   r=0;
	else
	   r=rand()% arcsize;
	//printf("\n r=%d",r);
	generations=1;
	//printf("running 13");
	for(j=0;j<MaxLen;j++)
	{
		for(dn=0;dn<d;dn++)
			current.z[j].x[dn]= archive[r].z[j].x[dn];


	}
	//current.len=archive[r].len;
	if(no_view>1)
		current.Membership=(double ***)malloc((no_view+1)*sizeof(double **));
	else
		current.Membership=(double ***)malloc(sizeof(double **));
	for(j=0;j<=no_view;j++)
	{
		current.Membership[j]=(double **)malloc(MaxLen*sizeof(double *));
		for(i=0;i<MaxLen;i++)
			current.Membership[j][i]=(double *)malloc(n*sizeof(double));
	}

	for(k=0;k<=no_view;k++)
	{
		for(i=0;i<MaxLen;i++)
		{
			for(j=0;j<n;j++)
				current.Membership[k][i][j]=archive[r].Membership[k][i][j];
		}
	}

	Obj_current=(double *)malloc(total_index*sizeof(double));
	Obj=(double *) malloc(total_index*sizeof(double));
	for(k=0;k<total_index;k++)
		Obj_current[k]=area[r][k];
	current_in_arc=1;
	pos=r;
	while (Temperature > TMin)
	{
		printf("\n In the main loop");
		for (generations = 0; generations < MAX_GEN; generations++)
		{

				/*mutation_string();  /*  from pool to new_pool */
				mutation();
				for(i=0;i<30;i++)
				{
						if(no_view>1)
						{
					        ComputeMembership(&new_pool,1);
							UpdateCenter(&new_pool,1);
printf("\n5\n");
							ComputeMembership(&new_pool,2);
printf("\n6\n");
							UpdateCenter(&new_pool,2);
							printf("\n1\n");
						    CombineCenter(&new_pool);
							printf("\n01\n");
						}
						else
						{
							ComputeMembership(&new_pool,0);
							UpdateCenter(&new_pool,0);
						}

			    }
			    //CombineCenter(&new_pool);
				for(k=0;k<total_index;k++)
				{
				    if(no_view>1)
				    {
						if (FitnessIndex[k] == 1) Obj[k] = ComputeFitnessPBM(&new_pool,1);
						else if (FitnessIndex[k] == 2) Obj[k] = ComputeFitnessPBM(&new_pool,2);
				        else if (FitnessIndex[k] == 3) Obj[k] = ComputeFitnessDAI(&new_pool);
				    }
				    else
				    {
				    	if (FitnessIndex[k] == 1) Obj[k] = ComputeFitnessXB(&new_pool,0);
						else if (FitnessIndex[k] == 2) Obj[k] = ComputeFitnessPBM(&new_pool,0);
				    }
			    }

				main_process();
		}
			Temperature=Temperature*alpha;
			printf("\n Temperature=%lf",Temperature);
	}
	printf("\n At end");

	printf("\narcsize=%d\n",arcsize);
	sprintf(c_file_name,"%s_plot_multiview.out",argv[3]);
	plt=fopen(c_file_name,"w+");
	for(ii=0;ii<arcsize;ii++)
	{
	               // ComputeMembership_sym(&archive[ii]);
	/*	if(no_view>1)
	    	CombineCenter(&archive[ii]);*/
        /*ComputeMembership(&archive[ii],0);
	    UpdateCenter(&archive[ii],0);
	   archive[ii].len=0;
	  printf("\n");
	    for (k = 0; k < MaxLen;k++)
	    {
	    	if(archive[ii].z[j].x[0] < DUMMY)
	    	{
	    		archive[ii].len++;
	    		printf("%lf",archive[ii].z[j].x[0]);
	    	}
	    }
	    printf("\n");*/
		
		calculate_len(&archive[ii],0);
    
		for(k=0;k<total_index;k++)
				fprintf(plt,"%lf ",area[ii][k]);

		fprintf(plt,"\n");

		sprintf(c_file_name,"%s_clus%d_multiview",argv[3],ii+1);

		printf("\nWriting file clus%d\n",ii+1);
		c_ptr = fopen(c_file_name,"w");
		fprintf(c_ptr,"%d  %d  %d\n",n,d,archive[ii].len);
		printclustering(c_ptr,&archive[ii]);
		fclose(c_ptr);
            /* minkow= minkwski(c_file_name);
                c_ptr=fopen(c_file_name,"a+");
                     fprintf(c_ptr,"\n Minkowski Score=%lf",minkow);
               fclose(c_ptr);*/


	}
	fclose(plt);
//	print_archive1();

}


void print_archive1()
{
	int i,j,k,m,h,dn,ii,*computed_class,Index2,flag,minpoint,Di3,jj,cl,count,p,minindex, point_in_cluster[NO_OF_ELM];
	double mk, *Di2,min,max;//,minkow,minminkow;
	int *flag1;
	char file1[30];
	FILE *c_ptr;
	//strcpy(file1,file);
	//strcat(file1,"_not_multicenter_amosa_mc_sym_XB_vga.txt");
	fpo=fopen(file1,"w+");
	if(fpo==NULL)
	{
		printf("\n Error is opening output file");
		exit(1);
	}
	printf("\n arcsize=%d",arcsize);
	for(i=0;i<arcsize;i++)
	{
		if ( (archive[i].Membership_int = (int **) malloc( (MaxLen)* sizeof (int *))) == NULL)
			exit(0);
	 	for (j=0;j< (MaxLen);j++)
	 		if ( (archive[i].Membership_int[j] = (int *) malloc( n* sizeof (int))) == NULL)
				exit(0);
	}

	for(i=0;i<arcsize;i++)
	{
		for(j=0;j<MaxLen;j++)
		{
			 if(archive[i].z[j].x[0]!=DUMMY)
			 {
				 for(p=0;p<d;p++)
				  	fprintf(fpo," %lf\t",archive[i].z[j].x[p]);
				 fprintf(fpo," \n");
			 }
		 }

	    for(j=0;j<(MaxLen);j++)
		{
		    for(k=0;k<n;k++)
			      archive[i].Membership_int[j][k]=0;
		}
	    count=0;

		for(j=0; j<archive[i].len; j++)  /* j represent the current cluster */
		{

			for (p=0;p<n;p++)
			{
			   if (archive[i].Membership[0][k][p]==1)
				{
					archive[i].Membership_int[count][p]=1;
					point_in_cluster[p]=count;
				}
			}
			count++;
	    }
	}

	 //  minminkow=99.9;
	for(i=0;i<arcsize;i++)
	{
		sprintf(file1,"%s_clus%d_amosa_mc_XB_sym_vga",file,i+1);

		printf("\nWriting file clus%d\n",i+1);
		c_ptr= fopen(file1,"w");
		h=0;
		fprintf(fpo,"\n Length %d: ",archive[i].len);


		for(ii=0;ii<217;ii++)
	    {
	      	for(j=0;j<181;j++)
	               	fprintf(c_ptr,"%d\t",(int) ((255.00/(count-1))*(point_in_cluster[ii*181+j])));

	       	fprintf(c_ptr,"\n");
	 	}

		fclose(c_ptr);


	}


		//fprintf(fpo,"\n Minimum Minkowski Score=%lf, solution number=%d",minminkow,minindex);
	fclose(fpo);

}

double find_dom1(double *f1,double *f2)
{
	   int i;
	   double amo=1;
	   for(i=0;i<total_index;i++)
	    {
	      if(fabs(f1[i]-f2[i])!=0)
	        {
	         	amo=amo*(f1[i]-f2[i])/(range_max[i]-range_min[i]);
		 	}

	    }
	   return(amo);
}

/* The below function is main process */
void main_process()
{
  int count1=0,count2=0,k;
  for(k=0;k<total_index;k++)
    {
         if(Obj_current[k]<=Obj[k])
	        count1++;

         if(Obj_current[k]>=Obj[k])
	       	count2++;
    }
  if(count1==total_index)
    {

		printf("\n Current dominates the newly generated one");
		process_current_dominates_new();
    }
  else if(count2==total_index)
    {
		printf("\n New solution dominates the current one");

		process_new_dominates_current();
    }
  else
    {
		printf("\n New and Current are non dominating to each other");
		process_new_current_nondominating();
    }
}


void mutation()
{

	double f,del;
	double rnd,rand_lap;
	int i,j,k,dn,pick1,m,pick2,j2,MutationOccurred = 1,l, index,flag=0,count,min_clus,pick,*picks,p;
	int **points_in_clus,max_clus,clus1,clus2,index2,jj;
	int *total_pt_in_clus;
	double *intra_distance,max_intra,s2,minimum,dist,sum,*sum1;
	double mut_scal=0.2, b=12.5;

	picks=(int*)calloc(MaxLen,sizeof(int));

	for(i=0;i<MaxLen;i++)
		picks[i]=DUMMY;

	for(i=0;i<MaxLen;i++)
	{
		for(dn=0;dn<d;dn++)
		{
			if(current.z[i].x[dn]!=DUMMY)
				picks[i] = current.z[i].x[dn];
		}
	}
	new_pool.len=0;
	for(j=0;j<MaxLen;j++)
	{
		for(dn=0;dn<d;dn++)
		{
		        // printf("\t %lf",current.z[j].x[dn]);
			new_pool.z[j].x[dn]= current.z[j].x[dn];
		}

		if(new_pool.z[j].x[0] < DUMMY)
			new_pool.len++;
	}
	/*printf("\n+++++In mutation before+++++\n");
	for(j=0;j<MaxLen;j++)
	{  if(new_pool.z[j].x[0]!=DUMMY)
            {
				for(dn=0;dn<d;dn++)
				{
							printf("%lf ",new_pool.z[j].x[dn]);
				}

				printf("\n");

			}

	}
	getchar();*/


	for(j=0;j<MaxLen;j++)
	{


            if(new_pool.z[j].x[0]!=DUMMY)
			{
                if(flip(PMUT)) {

	     		    for(dn=0;dn<d;dn++)
	     		    {
					   if (flip(MUTATE_NORMAL))
					   {


						   	do
						   	{
						   		pick = rand()%n;
						   		flag=0;
						   		for(p=0;p<MaxLen;p++)
						   		{
						   			if(pick==picks[p])
						   				flag=1;
						   		}
						   	}while(flag);

						   	new_pool.z[j].x[dn]= points[pick][dn];
					  	}
				    }
				}

            }

				else  if (flip(MUTATE_DELETE) && (new_pool.len>2))
				      {		//strcpy(mutat,"MUTATE_DELETE");
	     		    		for(dn=0;dn<d;dn++)
	         					new_pool.z[j].x[dn]=DUMMY;
				      }


			else  if(flip(MUTATE_INSERT))
	      		{	//strcpy(mutat,"MUTATE_INSERT");
					pick = rand()%n; /*select an element */
	     		    for(dn=0;dn<d;dn++)
	         			new_pool.z[j].x[dn]=points[pick][dn] ;
	      		}
                  }

		/*printf("\n+++++%s+++++\n",mutat);
		for(dn=0;dn<d;dn++)
		{
		       if(new_pool.z[j].x[0]!= DUMMY)
					printf("%lf ",new_pool.z[j].x[dn]);
		}

	printf("\n+++++In mutation after+++++\n");
	for(j=0;j<MaxLen;j++)
	{
		for(dn=0;dn<d;dn++)
		{
		       if(new_pool.z[j].x[0]!= DUMMY)
					printf("%lf ",new_pool.z[j].x[dn]);
		}

		printf("\n");

	}
	getchar();*/
}




double find_dom(double *f1,double *f2)
 {
   int i;
   double amo=1;
   for(i=0;i<total_index;i++)

      if(fabs(f1[i]-f2[i])!=0)

	        amo=amo*fabs(f1[i]-f2[i])/(range_max[i]-range_min[i]);/*amount of domination type 2*/

   return(amo);
}

void process_current_dominates_new(void)
{
	int dominated_by=0,i,count,j,k,u,v;
	double deldom=0.0,product=1,prob,ran,amount;
	deldom=0.0;
	for(i=0;i<arcsize;i++)
	{
	    count=0;
	    for(k=0;k<total_index;k++)
	    {
		        if(area[i][k]<=Obj[k])
					count++;
	    }
	    if(count==total_index)
		{
			 dominated_by++;
			 amount=find_dom(Obj,area[i]);
			 deldom=deldom+amount;
		}

	}

	if(current_in_arc==0)
	{
		amount=find_dom(Obj,Obj_current);
		deldom=deldom+amount;
		dominated_by++;
	}
	prob=1.0/(1.0+exp(deldom/(dominated_by*Temperature)));

	ran=frand();
	if(prob>=ran)
	{
	    for(j=0;j<MaxLen;j++)
	    {
	        current.z[j]=new_pool.z[j];
	    }

	    for(j=0;j<=no_view;j++)
	    {
		    for(u=0;u<MaxLen;u++)
			{
			    for(v=0;v<n;v++)
			   		current.Membership[j][u][v]=new_pool.Membership[j][u][v];
			}
	    }

	    for(k=0;k<total_index;k++)
		{
			Obj_current[k]=Obj[k];
		}

	    current_in_arc=0;
	}

}

void process_new_dominates_current(void)
{
	int dominated_by,dominates,i,j,k,h,count1,count2,loc,g,u,v;
	double minimum,product,prob,ran,amount;
	int min_point,dn;
	struct archive_element * archive1;
	double ** area1;
	int *flag_dom=(int *)malloc(arcsize*sizeof(int));
	if(current_in_arc==1) flag_dom[pos]=1;
		minimum=9999999;
	dominated_by=0;
	dominates=0;
	printf("\ndom1\n");
	for(i=0;i<arcsize;i++)
	{
		flag_dom[i]=0;
		count1=0;count2=0;
		for(k=0;k<total_index;k++)
		{
			if(area[i][k]<=Obj[k])
				count1++;
			if(area[i][k]>=Obj[k])
				count2++;
		}
		if(count1==total_index)
		{
		   dominated_by++;
		   product=find_dom(Obj,area[i]);

		   if(minimum>product)
			{
				minimum=product;
				min_point=i;
			}
		}
		if(count2==total_index)
		{
		   dominates++;
		   flag_dom[i]=1;
		}
	}
	// end of for
	if(current_in_arc==1)
		dominates--;
	printf("\ndom2\n");
	if(dominated_by==0 && dominates==0)
	{
		printf("\ndom3\n");
		if(current_in_arc==1)
		  loc=pos;
		else
		{
			

		    if(arcsize>=softl)
		    	clustering();
		    loc=arcsize;
		}

		for(j=0;j<MaxLen;j++)
		{
			for(dn=0;dn<d;dn++)
		   	{
			     archive[loc].z[j].x[dn]=new_pool.z[j].x[dn];
			     current.z[j].x[dn]=new_pool.z[j].x[dn];
		    }
		}
		printf("\ndom4\n");
		for(j=0;j<=no_view;j++)
		{
			for(u=0;u<MaxLen;u++)
			{
				for(v=0;v<n;v++)
				{
					archive[loc].Membership[j][u][v]=new_pool.Membership[j][u][v];
					//printf("%d \n",loc);
					fflush(stdout);
					current.Membership[j][u][v]=new_pool.Membership[j][u][v];
				}
			}
		}
		printf("\ndom5\n");
		for(k=0;k<total_index;k++)
		{
			area[loc][k]=Obj[k];
		  	if(range_max[k]< area[loc][k])
					range_max[k]=area[loc][k];
		  	if(range_min[k]> area[loc][k])
			   		range_min[k]=area[loc][k];
		  	Obj_current[k]=Obj[k];
		}
		current_in_arc=1;
		pos=loc;
	}
	//end of if
	else if(dominated_by==0 && dominates>0)
	{
		archive1=(struct archive_element *)malloc((arcsize)*sizeof(struct archive_element));
		area1=(double **)malloc((arcsize)*sizeof(double *));
		for(i=0;i<arcsize;i++)
			area1[i]=(double *)malloc(total_index*sizeof(double));
		printf("\ndom6\n");
		for(h=0;h<arcsize;h++)
		{
			for(j=0;j<MaxLen;j++)
			{
			  	for(dn=0;dn<d;dn++)
		        {
			      archive1[h].z[j].x[dn]=archive[h].z[j].x[dn];
			    }
			}
			if(no_view>1)
				archive1[h].Membership=(double ***)malloc((no_view+1)*sizeof(double **));
			else
				archive1[h].Membership=(double ***)malloc(sizeof(double **));
			for(k=0;k<=no_view;k++)
			{
				archive1[h].Membership[k]=(double **)malloc(MaxLen*sizeof(double *));

				for(i=0;i<MaxLen;i++)
						archive1[h].Membership[k][i]=(double *)malloc(n*sizeof(double));
			}

			for(i=0;i<=no_view;i++)
			{
				for(u=0;u<MaxLen;u++)
				{
				  for(v=0;v<n;v++)
					archive1[h].Membership[i][u][v]=archive[h].Membership[i][u][v];
				}
		    }

		    for(k=0;k<total_index;k++)
		    {
			   area1[h][k]=area[h][k];
			}
		}

		g=0;
		printf("\ndom7\n");
		for(h=0;h<arcsize;h++)
		{
		   if(flag_dom[h]==0)
		    {

				for(j=0;j<MaxLen;j++)
				{
			   		for(dn=0;dn<d;dn++)
			     		archive[g].z[j].x[dn]=archive1[h].z[j].x[dn];
			    }

				for(j=0;j<=no_view;j++)
				{
					for(u=0;u<MaxLen;u++)
					{
					   for(v=0;v<n;v++)
							archive[g].Membership[j][u][v]=archive1[h].Membership[j][u][v];
					}
				}

				for(k=0;k<total_index;k++)
				{
				    area[g][k]=area1[h][k];
				}
				g++;
		  	}
		}

		arcsize=g;
		printf("\ndom8\n");
		for(i=0;i<arcsize;i++)//edit
		{

		    free(area1[i]);
		}

		free(area1);
		if(arcsize>=softl)
			clustering();

		for(j=0;j<MaxLen;j++)
		{
			for(dn=0;dn<d;dn++)
		    {
			        current.z[j].x[dn]=new_pool.z[j].x[dn];
			        archive[arcsize].z[j].x[dn]=new_pool.z[j].x[dn];
			}
		}

		for(j=0;j<=no_view;j++)
		{
		   for(u=0;u<MaxLen;u++)
			{
			    for(v=0;v<n;v++)
			    {
			        archive[arcsize].Membership[j][u][v]=new_pool.Membership[j][u][v];
			        current.Membership[j][u][v]=new_pool.Membership[j][u][v];
			    }
			}
		}

		for(k=0;k<total_index;k++)
		{
		    area[arcsize][k]=Obj[k];
		    if(range_max[k]< area[arcsize][k])
		            range_max[k]=area[arcsize][k];
		    if(range_min[k]> area[arcsize][k])
		            range_min[k]=area[arcsize][k];
		    Obj_current[k]=Obj[k];
		}
		   arcsize++;
		   pos=arcsize-1;
		   current_in_arc=1;
	}
	else if(dominated_by>0 && dominates==0)
	{
	    prob=1.0/(1.0+exp(-minimum));
	    ran=frand();
	    printf("\ndom9\n");
	    if(prob>=ran)
	    {
	        for(j=0;j<MaxLen;j++)
		    {
				for(dn=0;dn<d;dn++)
		        {
		            current.z[j].x[dn]=archive[min_point].z[j].x[dn];
			    }
			}

			for(j=0;j<=no_view;j++)
			{
			    for(u=0;u<MaxLen;u++)
			    {
			        for(v=0;v<n;v++)
					    current.Membership[j][u][v]=archive[min_point].Membership[j][u][v];
			    }
			}

		    for(k=0;k<total_index;k++)
		    {
		        Obj_current[k]= area[min_point][k];
		    }
		    current_in_arc=1;
		    pos=min_point;
		}
		else
		{
		    for(j=0;j<MaxLen;j++)
		    {
		        for(dn=0;dn<d;dn++)
		        {
			   		current.z[j].x[dn]=new_pool.z[j].x[dn];
		        }
		    }

			for(j=0;j<=no_view;j++)
			{
			    for(u=0;u<MaxLen;u++)
				{
					for(v=0;v<n;v++)
						current.Membership[j][u][v]=new_pool.Membership[j][u][v];
				}
			}
		    for(k=0;k<total_index;k++)
		    {
				Obj_current[k]=Obj[k];
		    }
		    current_in_arc=0;
		}
		printf("\ndom10\n");
	}
}

void process_new_current_nondominating(void)
{
	int dominated_by=0,i,count1,count2,k,j,m,g,u,v;
	int dominates=0,dn;
	double deldom=0,product,prob,ran;
	struct archive_element *archive1;
	double **area1 ;
	int *flag_dom=(int *)malloc(arcsize*sizeof(int));
	for(i=0;i<arcsize;i++)
	{
		count1=0;count2=0;
		flag_dom[i]=0;

		for(k=0;k<total_index;k++)
		{

		   if(area[i][k]<=Obj[k]) count1++;
		   if(area[i][k]>=Obj[k]) count2++;
		}
		if(count1==total_index)
		{
			dominated_by++;
			product=find_dom(Obj,area[i]);

		        deldom=deldom+product;

		}
		if (count2==total_index)
		{
			dominates++;
			flag_dom[i]=1;
		}
	}
	//printf("Whats happen 1\n");
	fflush(stdout);
	if(dominated_by>0 && dominates==0)
	{
	    prob=1.0/(1.0+exp(deldom/(dominated_by*Temperature)));
	    ran=frand();
	    if(prob>=ran)
	    {
	        for(j=0;j<MaxLen;j++)
	        {
			    for(dn=0;dn<d;dn++)
			          current.z[j].x[dn]=new_pool.z[j].x[dn];
			}

			for(j=0;j<=no_view;j++)
			{
			    for(u=0;u<MaxLen;u++)
				{
				    for(v=0;v<n;v++)
				        current.Membership[j][u][v]=new_pool.Membership[j][u][v];
				}
			}

		    for(k=0;k<total_index;k++)
			{
			     Obj_current[k]=Obj[k];
			}
		    current_in_arc=0;
		}
	}
	else if(dominated_by==0 && dominates==0)
	{

		if(arcsize>=softl)
		    clustering();
		//printf("Whats happen 2\n");
	    fflush(stdout);
		for(j=0;j<MaxLen;j++)
		{
		    for(dn=0;dn<d;dn++)
		    {
		     	archive[arcsize].z[j].x[dn]=new_pool.z[j].x[dn];
		     	current.z[j].x[dn]=new_pool.z[j].x[dn];
		    }
		}

		for(j=0;j<=no_view;j++)
		{
			for(u=0;u<MaxLen;u++)
			{
				for(v=0;v<n;v++)
				{
					current.Membership[j][u][v]=new_pool.Membership[j][u][v];
					printf("%d \n",arcsize);
					fflush(stdout);
					archive[arcsize].Membership[j][u][v]=new_pool.Membership[j][u][v];
				}
			}
		}
		for(k=0;k<total_index;k++)
		{
			  area[arcsize][k]=Obj[k];
			  if(range_max[k]< area[arcsize][k])
				            range_max[k]=area[arcsize][k];
			  if(range_min[k]> area[arcsize][k])
				           range_min[k]=area[arcsize][k];
			  Obj_current[k]=Obj[k];
		}
		arcsize++;
		current_in_arc=1;
		pos=arcsize-1;
	    //printf("Whats happen 3\n");
	    fflush(stdout);
	}
	else if(dominated_by==0 && dominates>0)
	{
		archive1=(struct archive_element *)malloc((arcsize)*sizeof(struct archive_element));
		area1=(double **)malloc((arcsize)*sizeof(double));
		for(i=0;i<arcsize;i++)
		{
			if(no_view>1)
				archive1[i].Membership=(double ***)malloc((no_view+1)*sizeof(double **));
			else
				archive1[i].Membership=(double ***)malloc(sizeof(double **));
			for(k=0;k<=no_view;k++)
			{
		    	archive1[i].Membership[k]=(double **)malloc(MaxLen*sizeof(double *));
		    	for(j=0;j<MaxLen;j++)
					archive1[i].Membership[k][j]=(double *)malloc(n*sizeof(double));
			}
		}

		   for(k=0;k<arcsize;k++)
		     {
		       area1[k]=(double *)malloc(total_index*sizeof(double));
		     }
		   for(j=0;j<arcsize;j++)
		     {
		        for(m=0;m<MaxLen;m++)
		        {
		           for(dn=0;dn<d;dn++)
		           {
		                archive1[j].z[m].x[dn]=archive[j].z[m].x[dn];
		           }
		        }

		        for(k=0;k<=no_view;k++)
		        {
			       for(u=0;u<MaxLen;u++)
			       {
				        for(v=0;v<n;v++)
						    archive1[j].Membership[k][u][v]=archive[j].Membership[k][u][v];
				   }
				}

		       for(k=0;k<total_index;k++)
		       {
		            area1[j][k]=area[j][k];
		       }
			  // printf("Whats happen 4\n");
			   fflush(stdout);
		}
		g=0;
		for(i=0;i<arcsize;i++)
		{
			if(flag_dom[i]==0)
			{
				for(j=0;j<MaxLen;j++)
		  		{
		   			for(dn=0;dn<d;dn++)
		    		{
		      			archive[g].z[j].x[dn]=archive1[i].z[j].x[dn];
		     		}
		   		}

		   		for(k=0;k<=no_view;k++)
		   		{
					for(u=0;u<MaxLen;u++)
					{
						for(v=0;v<n;v++)
							archive[g].Membership[k][u][v]=archive1[i].Membership[k][u][v];
					}
				}

				for(k=0;k<total_index;k++)
		    		area[g][k]=area1[i][k];
				g++;
			//	printf("Whats happen 5\n");
				fflush(stdout);
			}
		}
		for(i=0;i<arcsize;i++)
		{
		   /* for(u=0;u<MaxLen;u++)
			{
				free(archive1[i].Membership[u]);
			}
		    free(archive1[i].Membership);*/

		    free(area1[i]);
		}
	//	free(archive1);
		free(area1);

		arcsize=g;
		/*if(arcsize>=POP_SIZE && arcsize<softl)
		{
			archive=(struct archive_element *)realloc(archive,sizeof(struct archive_element));
			archive[arcsize].Membership=(double **)malloc(MaxLen*sizeof(double *));
			area=(double **)realloc(area,sizeof(double));
			area[arcsize]=(double *)malloc(total_index*sizeof(double));
		}

		else */if(arcsize>=softl)
		     clustering();
		for(j=0;j<MaxLen;j++)
		{
		  for(dn=0;dn<d;dn++)
		    {
		       archive[arcsize].z[j].x[dn]=new_pool.z[j].x[dn];
		       current.z[j].x[dn]=new_pool.z[j].x[dn];
		    }
		}
		//printf("Whats happen 6\n");
		fflush(stdout);
		if(no_view>1)
			archive[arcsize].Membership=(double ***)malloc((no_view+1)*sizeof(double**));
		else
			archive[arcsize].Membership=(double ***)malloc(sizeof(double**));

		for(k=0;k<=no_view;k++)
		{
			archive[arcsize].Membership[k]=(double **)malloc(MaxLen*sizeof(double*));

			for(u=0;u<MaxLen;u++)
			{
				archive[arcsize].Membership[k][u]=(double *)malloc(n*sizeof(double));
				for(v=0;v<n;v++)
				{
					current.Membership[k][u][v]=new_pool.Membership[k][u][v];
					archive[arcsize].Membership[k][u][v]=new_pool.Membership[k][u][v];
		        }
			}

		}
		//printf("Whats happen 7\n");
		fflush(stdout);
		for(k=0;k<total_index;k++)
		{
			Obj_current[k]=Obj[k];
			area[arcsize][k]=Obj[k];
			if(range_max[k]< area[arcsize][k])
				    range_max[k]=area[arcsize][k];
			if(range_min[k]> area[arcsize][k])
				    range_min[k]=area[arcsize][k];
		}
		arcsize++;
		current_in_arc=1;
		pos=arcsize-1;
	}

}

/*void menu(void)

{


	printf("\n Enter 3 for ComputeFitnessPBM(pool)");
	printf("\n Enter 8 for ComputeFitnessDAI(pool)");

}*/



void read_file(char *file1,char *file2)
{
   char *inputfile;
   FILE *fp,*fp_views;
   int i=0,j=0,p=0,k=0;
 /*	printf("\n%s",file1);
	getchar();*/

   distance_matrix= (double***)calloc(2,sizeof(double**));
   for(i=0;i<2;i++)
   {
   		distance_matrix[i]= (double**)calloc(doc_count,sizeof(double*));
   		for(j=0;j<doc_count;j++)
   		{
   			distance_matrix[i][j]= (double*)calloc(doc_count,sizeof(double));
   		}
   }


   no_view=2;// 0 for single view and >=2 for multiview


   if(no_view>1)
   {

		fp=fopen(file1,"r");
		if(fp==NULL)
		{
			printf("\n%s",file1);
			perror("Error: ");
			exit(1);
		}
		for(i=0;i<doc_count;i++)
			for(j=0;j<doc_count;j++)
				fscanf(fp,"%lf",&distance_matrix[0][i][j]);// distance matrix of tf-idf
		fclose(fp);

		fp=fopen(file2,"r");
		if(fp==NULL)
		{
			printf("\n%s",file2);
			perror("Error1: ");
			exit(1);
		}
		for(i=0;i<doc_count;i++)
			for(j=0;j<doc_count;j++)
				fscanf(fp,"%lf",&distance_matrix[1][i][j]);// distance matrix of scp
		fclose(fp);
   }
   n=doc_count;
   d=1;
   views=(int*) malloc (2*sizeof(int));

   points=(double **)malloc(n*sizeof (double *));
   if(points==NULL)
   {
        printf("\n Error in memory allocation");
	 	exit(0);
   }
   for(i=0;i<n;i++)
   		points[i]=(double *)malloc(d*sizeof(double));
   for(i=0;i<doc_count;i++)
   		points[i][0]=i;

  /* min=(double *)malloc(d*sizeof(double));
   max=(double *)malloc(d*sizeof(double));
   for(j=0;j<d;j++)
   	{
		min[j]=DUMMY;

		max[j]=-DUMMY;
	}
    for(i=0;i<n;i++)
    {
		for(j=0;j<d;j++)
		{
			fscanf(fp,"%lf",&points[i][j]);
                        if(points[i][j]>max[j])  max[j]=points[i][j];
                        else if(points[i][j]<min[j]) min[j]=points[i][j];
		}
		fscanf(fp,"%d",&fclass[i]);

    }*/

    if(no_view>1)
    {
	  /*  for(j=0;j<d;j++)
	   		fscanf(fp_views,"%d",&views[j]);*/
    	views[0]=1;
    	views[1]=2;

    }
     /*for(i=0;i<doc_count;i++)
   {
   		
   		for(j=0;j<doc_count;j++)
   		{
   			printf("%lf \t",distance_matrix[1][i][j]);
   		}
   		printf("\n");
   }
   getchar();*/

}
/*The bellow function is responsible for computing the length of the chormosome */
int compute_length(struct archive_element chrom)
{
	int k, length=0;

	for(k=0;k<MaxLen;k++)
	if(chrom.z[k].x[0]!= DUMMY) length++;
	return (length);
}

int flip(double prob)
{
	double  i;
	i=frand();
	if((prob==1.0) || (i<prob))
	return(1);
	else return(0);
}

void UpdateCenter(struct archive_element *Chrom,int view)
{
	int Index1,dn,i,j=0,pos,count=0;
	double Value, weight=2.0,*arr,min,sum;

	arr= (double*)calloc(doc_count,sizeof(double));

	printf("\n((((((((((((((((before))))))))))))))\n");
	for (Index1 = 0; Index1 < MaxLen; Index1++)
	{

	   		printf("%lf\t",Chrom->z[Index1].x[0]);

	}
	Chrom->len=0;
	for (Index1 = 0; Index1 < MaxLen; Index1++)
	{
		if(Chrom->z[Index1].x[0]!=DUMMY && Chrom->z[Index1].x[0]>=0.0)
	   	{
		  	min=999999.9;
		  	count=0;
		  	Chrom->len++;
		  	//pos = DUMMY;

		  	for(j=0;j<doc_count;j++)
		  	{
		  		arr[j]=Chrom->Membership[view][Index1][j];
		  		if(Chrom->Membership[view][Index1][j]>0.0)
		  			count++;
		  	}
		  	printf("\ncount=%d\n",count);
		  	for(i=0;i<doc_count;i++)
		  	{

		  		if(arr[i]>0.0)
		  		{
		  			sum=0.0;
			  		for(j=0;j<doc_count;j++)
			  		{
			  			if(arr[j]>0.0 && i!=j)
			  			{
				  			if( view==2)
				  				sum+=distance_matrix[1][i][j];
				  			if(view==0 || view==1)
				  				sum+=distance_matrix[0][i][j];
				  		}
			  		}
			  		sum= sum/(double)count;
			  		if(min>sum)
			  		{
				  		//printf("\nmin=%lf sum=%lf\n",min,sum);
				  		min=sum;
				  		pos=i;
			  		}
			  	}

		  	}
		  	Chrom->z[Index1].x[0]=(double)pos;
		}
	}
	printf("\n((((((((((((((((after))))))))))))))\n");
	for (Index1 = 0; Index1 < MaxLen; Index1++)
	{

	   		printf("%lf\t",Chrom->z[Index1].x[0]);
	   		fflush(stdout);

	}
	//free memory-----------------------
	free(arr);
	//----------------------------------
 /*printf("\n---------------Update center-------------------\n");
	for (Index1 = 0; Index1 < MaxLen; Index1++)
	{
		if(Chrom->z[Index1].x[0]!=DUMMY && Chrom->z[Index1].x[0]>=0.0)
	   	{
		 	 	printf("\n%lf",Chrom->z[Index1].x[0]);
		}
	}
	printf("\n----------------------------------\n");*/

}
int count_points(int *a)
{
	int i=0,count=0;
	for(i=0;i<n;i++)
	{
		if(a[i]>0)
			count++;
	}
	return count;
}

int count_points(double *a)
{
	int i=0,count=0;
	for(i=0;i<n;i++)
	{
		if(a[i]>0.0)
			count++;
	}
	return count;
}

void calculate_len(struct archive_element *Chrom,int view)
{
	int count =0,i=0;
	for(i=0;i<MaxLen;i++)
	{
		if(count_points(Chrom->Membership[view][i])>0)
			count++;
			
	}
	Chrom->len=count;
printf("\nlength=%d\n",Chrom->len);
return;
}

void CombineCenter(struct archive_element *Chrom)
{
	int i,j,k,max1=0.0,*temp,pos,Index1=0,dn=0,count1=0,count2=0;
	temp = (int*)calloc(n,sizeof(int));
	for(i=0;i<MaxLen;i++)
		for(j=0;j<n;j++)
			Chrom->Membership[0][i][j]=0.0;

	for(i=0;i<MaxLen;i++)
	{
			if(count_points(Chrom->Membership[1][i])>0)
				count1++;
			if(count_points(Chrom->Membership[2][i])>0)
				count2++;
	}

	for(i=0;i<MaxLen;i++)
	{
		max1=0.0;
		for(j=0;j<MaxLen;j++)
		{
			for(k=0;k<n;k++)
			{
				if(count1>count2)
					temp[k]= (Chrom->Membership[1][i][k] && Chrom->Membership[2][j][k]);
				else
					temp[k]= (Chrom->Membership[1][j][k] && Chrom->Membership[2][i][k]);
			}
			if(max1<count_points(temp))
			{
				max1=count_points(temp);
				for(k=0;k<n;k++)
				{
					Chrom->Membership[0][i][k]=temp[k];
					temp[k]=0;
				}
			}
		}
	}
free(temp);
UpdateCenter(Chrom,0);
printf("\n2\n");
ComputeMembership(Chrom,0);
UpdateCenter(Chrom,0);
printf("\n02\n");

}

void WriteChromosome()

{
	int l,m,dn;
	for(m=0;m<POP_SIZE;m++)
	{
		printf("\n");
		for(l=0;l<MaxLen;l++)
		{

			printf("(");
			for(dn=0;dn<d;dn++)
				printf("%lf  ",pool[m].z[l].x[dn]);
			printf(") ");

		}
	}
}


void printarchive()
{
   int l,m,dn;
 // printf("\n In the initialization phase");
  for(m=0;m<arcsize;m++)
   {

    for(dn=0;dn<total_index;dn++)
      printf("\narea[%d][%d]=%lf\n",m,dn,area[m][dn]);

   }
}


void ComputeMembership(struct archive_element *Chrom,int view)
{
	int Index1, Index2, Index3, flag,minpoint,count,i,j;
	double Sum;
    double *Di2, Di3,min;
   /* Di2= (double*)calloc(MaxLen,sizeof(double));
    if(Di2==NULL)
    	exit(1);
    for (Index1 = 0; Index1 < n; Index1++)
	{
		
		for (Index2 = 0; Index2 < MaxLen; Index2++)
		{
                   Chrom->Membership[view][Index2][Index1]=0.0;
		    if(Chrom->z[Index2].x[0]!=DUMMY)
		    {
		    	fflush(stdout);
		    	printf("\nIndex2=%d Index1=%d view=%d Chrom->z[Index2].x[0]=%lf points[Index1][0]=%lf\n",Index2,Index1,view,Chrom->z[Index2].x[0],points[Index1][0]);
				fflush(stdout);
				i=(int)Chrom->z[Index2].x[0];
				j=(int)points[Index1][0];
				printf("\ni= %d j=%d Chrom->z[Index2].x[0]=%lf \n",i,j,Chrom->z[Index2].x[0]);
				fflush(stdout);
                Chrom->Membership[view][Index2][Index1]=0.0; // assuming view 1 membership is represented as Membership[1][0..class][0..n]
				if(view==0 || view==2)
					Di2[Index2] =  distance_matrix[1][i][j]; //FindDistance(Chrom->z[Index2].x,points[Index1],view);
                if(view==1)
                	Di2[Index2] =  distance_matrix[0][i][j];
                
			}
			printf("\niiindex2=%d index1=%d\n",Index2,Index1);
                fflush(stdout);
			
		}
		printf("\n++++++Di2 value for point %d++++++++\n",Index1);
		fflush(stdout);
		for (Index2 = 0; Index2 < MaxLen; Index2++)
	 	{
	 		printf("%lf\t",Di2[Index2]);
	 		fflush(stdout);
	  	}
	 	printf("\n");
	 	fflush(stdout);
       	min=9999999.99;
       	flag=0;
	 	for (Index2 = 0; Index2 < MaxLen && (!flag); Index2++)
	 	{
	 		printf("\nindex23 = %d",Index2);
			fflush(stdout);
		    if(Chrom->z[Index2].x[0]!=DUMMY)
		    {

                if(!flag)
                {
			       Chrom->Membership[view][Index2][Index1]=0.0;
		    	   if(Di2[Index2]<=0.0)
					{
						Chrom->Membership[view][Index2][Index1]=1.0;
                        flag=1;
			  			minpoint=Index2;
                        break;

					}
		       		else
					{
			  			Di3=Di2[Index2];
			  			if(min>Di3)
			   			{
 		     	     		min=Di3;
			     			minpoint=Index2;
			   			}
			 		}
                }
		 	}
	    }

      	Chrom->Membership[view][minpoint][Index1] = 1.0;



    }*/
     for (Index1 = 0; Index1 < n; Index1++)
	{
		min=9999999.99;
		flag=0;
		for (Index2 = 0; Index2 < MaxLen; Index2++)
		{
			Chrom->Membership[view][Index2][Index1]=0.0;
		    if(Chrom->z[Index2].x[0]!=DUMMY)
		    {
		    	i=(int)Chrom->z[Index2].x[0];
				j=(int)points[Index1][0];

				if(view==0 || view==2)
					Di3 =  distance_matrix[1][i][j]; //FindDistance(Chrom->z[Index2].x,points[Index1],view);
                if(view==1)
                	Di3 =  distance_matrix[0][i][j];
                if(!flag)
	            {
				    Chrom->Membership[view][Index2][Index1]=0.0;
			        if(Di3<=0.0)
					{
						Chrom->Membership[view][Index2][Index1]=1.0;
	                    flag=1;
				  		minpoint=Index2;
	                    break;

					}
			       	else
					{
				  		if(min>Di3)
				   		{
	 		     	   		min=Di3;
				    		minpoint=Index2;
				   		}
				 	}
	            }
			}
		}
		Chrom->Membership[view][minpoint][Index1] = 1.0;
	}	
    //free(Di2);
    Chrom->len=0;
    for (Index2 = 0; Index2 < MaxLen;Index2++)
    {
		printf("\nindex24 = %d\n",Index2);
			fflush(stdout);
		if(Chrom->z[Index2].x[0]!=DUMMY)
		{
		    Chrom->len++;
        	count=0;
         	for(Index1=0;Index1<n;Index1++){
            	if(Chrom->Membership[view][Index2][Index1]>0)
                  		count++;
            }
	   		Chrom->index1[Index2]=count;
        }
	}
	printf("\n++++chrom->len=%d+++++\n",Chrom->len);
	fflush(stdout);

		/*	printf("\n***********Membership for view %d\n**********************",1);
			fflush(stdout);
			for (int Index2 = 0; Index2 < MaxLen;Index2++)
		    {
				//if(Chrom->z[Index2].x[0]!=DUMMY && Chrom->z[Index2].x[0]>=0.0)
			//	{
				    for(int Index1=0;Index1<n;Index1++){
		            	printf("%lf ",Chrom->Membership[1][Index2][Index1]);

		            }
			   		printf("\n");
		      //  }
			}
			printf("\n***********XXXXXXXXX**********************\n");*/
//fflush(stdout);
}


void InitializePop()
{
	int i,j,m,l,dn,k,r,index,wrong,*g_pos, *Element;
	double delta;
	g_pos=(int*)calloc(MAX_LEN,sizeof(int));
	Element=(int*)calloc(NO_OF_ELM,sizeof(int));
	for(k=0;k<POP_SIZE;k++)
	{
	 	for(l=0;l<MaxLen;l++)
	 	{
	  		pool[k].index1[l]=0;
	  		for(dn=0;dn<d;dn++)
	    		pool[k].z[l].x[dn]=DUMMY;
	    }

	 	if (MaxLen == MinLen)
	 		r = MaxLen;
	 	else
	   		r=rand()%(MaxLen-MinLen)+MinLen;

	 	pool[k].len=r;
	    printf("\n len-%d\n",pool[k].len);
	 	r=rand()%MaxLen;
	 	g_pos[0]=r;
	//	printf("Running population number %d\n",k+1);
		for(i=1;i<pool[k].len;i++)
	  	{
	   		do{
	   				wrong=0;
	       		    r=rand()%MaxLen;
	       		    for(j=0;j<i;j++)
	         			if(g_pos[j]==r)
	         				wrong=1;
	     		}while(wrong==1);
	   		g_pos[i]=r;
	  	}
	 //	printf("\n In initialization \n");
	 	for(l=0;l<n;l++)
			Element[l] = 0;
	//	printf("Running 2");
	 	for(l=0;l<pool[k].len;l++)
	 	{
	    	do
	    	{
	       		r=rand()%n;
	    	}while (Element[r] == 1);
		    Element[r] = 1;
		    for(dn=0;dn<d;dn++)
		    {
				delta = frand();
				delta = 0.0;
				pool[k].z[g_pos[l]].x[dn]=points[r][dn];

			}
	 	}
	 //	printf("Running 3");

	}
	free(g_pos);
	free(Element);


}

int similarity_in_points()
{
	int i,j,total_similar=0,count=0,dn;
	for(i=0;i<n;i++)
	{
		for(j=(i+1);j<n;j++)
		{
			count=0;
			for(dn=0;dn<d;dn++)
			{
				if(points[i][dn]==points[j][dn])
					count++;
			}
			if(count==d)
				total_similar++;
		}
	}
	return(total_similar);
}

int see_similar(struct archive_element y)
{
	int i,j,count=0,flag=0,dn;
	for(i=0;i<MaxLen;i++)
	{

		for(j=(i+1);j<MaxLen;j++)
			{

				count=0;
				for(dn=0;dn<d;dn++)
				{
					if(y.z[i].x[dn]==y.z[j].x[dn])
						count++;
				}
				if(count==d)
				{
					flag=1;
					return(flag);
				}
			}


		}
	return(flag);
}

/* The following function will find the nondominated archive */
void form_nondominated_archive()
{
	int i,j,count1=0,count2=0,k,g,u,v,dn,ii,hill_climbno=2,f,count,jj;
	int *flag=(int *)malloc(POP_SIZE*sizeof(int));
	double Obj1;
	double **area1, *area2;
	//double **xnew;
	strcpy(func_name,"form_nondominated_archive");
	area1=(double **)malloc(POP_SIZE*sizeof(double));
	for(i=0;i<total_index;i++)
	{

		  range_min[i]=999999.99;
		  range_max[i]=0;
	}
	/*printf("\n******POOL from %s\n",func_name);
	for(i=0;i<POP_SIZE;i++)
	{
	    printf("\nPOOL Length= %d\n",pool[i].len);
		for(j=0;j<MaxLen;j++)
		{

			if(pool[i].z[j].x[0]!=DUMMY)
				for(k=0;k<d;k++)
					printf("%lf ",pool[i].z[j].x[k]);
		    printf("\n");
		}
		printf("\n");

	      // flag[i]=0;
	}
	getchar();*/

	for(i=0;i<POP_SIZE;i++)
	{
	    area1[i]=(double *)malloc(total_index*sizeof(double));
		for(j=0;j<total_index;j++)
		{
			if(no_view>1)
			{
				if (FitnessIndex[j] == 1) Obj1 = ComputeFitnessPBM(&pool[i],1);
				else if (FitnessIndex[j] == 2) Obj1 = ComputeFitnessPBM(&pool[i],2);
				else if (FitnessIndex[j] == 3) Obj1 = ComputeFitnessDAI(&pool[i]);
			}
			else
			{
				if (FitnessIndex[j] == 1) Obj1 = ComputeFitnessPBM(&pool[i],0);
				else if (FitnessIndex[j] == 2) Obj1 = ComputeFitnessPBM(&pool[i],0);

			}
			area1[i][j]=Obj1;

		}

	       flag[i]=0;
	}

//	printf("\n*********area1*******\n");
/*	for(i=0;i<POP_SIZE;i++)
	{
	   // area1[i]=(double *)malloc(total_index*sizeof(double));
		printf("\n");
		for(j=0;j<total_index;j++)
		{

			printf("%lf ",area1[i][j]);

		}


	}*/
	//getchar();
	area2=(double *)malloc(total_index*sizeof(double));
//	printf("\n IN Initialize ");

	for(ii=0;ii<POP_SIZE;ii++)
	{

	    for(jj=0;jj<hill_climbno;jj++)
	    {

	        for(i=0;i<MaxLen;i++)
	     	{
	            for(f=0;f<d;f++)
	            {
	                current.z[i].x[f]=pool[ii].z[i].x[f];

	            }
	     	}
	     	current.len=pool[ii].len;
	     	mutation();


	     /* 	for(i=0;i<MaxLen;i++)
	     	{
	     		if(new_pool.z[i].x[0]!=DUMMY)
	     		{
		     		printf("\n****new pool***** %d\n",new_pool.len);
		            for(f=0;f<d;f++)
		            {

		                	printf("%lf ",new_pool.z[i].x[f] );
		            }
	            	printf("\n");
	            }
	     	}
	     	//getchar();*/
		    if(no_view>1)
		    {
		      	for(i=0;i<30;i++)
				{
					track_count=i;
printf("\n003\n");
					ComputeMembership(&new_pool,1);
					UpdateCenter(&new_pool,1);
					printf("\n3\n");

				}
				for(i=0;i<30;i++)
				{
					track_count=i;
					ComputeMembership(&new_pool,2);
printf("\n003\n");
					UpdateCenter(&new_pool,2);
printf("\n03\n");

				}
				CombineCenter(&new_pool);
			}
			else
			{
				for(i=0;i<30;i++)
				{
					track_count=i;
					ComputeMembership(&new_pool,0);
					UpdateCenter(&new_pool,0);

				}
			}
			//perror("Error in ii: ");
	           printf("\nout mutation\n");
		    for(j=0;j<total_index;j++)
		    {
				if(no_view>1)
				{
			        if (FitnessIndex[j] == 1) Obj1 = ComputeFitnessPBM(&new_pool,1);
					else if (FitnessIndex[j] == 2) Obj1 = ComputeFitnessPBM(&new_pool,2);
					else if (FitnessIndex[j] == 3) Obj1 = ComputeFitnessDAI(&new_pool);
				}
				else
				{
					if (FitnessIndex[j] == 1) Obj1 = ComputeFitnessXB(&new_pool,0);
					else if (FitnessIndex[j] == 2) Obj1 = ComputeFitnessPBM(&new_pool,0);
				}

		        area2[j]=Obj1;
			    printf("area2[%d] %lf ",j,area2[j]);
		    }

		    count=0;
		    for(i=0;i<total_index;i++)
		    {
		        if(area1[ii][i]>=area2[i])
		            count++;
		    }
		    if(count==total_index)
		    {
		        for(i=0;i<MaxLen;i++)
			    {
		            for(f=0;f<d;f++)
		            {
		                pool[ii].z[i].x[f]=new_pool.z[i].x[f];
		            }
			    }

				for(j=0;j<total_index;j++)
			    {
				  	area1[ii][j]=area2[j];
				}

				for(u=0;u<MaxLen;u++)
			    {
			        for(v=0;v<n;v++)
				        pool[ii].Membership[0][u][v]=new_pool.Membership[0][u][v];
			    }
			}

		}


	}
	 printf("\n At the end of initialize solution\n\n");
	for(i=0;i<POP_SIZE;i++)
	{
		if(flag[i]==0)
			{
			for(j=i+1;j<POP_SIZE;j++)
				{
				if(flag[i]==0)
				{
				if(flag[j]==0)
				{
					count1=0;count2=0;
					for(k=0;k<total_index;k++)
					{
						if(area1[i][k]>=area1[j][k])
							count1++;
						if(area1[i][k]<=area1[j][k])
							count2++;
					}
					if(count1==total_index)
						flag[i]=1;
					else if(count2==total_index)
						flag[j]=1;

				}
				}
			}
		}
	}

	area=(double **)malloc((softl)*sizeof(double));

	for(j=0;j<(softl);j++)

	{
		area[j]=(double *)malloc(total_index*sizeof(double));

	}
	g=0;
	// fprintf(fpo,"\n POP_SIZE=%d",POP_SIZE);

	for(i=0;i<POP_SIZE;i++)

	{
	// fprintf(fpo,"\n flag[%d]=%d\n",i,flag[i]);

		if(flag[i]==0)
		{
		//archive[g].len=pool[i].len;
		// fprintf(fpo,"\n length=%d",archive[g].len);

		    for(k=0;k<MaxLen;k++)
		    {
		       for(dn=0;dn<d;dn++)
		           {
			    archive[g].z[k].x[dn]=pool[i].z[k].x[dn];
			     //printf("\n  %lf",archive[g].z[k].x[dn]);
			   }
		    }

		    for(k=0;k<=no_view;k++)
		    {
			   for(u=0;u<MaxLen;u++)
			    {
			        for(v=0;v<n;v++)
				    	archive[g].Membership[k][u][v]=pool[i].Membership[k][u][v] ;
			    }
			}

		   for(j=0;j<total_index;j++)

		     {
			 area[g][j]=area1[i][j];
			 if(range_max[j]< area[g][j])
			   range_max[j]=area[g][j];
		          if(range_min[j]> area[g][j])
			   range_min[j]=area[g][j];
	//	        printf("\n area[%d][%d]=%lf",g,j,area[g][j]);
		     }

		 g++;
	  }
	}
//	printf("\n Size of g is:%d",g);
	arcsize=g;
	free(flag);
	free(area2);
	for(i=0;i<POP_SIZE;i++)
	{
	    free(area1[i]);
	}
	free(area1);


}

void adjacency_matrix( struct archive_element *c,int view,int **adj)
{
	int i=0,j=0,k=0;
	for(i=0;i<MaxLen;i++)
	{
		for(j=0;j<n;j++)
		{
			for(k=0;k<n;k++)
			{
				//if(k==j)
				//+	continue;
				if(c->Membership[view][i][k]>0.0 && c->Membership[view][i][j]>0.0){

					adj[j][k]=1;

				}
			}
		}
	}

}


double ComputeFitnessDAI(struct archive_element *c)
{
	int ***adj,i=0,j=0,k=0,l=0,m=0;
	double DAI=0.0,agree=0.0,disagree=0.0;

	adj= (int***)calloc((no_view+1),sizeof(int**));
	//printf("\nIn DAI1\n");
	//getchar();
	for(i=0;i<=no_view;i++)
	{
		adj[i]= (int**)calloc(n,sizeof(int*));
		for(j=0;j<n;j++)
			adj[i][j]= (int*)calloc(n,sizeof(int));
	}
	//printf("\nIn DAI2\n");
	//getchar();
	for(i=1;i<=no_view;i++)
		adjacency_matrix(c,i,adj[i]);

	//printf("\nIn DAI3\n");
	for(i=1;i<=no_view;i++)
	{
		for(j=1;j<=no_view;j++)
		{
			if(i==j)
				continue;
			for(k=0;k<n;k++)
			{
				for(l=0;l<n;l++)
					adj[0][k][l]= (adj[i][k][l] && adj[j][k][l]); // works for 2 view problem
			}

		}
	}

	/*for(i=0;i<=no_view;i++)
	{
		printf("\nadjacency matrix of view %d********\n",i);
		for(k=0;k<n;k++){
			for(j=0;j<n;j++)
			 printf("%d ",adj[i][k][j]);
			printf("\n");
		}
		printf("\n");
		getchar();
	}*/

	for(k=0;k<n;k++)
	{
		for(l=0;l<n;l++)
			if(adj[0][k][l]>0)
				agree++;
			else
				disagree++;
	}
	printf("\nagree=%lf disagree=%lf",agree,disagree);
	//free memory-----------------------
	for(i=0;i<=no_view;i++)
	{
		for(j=0;j<n;j++)
			free(adj[i][j]);
		free(adj[i]);
	}
	free(adj);
	//----------------------------------
	if(agree<=0.0)
		return 100000;
	DAI= disagree/agree;
	printf("\nDAI=%lf\n",DAI);


	return DAI;
}

double ComputeFitnessPBM(struct archive_element *c,int view)
{
	int i,j,k,dn,len,cl,m,p;
	//printf("\nIN pBM\n");
	//getchar();
	double ed,ed1,inter,max,min=9999999.999,fitns,product;
	double s2,sum2,weight=2.0;
	ed=0.0;
	len=0;

	for(i=0;i<MaxLen;i++)
	{
		if (c->z[i].x[0] != DUMMY)
		{ //ask

			len++;
			for(j=0;j<n;j++)
			{

					s2=0.0;
					for (dn=0;dn<d;dn++)
					{
						/*if(views[dn]==view)
							s2+=pow((c->z[i].x[dn]-points[j][dn]),2.0);//what points does: two dimensional
						if(view==0)
		                    s2+=pow((c->z[i].x[dn]-points[j][dn]),2.0);      */
		                m= (int)c->z[i].x[dn];
		            	p= (int)points[j][dn];

		                if(view==0 || view==2)
		                	s2= distance_matrix[1][m][p];
		                else
		                	s2= distance_matrix[0][m][p];
					}


					ed=ed + s2*pow(c->Membership[view][i][j],weight);

			}
	    }
	}


	inter=0.0;
	for(i=0;i<MaxLen;i++)
	{
		if (c->z[i].x[0] != DUMMY)
		{
			for(k=0;k<MaxLen;k++)
			{
				 if (c->z[k].x[0] != DUMMY)
				 {
					if (i!= k)
					{    /* valid center */
						sum2=0.0;
						for(dn=0;dn<d;dn++)
						{
							/*if(views[dn]==view)
								product=pow((c->z[k].x[dn]-c->z[i].x[dn]),2.0);
							else
							{
								if(view==0)
								{
									product=pow((c->z[k].x[dn]-c->z[i].x[dn]),2.0);
								}
								else
									product=0;

							}*/
								product =0.0;
							m= (int)c->z[k].x[dn];
							p= (int)c->z[i].x[dn];
							if(view==0 || view==2)
			                	product= distance_matrix[1][m][p];
			                if(view==1)
			                	product= distance_matrix[0][m][p];

							sum2=sum2+product;
						}
						if(sum2>inter)
							inter=sum2;
					}
				}
			}
		}
	}

	if(len>=2 && inter > 0.0)
		fitns=(inter)/(len*ed);

	else
		fitns=0.000001;
	printf("\nPBM fitness= %lf\n",fitns);
	//getchar();
	if(fitns==0.0)
		fitns=0.000001;
	return(1.0/fitns);
}

double ComputeFitnessXB(struct archive_element *c,int view)
{
	int i,j,k,dn;
	double s2,sum2,sum=0.0,min=9999999.99,dd,temp,nsum;
	double XB=0.0;
	c->len=0;
	for(i=0;i<MaxLen;i++)
	{
		if (c->z[i].x[0] != DUMMY)
		{
			c->len++;
			for(k=0;k<MaxLen;k++)
			{
			 	if (c->z[k].x[0] != DUMMY)
			 	{
						if (i!= k)
						{
							sum2=0.0;
							for (dn=0;dn<d;dn++)
							{
								if(view==views[dn])
									dd=pow((c->z[k].x[dn]-c->z[i].x[dn]),2.0);
								else
								{
									if(view==0)
										dd=pow((c->z[k].x[dn]-c->z[i].x[dn]),2.0);
									else
										dd=0.0;
								}
								sum2=sum2+dd;
							}
							if(sum2<min)
								min=sum2;
	    				}
	   			}
			}
		}
	}
	//printf("\n len=%d",c->len);
	//printf("\n minimum separation=%lf",sum2);
	if (min == 0)/*two clusters centers same*/
	{
		XB=9999999999.999;
		 //printf("\n%lf\n",XB);
		return(XB);
	}
	sum=0.0;nsum=0.0;
	for(i=0;i<MaxLen;i++)
	{
		if (c->z[i].x[0] != DUMMY)
		{
			for(j=0;j<n;j++)
			{
				s2=0.0;
				for (dn=0;dn<d;dn++)
				{
					if(view==views[dn])
						s2+=pow((c->z[i].x[dn]-points[j][dn]),2.0);
					if(view==0)
						s2+=pow((c->z[i].x[dn]-points[j][dn]),2.0);
				}
				sum += c->Membership[view][i][j] * c->Membership[view][i][j] * s2;
			}
	  	}
	 }
	if(c->len>=2)
	{
		XB = sum/(n * min);
		if (XB == 0)
			XB = 0.0000001;

	}
	else
		XB=99999.99;
	//printf("\nXB=%lf\n",XB);
	return (XB);
}


double FindDistance(double *x, double *y,int view)
{
	double distance=0.0,sum=0.0;
	int i,flag=0;
	for(i=0;i<d;i++)
       {
			if(view==views[i]) // considering only those dimensions belonging to that particular view
				sum=sum+pow((x[i]-y[i]),2);
			if(view==0)
				sum=sum+pow((x[i]-y[i]),2);
		}

	distance=sqrt(sum);

	return(distance);
}

void writechromosome(struct archive_element y)
{
	/*int i,j;
	printf("\n\n");
	for(i=0;i<MaxLen;i++)
	{
		printf("(");
		for(j=0;j<d;j++)
			{
			printf("%lf",y.z[i].x[j]);
			}
		printf(")");
	}*/
}

void printclustering(FILE *c_ptr,struct archive_element *y)//important
{
  int i,j,cl,*clas;
  float max;
  clas=(int*)malloc(n*sizeof(int));
  for (j=0;j<n;j++)
    {
    	max = 0.0;
	  	for (i=0;i<MaxLen;i++)
		{
			if (y->Membership[0][i][j]>max)
			{
		    	max = y->Membership[0][i][j];
		    	cl=i;
		  	}
		}
    	for (i=0;i<d;i++)
			fprintf(c_ptr,"%f   ",points[j][i]);
    	fprintf(c_ptr,"%d\n",cl+1);
    	clas[j]=cl+1;
    }
    free(clas);
	

}
void clustering()
{
    int i,j,k,l,h,p,g2,g,pp,dn,uu,vv,jj,kk,hh;
    double **area1;
	struct archive_element * archive1;
	int *point1, *point2;
	int u=0,v=0,w,m;
	double *dist,min, **distance;
    int no_clus=arcsize;
	int *flag;
    int *cluster=(int *)malloc(sizeof(int)*(softl+1));
    double **arc2;

    archive1=(struct archive_element *)malloc((softl+2)*sizeof(struct archive_element));
	area1=(double **)malloc((softl+2)*sizeof(double *));
	for(i=0;i<(softl+2);i++)
		area1[i]=(double *)malloc(total_index*sizeof(double));
	for(h=0;h<arcsize;h++)
	{
		if(no_view>1)
			archive1[h].Membership=(double ***)malloc((no_view+1)*sizeof(double **));
		else
			archive1[h].Membership=(double ***)malloc(sizeof(double **));
		for(j=0;j<=no_view;j++)
		{
			archive1[h].Membership[j]=(double **)malloc(MaxLen*sizeof(double *));
			for(i=0;i<MaxLen;i++)
			    archive1[h].Membership[j][i]=(double *)malloc(n*sizeof(double));
		}

	 }

	fflush(stdout);
    point1=(int *)malloc(sizeof(int)*(softl+1));
    point2=(int *)malloc(sizeof(int)*(softl+1));
    dist=(double *)malloc(sizeof(double)*softl);

    distance=(double **)malloc((softl+1)*sizeof(double *));
	k=arcsize;
    for(i=0;i<(softl+1);i++)
        distance[i]=(double *)malloc((softl+1)*sizeof(double));
    for(i=0;i<arcsize;i++)
        cluster[i]=i;
    for(i=0;i<k;i++)
    {
        distance[i][i]=2000000;
        for(j=i+1;j<k;j++)
        {
		    distance[i][j]=0.0;
            for(p=0;p<total_index;p++)
		        distance[i][j]=distance[i][j]+pow((area[i][p]-area[j][p]),2);
		    distance[j][i]=sqrt(distance[i][j]);
        }

    }
  printf("Whats happen 10\n");
	fflush(stdout);

    flag=(int *)malloc((softl+1)*sizeof(int));
    while(no_clus>hardl)
    {
        min=2000000;
        for(i=0;i<arcsize;i++)
            flag[i]=0;
        for(i=0;i<k;i++)
        {
            for(j=0;j<k;j++)
            {
                if(j!=k)
                {
                    if(min>distance[i][j])
                    {
                        min=distance[i][j];
                        u=i;
                        v=j;
                    }
                }
            }
        }
                if(cluster[u]==u && cluster[v]==v)
                   {
                     cluster[u]=v;
                     cluster[v]=u;
                    }
                 else if(cluster[u]==u)
                    {
                      j=cluster[v];
                      while(cluster[j]!=v)
                          {
                            j=cluster[j];
                           }
                       cluster[j]=u;
                       cluster[u]=v;
                     }
                  else if(cluster[v]==v)
                     {
                        j=cluster[u];
                        while(cluster[j]!=u)
                          {
                            j=cluster[j];
                          }
                         cluster[j]=v;
                         cluster[v]=u;
                       }
                   else
                      {
                         j=cluster[u];
                         while(cluster[j]!=u)
                            {
                              j=cluster[j];
                             }
                          cluster[j]=v;
                          p=cluster[v];
                          while(cluster[p]!=v)
                             {
                               p=cluster[p];
                              }
                           cluster[p]=u;
                      }
                    no_clus=no_clus-1;
                    g=0;
                    point1[g]=u;
                    j=cluster[u];
                    while(j!=u)
                        {
                          g++;
                          point1[g]=j;
                          j=cluster[j];
                         }
                     for(i=0;i<=g;i++)
                         {
                            w=point1[i];
                            flag[w]=1;
                            for(j=i+1;j<=g;j++)
                              {

                                m=point1[j];
                                flag[m]=1;
                                distance[m][w]= distance[w][m]=2000000;

                               }
                          }
                       for(i=0;i<arcsize;i++)
                           {
                             if(flag[i]==0)
                              {  /*see its's end bracket*/
                               if(cluster[i]==i)
                                  {
                                     w=point1[0];
                                     min=distance[w][i];
                                     for(j=1;j<=g;j++)
                                        {
					   m=point1[j];
                                           if(min>distance[m][i])
                                                min=distance[m][i];
                                           }
                                       for(j=0;j<=g;j++)
                                          {
                                            m=point1[j];
                                            distance[m][i]=min;
                                           }
                                        flag[i]=1;
                                     }
                                  else
                                     {
                                         g2=0;
                                         point2[g2]=i;j=cluster[i];
                                         while(j!=i)
                                             {
                                               g2++;
                                               point2[g2]=j;
                                               j=cluster[j];
                                              }
                                          w=point1[0];
                                          m=point2[0];
                                          min=distance[w][m];
                                          for(j=0;j<=g;j++)
                                            {
                                              w=point1[j];
                                              for(p=0;p<=g2;p++)
                                               {
                                                m=point2[p];
                                                if(min>distance[w][m])
                                                    min=distance[w][m];
                                                }
                                              }
                                          for(j=0;j<=g;j++)
                                               {
                                                 for(p=0;p<=g2;p++)
                                                      {
                                                        w=point1[j];
                                                        m=point2[p];
                                                        distance[m][w]=distance[w][m]=min;
                                                        flag[m]=1;
                                                       }
                                                    }

                                          }
                                      }
                                 }
                             }
			 for(hh=0;hh<arcsize;hh++)
	                    {

		                for(jj=0;jj<MaxLen;jj++)
		                 {
		                    for(dn=0;dn<d;dn++)
		                          archive1[hh].z[jj].x[dn]=archive[hh].z[jj].x[dn];
		                  }
		               for(kk=0;kk<=no_view;kk++)
		                 for(uu=0;uu<MaxLen;uu++)
		                   {
		                     for(vv=0;vv<n;vv++)
			                  archive1[hh].Membership[kk][uu][vv]=archive[hh].Membership[kk][uu][vv];
		                   }
		                  for(uu=0;uu<total_index;uu++)
		                     {
		                      area1[hh][uu]=area[hh][uu];
		                    }

	                       }
	printf("Whats happen 11\n");
fflush(stdout);

                  for(i=0;i<arcsize;i++)
                    {
                     flag[i]=0;
                     }
                  k=0;
                  for(i=0;i<arcsize;i++)
                    {
		     if(flag[i]==0)
		     {
                     if(cluster[i]!=i)
                       {
                        g=0;
                        point1[g]=i;
                        flag[i]=1;
                        j=cluster[i];
                        while(j!=i)
                          {
                           g++;
                           point1[g]=j;
                           flag[j]=1;
                           j=cluster[j];

                          }

                        for(j=0;j<=g;j++)
                          {
                           dist[j]=0;
                           w=point1[j];
                           for(p=0;p<=g;p++)
                              {
			       if(p!=j)
			       {
                                  m=point1[p];
			         for(pp=0;pp<total_index;pp++)
                                   dist[j]=dist[j]+pow((area[w][pp]-area[m][pp]),2);
				 dist[j]=sqrt(dist[j]);
                              }
			     }
                           }

                        min=dist[0];
                        w=point1[0];
                        for(j=1;j<=g;j++)
                          {
                           if(min>dist[j])
                            {
                             min=dist[j];
                             w=point1[j];
                            }
                          }

		for(jj=0;jj<MaxLen;jj++)
		{
		   for(dn=0;dn<d;dn++)
		     archive[k].z[jj].x[dn]=archive1[w].z[jj].x[dn];
		}

		for(uu=0;uu<MaxLen;uu++)
		{
		   for(vv=0;vv<n;vv++)
			archive[k].Membership[0][uu][vv]=archive1[w].Membership[0][uu][vv];
		}
		for(uu=0;uu<total_index;uu++)
		{
		    area[k][uu]=area1[w][uu];
		}
		k++;

            }
           else
               {
	            for(jj=0;jj<MaxLen;jj++)
		       {
		           for(dn=0;dn<d;dn++)
		          archive[k].z[jj].x[dn]=archive1[i].z[jj].x[dn];
		      }

		     for(uu=0;uu<MaxLen;uu++)
		       {
		           for(vv=0;vv<n;vv++)
			      archive[k].Membership[0][uu][vv]=archive1[i].Membership[0][uu][vv];
		       }
		     for(uu=0;uu<total_index;uu++)
		        {
		           area[k][uu]=area1[i][uu];
		        }

                     k++;
                  }
                }
	      }
              arcsize=k;
  //            printf("\n arcsize=%d",arcsize);
//	      printf("\n afterclustering:");
	      /*for(i=0;i<arcsize;i++)
              {
                printf("\n %f",area[i][0]);
                printf("\t %f",area[i][1]);
              }*/
          }
double minkwski(char *name)
{
	int i,j,k,sol[1000],st,count1=0,count2=0,count3=0; //changes depending on size of datapoints
	double min, mink,as,bs;
	char ch,name1[30];
	FILE *ifpa;
	if((ifpa=fopen(name,"r"))==NULL)
	{
	printf("\n File not found");
	exit(9);
	}
	fscanf(ifpa,"%d %d %d", &i,&j,&k);
	for(i=0;i<n;i++)
	{
           for(j=0;j<d;j++)
            fscanf(ifpa,"%lf",&as);

		fscanf(ifpa,"%d",&sol[i]);
	//	printf("%d\t",sol[i]);
	}
   fclose(ifpa);
	for(i=0;i<(n-1);i++)
	{
		for(j=i+1;j<n;j++)
		{
			//if(i!=j)
			//{
			if(fclass[i]==fclass[j] && sol[i]==sol[j])
			count1++;
			if(fclass[i]==fclass[j] && sol[i]!=sol[j])
			count2++;
			if(fclass[i]!=fclass[j] && sol[i]==sol[j]){
		//	printf("It happens\n");
		//	printf("tru[%d]=%d, tru[%d]=%d, sol[%d]=%d, sol[%d]=%d \n",i,tru[i],j,tru[j],i,sol[i],j,sol[j]);
			count3++;
                    }
			//}
		}
	}
	//printf("count1=%d, count2=%d, count3=%d\n",count1,count2,count3);
	as=count2+count3;
	bs=count1+count2;
	min=(as/bs);

	mink=sqrt(min);
	printf("Minkoski score is = %.4lf",mink);
	return(mink);
}

