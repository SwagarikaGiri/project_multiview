#include<limits.h>
#include "distance.h"
#define MAX(x,y) ((x) > (y) ? (x) : (y))

double ComputeFitnessSilhouette(double **Membership,double **points,int no_cluster,int no_points,int dim_points)
{
	int i,j,k,l,count_inclus=0,count_outclus=0;
	double min=(double)LONG_MAX,sum=0.0,b_i=0.0,a_i=0.0,avg=0.0, s_i=0.0, s_k=0.0, SI=0.0;
	
	for(i=0;i<no_cluster;i++)
	{
		count_inclus=0;
		s_k=0.0;
		for(j=0;j<no_points;j++) //counting total number of points in the cluster C_i
			if(Membership[i][j]==1.0)
				count_inclus++;
		if(count_inclus>0)
		{
		
			for(j=0;j<no_points;j++) 
			{
				sum=0.0;
				min=(double)LONG_MAX;
				if(Membership[i][j]==1.0)
				{	
					// Measuring intra cluster distance
					
					for(k=0;k<no_points;k++) // Measuring intra cluster distance
					{
						if(Membership[i][k]==1.0)
						{
							sum+=FindDistance(points[j],points[k],dim_points);
						}
						
					}
					if(count_inclus>1)
						a_i = sum/((double)(count_inclus-1));
					else
						a_i = 0.0;
					
					//----------------------------------------------------

					// Measuring inter cluster distance
					
					for(k=0;k<no_cluster;k++) // Measuring inter cluster distance
					{
						count_outclus=0;
						sum=0.0;
						if(i!=k)
						{
							for(l=0;l<no_points;l++) //counting total number of points in the cluster C_k where C_k != C_i
							{
								if(Membership[k][l]==1.0)
								{
									count_outclus++;
								}
							}
							for(l=0;l<no_points;l++) //counting total number of points in the cluster C_k where C_k != C_i
							{
								if(Membership[k][l]==1.0)
								{
									sum+= FindDistance(points[l],points[j],dim_points);
								}
							}
							if(count_outclus > 0)
								b_i = sum/(double)count_outclus;
							else
								b_i = (double)LONG_MAX;

							if(min > b_i)
								min=b_i;

						}

					}
					b_i = min; // minimum cluster distance is stored in b_i
					
					if(b_i >= (double)LONG_MAX)
						b_i = 0.0;

					//---------------------------------------------------------------

					s_i = (b_i - a_i)/(MAX(a_i,b_i));
					s_k += s_i;
				}
			}
			s_k = s_k/(double)count_inclus;

			SI += s_k;
		}
	}

	SI = SI/(double)no_cluster;

	return SI;
}

double ComputeFitnessSilhouette(int **Membership,double **points,int no_cluster,int no_points,int dim_points)
{
	int i,j,k,l,count_inclus=0,count_outclus=0;
	double min=(double)LONG_MAX,sum=0.0,b_i=0.0,a_i=0.0,avg=0.0, s_i=0.0, s_k=0.0, SI=0.0;
	
	for(i=0;i<no_cluster;i++)
	{
		count_inclus=0;
		s_k=0.0;
		for(j=0;j<no_points;j++) //counting total number of points in the cluster C_i
			if(Membership[i][j]==1)
				count_inclus++;
		if(count_inclus>0)
		{
		
			for(j=0;j<no_points;j++) 
			{
				sum=0.0;
				min=(double)LONG_MAX;
				if(Membership[i][j]==1)
				{	
					// Measuring intra cluster distance
					
					for(k=0;k<no_points;k++) // Measuring intra cluster distance
					{
						if(Membership[i][k]==1)
						{
							sum+=FindDistance(points[j],points[k],dim_points);
						}
						
					}
				
					if(count_inclus>1)
						a_i = sum/((double)(count_inclus-1));
					else
						a_i = 0.0;

									
					//----------------------------------------------------

					// Measuring inter cluster distance
					
					for(k=0;k<no_cluster;k++) // Measuring inter cluster distance
					{
						count_outclus=0;
						sum=0.0;
						if(i!=k)
						{
							for(l=0;l<no_points;l++) //counting total number of points in the cluster C_k where C_k != C_i
							{
								if(Membership[k][l]==1)
								{
									count_outclus++;
								}
							}
							for(l=0;l<no_points;l++) //counting total number of points in the cluster C_k where C_k != C_i
							{
								if(Membership[k][l]==1)
								{
									sum+= FindDistance(points[l],points[j],dim_points);
								}
							}
							
							if(count_outclus > 0)
								b_i = sum/(double)count_outclus;
							else
								b_i = (double)LONG_MAX;

							if(min > b_i)
								min=b_i;

						}

					}
					b_i = min; // minimum cluster distance is stored in b_i

					if(b_i >= (double)LONG_MAX)
						b_i = 0.0;
					
					//---------------------------------------------------------------

					s_i = (b_i - a_i)/(MAX(a_i,b_i));
				
					s_k += s_i;

				}
			}
			
			s_k = s_k/(double)count_inclus;
			
			SI += s_k;
		}
	}

	SI = SI/(double)no_cluster;

	return SI;
	
}