#include<math.h>
double FindDistance(double *x, double *y,int dim)
{
	double distance=0.0,sum=0.0;
	int i,flag=0;
	for(i=0;i<dim;i++)
       sum=sum+pow((x[i]-y[i]),2);
	distance=sqrt(sum);
	return(distance);
}