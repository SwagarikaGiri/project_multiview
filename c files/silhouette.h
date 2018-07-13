#ifndef SILHOUETTE_H_   /* Include guard */
#define SILHOUETTE_H_

double ComputeFitnessSilhouette(double **Membership,double **points,int no_cluster,int no_points,int dim_points);  /* Function declaration */
double ComputeFitnessSilhouette(int **Membership,double **points,int no_cluster,int no_points,int dim_points);

#endif // SILHOUETTE_H_