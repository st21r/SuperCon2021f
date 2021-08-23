# include <stdio.h>
# include <math.h>
# include "cssl.h"
# include <stdlib.h>
# include <omp.h>
# include "mpi.h"
# include "sc21.h"

# define maxstep_greedy 100000
# define maxreplica 48
# define seed0 12345

void const_network_random( int N_LINK, int ix, int id);

double S[N_GROUP];
double I[N_GROUP];
double R[N_GROUP];

double S0[N_GROUP];
double I0[N_GROUP];
double R0[N_GROUP];

void sim(int n_step);

void swap_C_random( int *i1, int *j1, int *i2, int *j2, double r[4], int *iswap);

void swap_C_targets( int i1, int j1, int i2, int j2);

double l2_norm( double p1[N_GROUP], double p2[N_GROUP]);

int main(int argc, char **argv){

  int seed[maxreplica];

  /*  --   for MPI  ----*/ 

  int myid,n_procs;

  MPI_Status status;

  MPI_Init( &argc, &argv );
  SC_input();

  MPI_Comm_size(MPI_COMM_WORLD, &n_procs); 
  MPI_Comm_rank(MPI_COMM_WORLD, &myid); 

  if (myid == 0){

    int n;

    int ix,icon;
    double dwork[8];
    double r[maxreplica];
  
    ix=seed0;
    c_dm_vranu5(&ix, r, n_procs, 0, dwork, &icon);
    for (n=0; n < n_procs; n++){
      seed[n]=abs((int)(r[n]*2147483647));
      /*printf("#seed[%d]=%d\n",n,seed[n]);*/
    };
  };

  MPI_Bcast(&seed,n_procs,MPI_INT,0,MPI_COMM_WORLD);

  /* inference simulation*/

  const_network_random(N_LINK,seed[myid],myid);

  int i;
  for (i=0;i< N_GROUP;i++){
    S0[i]=(double) N[i];
    I0[i]=0.0;
    R0[i]=0.0;
  };

  I0[0]=1.0;
  S0[0]=(double) (N[0]-1);
  R0[0]=0.0;

  sim(T);

  double l20,l2,dl2;

  l20=l2_norm(I_PROB,I);

  int istep_greedy;
  int i1,i2,j1,j2;

  int ix,icon;
  double dwork[8];
  double r[4];

  int iswap;

  ix=seed[myid];

  istep_greedy=0;
  while (istep_greedy < maxstep_greedy){

    iswap=0;
    while (iswap == 0){
      c_dm_vranu5(&ix, r, 4,  (long) 0, dwork, &icon);
      swap_C_random(&i1,&j1,&i2,&j2,r,&iswap);
    };

    for (i=0;i< N_GROUP;i++){
      S0[i]=(double) N[i];
      I0[i]=0.0;
      R0[i]=0.0;
    };
    I0[0]=1;
    S0[0]=(double) (N[0]-1);
    R0[0]=0.0;

    sim(T);

    l2=l2_norm(I_PROB,I);
    dl2=l2-l20;

    if (dl2 > (double) 0.0){

      swap_C_targets(i1,j1,i2,j2);

      l2=l20;

    }else{
      l20=l2;     
    };

    /*          printf("%d %d %lf \n",myid,istep_greedy,l2);*/

    istep_greedy++;

  };

  /*
    int M_LINK;
    int j;
    M_LINK=0;
    for (i=0;i< N_GROUP;i++){
    for (j=i;j< N_GROUP;j++){
    if (C[i][j]==1){
    M_LINK=M_LINK+1;
    };
    };
    };
    printf ("#myid,M_LINK,N_LINK=%d %d %d \n",myid,M_LINK,N_LINK);
  */

  double *l2box;
  l2box=(double *)malloc(n_procs*sizeof(double));

  MPI_Gather(&l2, 1, MPI_DOUBLE, l2box, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int n_min;

  if (myid==0){

    int n;
    double l2min;

    l2min=l2box[0];
    n_min=0;
    for (n=0; n < n_procs; n++){
      /*   printf("#l2[%d]=%lf\n",n,l2box[n]);*/
      if (l2min > l2box[n]){
	l2min=l2box[n];
	n_min=n;
      };
    };
    
    /*  printf("#l2min,n_min= %lf %d\n",l2min,n_min);*/

  };

  /*  for (i=0; i < N_GROUP; i++){printf("%d %d %f %f\n",myid,i,I[i],I_prob[i]);};*/

  MPI_Bcast(&n_min,1,MPI_INT,0,MPI_COMM_WORLD);

  int itag=0;
  if (n_min != 0){
    if (myid == n_min){
      MPI_Send(C,(N_GROUP*N_GROUP),MPI_INT,0,itag,MPI_COMM_WORLD);
    }else{
      if (myid ==0){
	MPI_Recv(C,(N_GROUP*N_GROUP),MPI_INT,n_min,itag,MPI_COMM_WORLD, &status);
      };
    };
  };

  if (myid == 0){ 
    SC_output();
  };

  MPI_Finalize();

  return 0;

};


void const_network_random( int N_LINK, int ix, int id){

  int i,j;

  for (i=0;i< N_GROUP;i++){
    for (j=0;j< N_GROUP;j++){
      C[i][j]=0;
    };
  };

  int m;

  int icon;
  double dwork[8];

  double r[2];

  m=0;
  while (m < N_LINK){

    c_dm_vranu5(&ix, r, 2,   (long)   0, dwork, &icon);

    i=(int) (r[0]*N_GROUP);
    j=(int) (r[1]*N_GROUP);

    if(C[i][j]==0 && i != j){
      m=m+1;
      C[i][j]=1;
      C[j][i]=1;

    };
  };

};

void sim(int n_step){

  int i,j,istep;

  double contact;

#pragma omp parallel
  {

    for (istep=0; istep< n_step; istep++){

#pragma omp for private(j,contact)
      for (i=0; i < N_GROUP; i++){ 

	contact=BETA*I0[i];
	for (j=0; j < N_GROUP; j++){ 
	  contact=contact+BETA2*C[i][j]*I0[j];
	};

	S[i]=S0[i]-contact*S0[i];
	I[i]=I0[i]+contact*S0[i]-GAMMA*I0[i];
	R[i]=R0[i]+GAMMA*I0[i];

      };

#pragma omp for
      for (i=0; i < N_GROUP; i++){ 
	S0[i]=S[i];
	I0[i]=I[i];
	R0[i]=R[i];
      };

    };
  };
};

    

double l2_norm( double p1[N_GROUP], double p2[N_GROUP]){

  int i;
  double l2;

  l2=(double) 0.0;
  for (i=0; i< N_GROUP; i++){
    l2=l2+(p1[i]-p2[i])*(p1[i]-p2[i]);
  };
  return(l2);
  
};


void swap_C_random(int *i1_out, int *j1_out, int *i2_out, int *j2_out,  double r[4], int *iswap){

  int i1,i2,j1,j2;

  *iswap=0;

  i1=(int) (r[0]*N_GROUP);
  j1=(int) (r[1]*N_GROUP);
  i2=(int) (r[2]*N_GROUP);
  j2=(int) (r[3]*N_GROUP);


  int ic;
  if (i1 !=j1 && i2 != j2 && C[i1][j1]+C[i2][j2]==1){

    ic=C[i1][j1];
    C[i1][j1]=C[i2][j2];
    C[i2][j2]=ic;

    C[j1][i1]=C[i1][j1];
    C[j2][i2]=C[i2][j2];

    *iswap=1;
  };

  *i1_out=i1;
  *j1_out=j1;
  *i2_out=i2;
  *j2_out=j2;

};

void swap_C_targets( int i1, int j1, int i2, int j2){

  int ic;

  ic=C[i1][j1];
  C[i1][j1]=C[i2][j2];
  C[j1][i1]=C[j2][i2];
  C[i2][j2]=ic;
  C[j2][i2]=ic;

};

