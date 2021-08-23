# define N_GROUP 100

const double BETA=0.0002;
const double BETA2=0.000001;
const double GAMMA=0.1;

const int T=200;

void SC_input();
void SC_output();

int N_LINK;

int N[N_GROUP];
double I_PROB[N_GROUP];

int C[N_GROUP][N_GROUP];

double TIME0;

void SC_input(){
    int rank = 0;
#ifdef MPI_VERSION
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if(0 == rank){
        scanf("%d",&N_LINK);
        int i;
        for (i=0; i< N_GROUP; i++){
            scanf("%d\n",&N[i]);
        }
        for (i=0; i< N_GROUP; i++){
            scanf("%lf",&I_PROB[i]);
        }
    }
#ifdef MPI_VERSION
    MPI_Bcast(&N_LINK, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(N, N_GROUP, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(I_PROB, N_GROUP, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    TIME0=omp_get_wtime();
}


void SC_output(){

  printf("# elapsed time= %f \n",omp_get_wtime()-TIME0);
  
  int i,j;
  for (i=0; i < N_GROUP; i++){
    for (j=i+1;  j < N_GROUP; j++){
      printf("%d ",C[i][j]);};
  };
  printf("\n");
  fflush(NULL);
  
};

