#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <sys/time.h>

double *allocMatrix( size_t n)
{
   double *m;
   m = (double *)malloc( n*n*sizeof(double));
   if( m==NULL) {
      perror( "failed to allocate matrix; ");
   }
   return m;
}

void print( double *out, size_t n, int rows)
{
   size_t i,j,maxn, maxl;
   maxn = (n < 20 ? n : 20);
   maxl = (rows < 20 ? rows : 20);
   for( i=0; i<maxl; i++) {
      printf( "|");
      for( j=0; j<maxn; j++) {
         printf( " %7.2f", out[i*n+j]);
      }
      if( maxn < n) {
         printf( "...|\n");
      } else {
         printf( "|\n");
      }
   }
   if( maxl < rows) {
         printf( "...\n");
      }
}
void init( double *out, size_t n)
{
   size_t i,j;

   for( i=0; i<n; i++) {
      for( j=0; j<n; j++) {
         out[i*n+j] = 0;
      }
   }

}

//Cuts up a given matrix given the start and end values. 
double* get_part(double *out, size_t n, int start, int end)
{
    int len = end - start;
    int i = 0;
    double *m;
    
    m = (double *)malloc(len*sizeof(double));
    
    
    while (i < len){
        
        m[i] = out[start+i];
        i++;
        
    }

    return m;
}

//Returns the amount of rows that a given rank has to calculate, the final rank gets the remainder
int get_row_number(size_t n, int total_nodes, int node_number){
    int rows_per_node = floor((n-2)/total_nodes);
    int remainder = (n-2) % total_nodes;
    if (node_number == total_nodes-1){
        return rows_per_node+remainder+2;
    }
    return rows_per_node + 2;
}


//Splits the full matrix into the different parts that go to the different ranks, the final rank gets the remainder
double* split(double *out, size_t n, int total_nodes, int node_number)
{
    
    int rows_per_node = floor((n-2)/total_nodes);
    int remainder = (n-2) % total_nodes;
    
    if (node_number==0){
        double* m = get_part(out, n, 0, (rows_per_node+2)*n);
        return m;
        
    }else if (node_number==total_nodes-1){
        double* m = get_part(out, n, n*n-(rows_per_node+remainder+2)*n, n*n);
        return m;
    }
    
    double* m = get_part(out, n, node_number*rows_per_node*n, ((node_number+1)*rows_per_node+1)*n+n);
    return m;
}

//Updates the array slice with the values obtained from the other ranks
double* update(double *out, size_t n, double *update_start, double *update_end,int total_nodes, int node_number) {
    
    int length = get_row_number(n, total_nodes, node_number);
    
    if (node_number != 0) {
        for (int i=0; i < n; i++) {
            out[i] = update_start[i];
            
        }
    }
    if (node_number != total_nodes-1) {
        for (int j=0; j<n; j++){
            out[(length*n-n)+j] = update_end[j];
            
        }
        
    }
    return out;
}

void relax( double *in, double *out, size_t n, int total_nodes, int node_number)
{
   size_t i,j;
   int rows = get_row_number(n, total_nodes, node_number);
    
   for( i=1; i<rows-1; i++) {
      for( j=1; j<n-1; j++) {
         out[i*n+j] = 0.25*in[(i-1)*n+j]      // upper neighbour
                      + 0.25*in[i*n+j]        // center
                      + 0.125*in[(i+1)*n+j]   // lower neighbour
                      + 0.175*in[i*n+(j-1)]   // left neighbour
                      + 0.2*in[i*n+(j+1)];    // right neighbour
      }
      
   }
  
}

//Copies the values of one array into another array
void copy(double *in, double *out, size_t n, int total_nodes, int node_number){
    int rows = get_row_number(n, total_nodes, node_number);
    for (int i = 0; i< rows*n; i++) {
        out[i] = in[i];
        
    }
    
    
}

void get_n_elements(double *in, double *out, size_t n, int begin, int end){
    
    for (int i = 0; i < (end-begin); i++){
        out[i] = in[i+begin];
        
    }
}

void update_slice(double *in, double *out, size_t n, int begin, int end) {
    
    for (int i = 0; i < (end-begin); i++){
        out[i+begin] = in[i];
        
    }
    
}


int main (int argc, char *argv[])
{
   struct timespec t1, t2, t3, t4;
   double *a,*b, *c, *tmp;
   size_t n=0;
   int i;
   int max_iter;
   double* m;
   
   
    
   if( argc < 3) {
      printf("call should have two arguments \"%s <n> <iter>\"\n", argv[0]);
      exit(1);
   }
   if( sscanf( argv[1], "%zu", &n) != 1) {
      printf("non size_t value for matrix size\n");
      exit(1);
   }
  
   if( sscanf( argv[2], "%d", &max_iter) != 1) {
      printf("non int value for # iterations\n");
      exit(1);
   }
    
   //Sets up the MPI structure and gets the total ranks and the rank of the specific core. 
   int rank = 0; 
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   int total_ranks = 1; 
   
   MPI_Comm_size(MPI_COMM_WORLD, &total_ranks);
   if (rank == 0 && clock_gettime(CLOCK_MONOTONIC, &t1) != 0) {
      return -1;
   }
    
   printf("Running on %i Cores, reporting from node %i \n", total_ranks, rank);
    
   //Allocates memory for the different matrices used. 
   b = (double *)malloc(get_row_number(n, total_ranks, rank)*n*sizeof(double));
   a = allocMatrix( n);
   c = allocMatrix(n);
   if (rank ==0){
       
       init( a, n);
       init( c, n);

       a[n/4] = 100.0;
       c[n/4] = 100.0;

       a[(n*3)/4] = 1000.0;
       c[(n*3)/4] = 1000.0;

       printf( "size   : n = %zu => %d M elements (%d MB)\n",
               n, (int)(n*n/1000000), (int)(n*n*sizeof(double) / (1024*1024)));
       printf( "iter   : %d\n ", max_iter);
       print(a, n, n);
   }
     

   if (rank == 0) {
       
       for (i = 1; i< total_ranks; i++) {
           //Send the slices to the right ranks
           m = split(a, n, total_ranks, i);
           MPI_Send(m, get_row_number(n, total_ranks, i)*n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
           
       }
       
       m = split(a, n, total_ranks, 0);
       
   }else {
       //Receive the matrix slices from node 0
       m = (double *)malloc(get_row_number(n, total_ranks, rank)*n*sizeof(double));
       MPI_Recv(m, get_row_number(n, total_ranks, rank)*n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
       
     
   }
   
   //Copy the original slice so that the input and output of the relax is correct. 
   copy(m,b, n, total_ranks, rank);
   if (rank == 0 && clock_gettime(CLOCK_MONOTONIC, &t2) != 0) {
      return -1;
   }
   for(i=0; i<max_iter; i++) {
      tmp = m;
      m = b;
      b = tmp;
      relax(m, b, n, total_ranks, rank);
      
      if (rank != 0) {
          
          //Create the two slices that will be sent to the other ranks
          double* begin = (double *)malloc(n*sizeof(double));
          get_n_elements(b, begin, n, n, 2*n);
                  
          //Send the two slices that changed to the ranks that care about them
          //printf("%i sending to %i \n", rank, rank-1);
          MPI_Send(begin, n, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD); 
                   
          //And receive the updates slices from the other ranks
          //printf("%i receiving from %i \n", rank, rank-1);
          MPI_Recv(begin, n, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
         
          //Update the slices
          update_slice(begin, b, n, 0, n);
         
      }
      if (rank != total_ranks-1){
         
          //Create the two slices that will be sent to the other ranks         
          double* end = (double *)malloc(n*sizeof(double));
          double* end2 = (double *)malloc(n*sizeof(double));
          get_n_elements(b, end, n, (get_row_number(n, total_ranks, rank)-2)*n, (get_row_number(n, total_ranks, rank)-1)*n);
          
          //And receive the updates slices from the other ranks 
          //printf("%i receiving from %i \n", rank, rank+1);
          MPI_Recv(end2, n, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          
          //Send the two slices that changed to the ranks that care about them       
          //printf("%i sending to %i \n", rank, rank+1);
          MPI_Send(end, n, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
          
         
          
          //Update the slices
          update_slice(end2, b,n,  (get_row_number(n, total_ranks, rank)-1)*n, (get_row_number(n, total_ranks, rank)*n));
         
          
      }
   }
    if (rank == 0 && clock_gettime(CLOCK_MONOTONIC, &t3) != 0) {
      return -1;
   }
   if (rank != 0) {
       //Send the final matrix slices back to node 0
       MPI_Send(b, get_row_number(n, total_ranks, rank)*n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
       free(b);
       MPI_Finalize();
       
   }else {
        update_slice(b, c, n, 0, get_row_number(n, total_ranks, 0)*n);
      
        for (i = 1; i< total_ranks; i++) {
           
           //Receive the matrix slices and put the final matrix back together
           if (i == total_ranks-1) {
                b = (double *)malloc(get_row_number(n, total_ranks, i)*n*sizeof(double));
           }
           MPI_Recv(b, get_row_number(n, total_ranks, i)*n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         
           if (i == total_ranks-1){
               update_slice(b, c, n, (n*n - get_row_number(n, total_ranks, i)*n) , n*n);
           }else {
               update_slice(b, c, n, (get_row_number(n, total_ranks, i)-2)*n*i , ((get_row_number(n, total_ranks, i)-2)*(i+1)+2)*n);
           }
       }
       
       
       //Printing the finalized matrix
       printf("\n");
       print(c, n, n);
       printf("\n");
       MPI_Finalize();
       
       
       double timeSpendRelax = ((double)t3.tv_sec - (double)t2.tv_sec)
                                + ((double)t3.tv_nsec - (double)t2.tv_nsec) / 1000000000.0;

       double gflop = (5.0 * (double)n * (double)n * (double)max_iter) / 1000000000.0;
       if (clock_gettime(CLOCK_MONOTONIC, &t4) != 0) {
          return -1;
       }
       double timeSpendTotal = ((double)t4.tv_sec - (double)t1.tv_sec)
      + ((double)t4.tv_nsec - (double)t1.tv_nsec) / 1000000000.0;
       printf("Core_count: %i\n", total_ranks);
       printf("Matrix_size: %zu\n", n);
       printf("Iter   : %d\n", max_iter);
       printf("Size: %.6f\n", gflop / timeSpendRelax);
       printf("Time spend relax: %.6f seconds\n", timeSpendRelax);
       printf("Time spend total: %.6f seconds\n", timeSpendTotal);
       printf("Computation: %.6f GFLOP/s\n", gflop / timeSpendRelax);
   }
   
   
}