#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))

int main (int argc, char *argv[])
{
   unsigned long int    count;        /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   unsigned long int    first;        /* Index of first multiple */
   int   local_first;
   unsigned long int    global_count = 0; /* Global prime count */
   unsigned long long int    high_value;   /* Highest value on this proc */
   unsigned long int    i;
   int    id;           /* Process ID number */
   unsigned long int    index;        /* Index of current prime */
   unsigned long long int    low_value;    /* Lowest value on this proc */
   char  *marked;       /* Portion of 2,...,'n' */
   char  *local_prime_marked;
   unsigned long long int    n;            /* Sieving from 2, ..., 'n' */
   int    p;            /* Number of processes */
   unsigned long int    proc0_size;   /* Size of proc 0's subarray */
   unsigned long int    prime;
   unsigned long int  local_prime;        /* Current prime */
   unsigned long int    size;         /* Elements in 'marked' */
   unsigned long int  local_prime_size;


   MPI_Init (&argc, &argv);

   /* Start the timer */

   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();

   if (argc != 2) {
      if (!id) printf ("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit (1);
   }

   n = atoll(argv[1]);

   /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */

   /* Add you code here  */

   low_value = 2 + id * (n - 1) / p;
   high_value = 1 + (id + 1) * (n - 1) / p;

   // make sure low_value is odd
   low_value = (low_value % 2) == 0 ? low_value : low_value + 1;
   // make sure high_value is odd
   high_value = (high_value % 2) == 0 ? high_value : high_value - 1;
   
   // total number except even
   size = (high_value - low_value) / 2 + 1;
   
   local_prime_size = (int)sqrt((double)(n)) - 1;

   /* Bail out if all the primes used for sieving are
      not all held by process 0 */

   proc0_size = (n/2 - 1) / p;

   if ((2 + proc0_size) < (int) sqrt((double) n/2)) 
   {
       if (!id) printf("Too many processes\n");
       MPI_Finalize();
       exit(1);
   }

   /* Allocate this process's share of the array. */

   marked = (char *) malloc(size);
   local_prime_marked = (char *) malloc(local_prime_size);

   if (marked == NULL) 
   {
       printf("Cannot allocate enough memory\n");
       MPI_Finalize();
       exit(1);
   }

   // find all prime numbers within the sqrt of n   
   for (i = 0; i < local_prime_size; i++) local_prime_marked[i] = 0;
   index = 0;
   local_prime = 2;
   do
   {
       local_first = local_prime * local_prime - 2;
       for (i = local_first; i < local_prime_size; i += local_prime) local_prime_marked[i] = 1;
       while (local_prime_marked[++index] == 1);
       local_prime = 2 + index;
   } while (local_prime * local_prime <= n);
   
   // set the size of block
   unsigned long int size_block = 500000;
   // calculate the number of block
   unsigned long int num_block = size / size_block;
   num_block = ((size % size_block) == 0) ? num_block : num_block + 1;
   
   unsigned long int orignal_low_value = low_value;
   unsigned long int orignal_high_value = high_value;
   // let low_value be the low value of block
   low_value = low_value;
   // let high_value be the high value of block
   high_value = low_value + 2 * (size_block - 1);
   
   // find all prime numbers within n
   for (i = 0; i < size; i++) marked[i] = 0;
   while (num_block--)
   {   
       index = 0;
       prime = 3;  
       do {
           if (prime * prime > low_value)
           {   
               first = (prime * prime - low_value) / 2;
           }   
           else 
           {    
               if (!(low_value % prime)) 
               {
                  first = 0;
               }
               else
               {   
                  first = (low_value / prime + 1) * prime;
                  // make sure first is odd
                  first = ((first - low_value) % 2) == 0 ? first : first + prime;
                  // get the index of array
                  first = (first - low_value) / 2;
               }   
           }
           for (i = first + (low_value - orignal_low_value) / 2; i <= (high_value - orignal_low_value) / 2; i += prime) 
              marked[i] = 1;
           //if (!id) {
           while (local_prime_marked[++index]);
           // get next prime
           prime = index + 2;
           //}
           //if (p > 1) MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
       } while (prime * prime <= high_value); 
       
       // calculate the low_value and high_value for next block
       orignal_low_value = low_value;
       low_value = high_value + 2;
       high_value = low_value + 2 * (size_block - 1);
       if (orignal_high_value < high_value)
       {   
          high_value = orignal_high_value;
       } 
   }
   
   count = 0;
   for (i = 0; i < size; i++)
       // calculate the number of prime
       if (!marked[i]) count++;
   if (p > 1)
       MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   
   global_count++;
   
   /* Stop the timer */

   elapsed_time += MPI_Wtime();
   
   /* Print the results */

   if (!id) {
      printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count, elapsed_time, p);

   }
   MPI_Finalize ();
   return 0;
}

