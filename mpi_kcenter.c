#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <mpi.h>

#define MAX_POINTS 200

// read the file of points
int read_file (double* points, char* filename) {

    int n = 0;
    double next[2];
    FILE* file_ptr;
    file_ptr = fopen(filename,"r");
    if (file_ptr == NULL) {
        printf ("error : could not open file %s for reading\n",filename);
        exit(1);
    }
    while (fscanf (file_ptr,"%lf %lf",next,next+1) == 2) {
        if (n < MAX_POINTS) {
            points[2*n] = next[0];
            points[2*n+1] = next[1];
            n += 1;
        } else {
            printf ("Too many points in file %s\n",filename);
            fclose (file_ptr);
            exit(1);
        }
    }
    fclose (file_ptr);
    return n;
}

// calculate the distance squared between two points
double calc_dist_sq (double* u, double* v) {
    double diff_x = u[0] - v[0];
    double diff_y = u[1] - v[1];
    return (diff_x*diff_x + diff_y*diff_y);
}

// calculate the cost squared of a given set of center locations
double center_cost_sq (double* dist_sqs, int n, int* centers, int k, double min_cost_sq) {
    double cost_sq = 0;
    for (int i=0;i<n;i++) {
        double min_dist_sq = DBL_MAX;
        for (int j=0;j<k;j++) {
            double dist_sq = dist_sqs[i*n+centers[j]];
            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
            }
        }
        // check to see if we can abort early
        if (min_dist_sq > min_cost_sq) {
            return min_dist_sq;
        }
        if (min_dist_sq > cost_sq) {
            cost_sq = min_dist_sq;
        }
    }
    return cost_sq;
}

int main (int argc, char** argv) {

    MPI_Init (&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    // check command line arguments : filename, k, m, and s
    if (argc < 5) {
        printf ("Command usage : %s %s %s %s %s\n",argv[0],"filename","k","m","s");
        return 1;
    }

    // the points array
    double points[2*MAX_POINTS];

    // read dataset
    int n = read_file (points,argv[1]);

    // calculate the distances squared table
    double dist_sqs[n*n];
    for (int i= 0; i<n; i++) {
        for (int j=0;j<n;j++) {
            dist_sqs[i*n+j] = calc_dist_sq (points+2*i,points+2*j);
        }
    }

    // get k and m from the command line
    int k = atoi(argv[2]);
    int m = atoi(argv[3]);

    // seed the random number generator using command line seed
    srandom(atoi(argv[4]) + rank);

    // start the timer
    double start_time;
    start_time = MPI_Wtime();

    // check the cost_sq of m random sets of k centers
    int centers[k];
    int optimal_centers[k];
    double min_cost_sq = DBL_MAX;
    double cost_sq;
    for (int i = 0; i < m; i ++) {
        for (int j=0;j<k;j++) {
            centers[j] = random() % n;
        }
        cost_sq = center_cost_sq(dist_sqs,n,centers,k,min_cost_sq);
        if (cost_sq < min_cost_sq) {
            min_cost_sq = cost_sq;
            for (int j=0;j<k;j++) {
                optimal_centers[j] = centers[j];
            }
        }
    }

    if (rank == 0) {
        int rank_optimal[k];
        MPI_Status status;
        for (int source = 1; source < size; source++) {
            MPI_Recv(rank_optimal, k, MPI_INT, source, 0, MPI_COMM_WORLD, &status);
            double cost_sq = center_cost_sq(dist_sqs, n, rank_optimal, k, min_cost_sq);
            if (cost_sq < min_cost_sq) {
                min_cost_sq = cost_sq;
                for (int i = 0; i < k; i++) {
                    optimal_centers[i] = rank_optimal[i];
                }
            }
        }
    } else {
        MPI_Send(optimal_centers, k, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    // stop the timer
    double end_time;
    end_time = MPI_Wtime();

    if (rank == 0) {
        // print wall time
        printf ("wall time used = %.4f sec\n",(end_time-start_time));

        // print out the total number of k-tuples checked
        printf("total number of %d-tuples checked = %lld\n", k, (long long)(m * size));


        // print the approximate minimal cost for the k-center problem
        printf ("approximate minimal cost = %.4lf\n",sqrt(min_cost_sq));

        // print an approx optimal solution to the k-center problem
        printf ("approximate optimal centers : ");
        for (int j=0;j<k;j++) {
            printf ("%d ",optimal_centers[j]);
        }
        printf ("\n");
    }

    MPI_Finalize();
    return 0;
}
