// Coding by Oscar
//
//
//
#include "mpi.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>

// Global variables
long a = 1664525;
long m = 4294967296;
long c = 1013904223;

/*
Method: Modulo exponentiation
Formula : x**n mod p
*/
long calExpo(long x, long n, long p)
{
    long y = 1;
    for (int i = 0; i < n; ++i)
    {
        y = ( y * x) % p;
    }
    return y;
}

/*
Method: Calculate the constant A used for frog leaping
Formula: A = a**k mod m
*/
long calA(long a, int k, long m)
{
    return calExpo(a, k, m);
}

/*
Method: Calculate the constant C used for frog leaping
Formula: C = c*(a**(k-1) + a**(k-2) +.....+ a**1 + a**0) mod m
*/
long calC(long a, int k, long m)
{
    long sumf = 0;
    for (int i = k; i >= 0; --i)
    {
        sumf += calExpo(a, i, m);
    }
    long res = (c * (sumf % m)) % m;
    return res;
}

/*
Method: Generate the random number in parallel evironment
Formula: x(i+k) = (A * xi + C) mod m
Variable Desciption:
world_size: the numbers of process
length: the total length of random number calculated
processorNo: the index of processor
*/

double generateAndCRandomNumInParallel(int seed, int world_size, int length, int processorNo)
{
    double sumf = 0;
    long i_prev = seed; // assign a seed
    long A = calA(a, world_size, m);
    long C = calC(a, world_size, m);
    for (int i = processorNo; i < length - world_size; i = i + world_size)
    {
        long i_next = (A * i_prev + C) % m;
        double x_rand = float(i_next) / (m - 1); // Scale to a float in the range 0-1
        i_prev = i_next;
        double fx = sqrt(1 - x_rand * x_rand);
        sumf += fx;
    }
    double area = sumf / ((length - world_size)/ world_size);
    double pi_calc = 4 * area;
    return pi_calc;
}

/*
Method: Generate the random number sequencially
Formula: x(i+1) = (a * xi + c) mod m
Variable Desciption:
world_size: the numbers of process
length: the total length of random number calculated
array[]: output parameters which is used for saving the random numbers generated lately
*/

long* generateRandomNumInSequence(const long seed, const long length)
{
    long *sequenceSeed = new long[length];
    long i_prev = seed;
    long i_next = 0;
    for (int i = 0; i < length;++i)
    {
        i_next = (a * i_prev + c) % m; // Next Random integer in the sequence
        sequenceSeed[i] = i_next;
        i_prev = i_next;
    }
    return sequenceSeed;
}




int main(int argc, char *argv[])
{
    double t_start, t_start_p1, t_start_p2, t_end;
    double t_part1, t_part2;
    double result = 0, ret0 = 0, ret1 = 0;
    const long length = atol(argv[1]);
    const long seed = atol(argv[2]);
    int world_size;
    int rank;

    MPI_Status Stat;// status variable, so operations can be checked
    MPI_Init(&argc, &argv);
    t_start = MPI_Wtime(); // set the start time
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        std::cout<<"------------------------Parallel Version------------------------"<<std::endl;
        std::cout<<"================================================================"<<std::endl;
        // 1.Generate seeds of corresponding numbers with the numbers of processors
        long *arraySeed = new long[world_size];
        arraySeed = generateRandomNumInSequence(seed, world_size);
        //fprintf(stdout, "processor[%d] intital seed: %ld\n", rank, arraySeed[rank]);

        // 2. Master sends seed to all the slave processes
        for (int i = 1; i < world_size; i++)
        {
            fprintf(stdout, "processor[%d] intital seed: %ld\n", i, arraySeed[i]);
            MPI_Send(&arraySeed[i], 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
        }

        // 3. Master calculates its own partial pi
        ret0 = generateAndCRandomNumInParallel(arraySeed[0], world_size, length, rank );
        //fprintf(stdout, "processor[%d] return value: %.10f\n", rank, ret0);
        result += ret0;
        t_start_p2 = MPI_Wtime();

        // 7. Receive nodes from all nodes
        for (int i = 1; i < world_size; i++)
        {
            MPI_Recv(&ret1, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Stat);
            result += ret1;
        }

        result = result / world_size;
        t_end = MPI_Wtime();    //set the end time
        //t_part2 = t_end - t_start_p2;
        //fprintf(stdout, "----------Time before master's receiving return value: %.10f----------------------------------\n", t_start_p2);
        //fprintf(stdout, "----------End time: %.10f----------------------------------\n", t_end);
        //fprintf(stdout, "----------Time of sequential part 2: %.10f----------------------------------\n", t_part2);


        fprintf(stdout, "----------Time elapsed: %.10f----------------------------------\n", t_end - t_start);
        fprintf(stdout, "----------PI value calculated by %ld random numbers: %.10f----------------------------------\n", length, result);
        std::cout <<"================================================================"<<std::endl;

    }
    else // this is not the master
    {
        long slaveseed;
        // 4. slave processor wait for receiving the original seed
        MPI_Recv(&slaveseed,1 , MPI_LONG, 0, 0, MPI_COMM_WORLD, &Stat);
        //t_start_p1 = MPI_Wtime();
        //fprintf(stdout, "**********slave time after receving %.10f \n", t_start_p1);
        //t_part1 = t_start_p1 - t_start;
        //fprintf(stdout, "----------Start time:  %.10f----------------------------------\n", t_start);
        //fprintf(stdout, "----------Time after all slaves' receiving: %.10f----------------------------------\n", t_start_p1);
        //fprintf(stdout, "----------Time of sequential part 1: %.10f----------------------------------\n", t_part1);
        // 5. slave processor calculate their own pi
        ret1 = generateAndCRandomNumInParallel(slaveseed, world_size, length,  rank);
        //fprintf(stdout, "processor[%d] return value: %.10f\n", rank, ret1);

        // 6. slave processor send the pi value back to master processor
        MPI_Send(&ret1, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    //fprintf(stdout, "---------- %.10f %.10f %.10f %.10f----------------------------------\n", t_part1,t_part2 ,t_end, t_start );
    //t_end = MPI_Wtime();
    //fprintf(stdout, "----------Time percetage of sequential part of the whole program: %.10f----------------------------------\n", (t_part1+t_part2)/(t_end - t_start) );
    MPI_Finalize();
    fflush(stdout);
    return 0;
}


