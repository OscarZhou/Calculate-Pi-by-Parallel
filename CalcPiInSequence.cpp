// Coding by Oscar
//
//
//
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>

// Global variables
long a = 1664525;
long m = 4294967296;
long c = 1013904223;

/*
 * Method: Generate the random number sequencially
 * formula: x(i+1) = (a * xi + c) mod m
 * Variable DEscription:
 */
double generateAndCalculatePiInSequence(const long seed, const long length)
{

	double sumf = 0;
	long i_prev = seed;
	long i_next = 0;
    for (long i = 0; i < length;++i)
	{
		i_next = (a * i_prev + c) % m;
        double x_rand = float(i_next) / (m - 1);
		i_prev = i_next;
		double fx = sqrt(1 - x_rand * x_rand);
		sumf += fx;
	}
	double area = sumf / length;
	double pi_calc = 4 * area;
	return pi_calc;
}


int main(int argc, char* argv[])
{
    double result = 0;
    const long length = atol(argv[1]);
    const long initial_seed = atol(argv[2]);

    long t_start = clock();
	result = generateAndCalculatePiInSequence(initial_seed, length);
    long t_end = clock();

    std::cout<<"------------------------Sequential Version------------------------"<<std::endl;
    std::cout<<"================================================================"<<std::endl;
    std::cout << "----------Start time: "<<t_start <<"----------------------------------"<<std::endl;
    std::cout << "----------End time: "<<t_end <<"----------------------------------"<<std::endl;
    std::cout << "----------Time elapsed: "<< (t_end - t_start)/double(CLOCKS_PER_SEC)<<" s------------------------"<<std::endl;
    std::cout << "----------PI value calculated by "<< length<<" random numbers: "<< result<<"--------------------"<<std::endl;
    std::cout <<"================================================================"<<std::endl;
	return 0;
}
