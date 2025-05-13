#include <chrono> 
#include <iostream> 

// Note: include -fopenmp tag while compiling

int sum_serial(int n) 
{ 
	int sum = 0; 
	for (int i = 0; i <= n; ++i) { 
		sum += i; 
	} 
	return sum; 
} 

// Parallel programming function 
int sum_parallel(int n) 
{ 
	int sum = 0; 
    #pragma omp parallel for reduction(+ : sum) 
	for (int i = 0; i <= n; ++i) { 
		sum += i; 
	} 
	return sum; 
} 

// Driver Function 
int main() 
{ 
	const int n = 100000000; 

	auto start_time 
		= std::chrono::high_resolution_clock::now(); 

	int result_serial = sum_serial(n); 

	auto end_time 
		= std::chrono::high_resolution_clock::now(); 

	std::chrono::duration<double> serial_duration 
		= end_time - start_time; 

	start_time = std::chrono::high_resolution_clock::now(); 

	int result_parallel = sum_parallel(n); 
	end_time = std::chrono::high_resolution_clock::now(); 
	std::chrono::duration<double> parallel_duration 
		= end_time - start_time; 

	std::cout << "Serial result: " << result_serial 
			<< std::endl; 
	std::cout << "Parallel result: " << result_parallel 
			<< std::endl; 
	std::cout << "Serial duration: "
			<< serial_duration.count() << " seconds"
			<< std::endl; 
	std::cout << "Parallel duration: "
			<< parallel_duration.count() << " seconds"
			<< std::endl; 
	std::cout << "Speedup: "
			<< serial_duration.count() 
					/ parallel_duration.count() 
			<< std::endl; 

	return 0; 
}
