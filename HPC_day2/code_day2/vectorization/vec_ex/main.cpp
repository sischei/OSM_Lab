#include <time.h>
#include <iostream>

int main( int argc, const char* argv[] ) {

	int n = 2048;

	double a[n];
	double b[n];
	float a1[n];
	float b1[n];

	for (int i = 0; i < n; i++) {
		a[i] = 3.0;
		b[i] = 7.0;
		a1[i] = 3.0;
		b1[i] = 7.0;
	}

	clock_t t = clock();
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			a[i] += b[i];
		}
	}
	t = (clock() - t);

	clock_t t1 = clock();
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			a1[i] += b1[i];
		}
	}
	t1 = (clock() - t1);

	std::cout <<  "CPU Clock [ double , float , double/float ] : [" << t << "," << t1 << "," << (double(t) / double(t1)) << "]" << std::endl;
	return 1;
}