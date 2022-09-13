#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

int main() {

	double tau0, tau1;
	tau0 = time(NULL);
	double m = 1, k = 1, x = 0, v = 1, t_max = 100, dt = 0.1;

	vector<double> x_list, v_list, t_list;

	for (int t = 0; t <= t_max; t = t + dt) {
		double a;
		t_list.push_back(t);
		x_list.push_back(x);
		v_list.push_back(v);
		a = -k * x / m;
		x = x + v * dt;
		v = v + a * dt;
	}
	tau1 = time(NULL);
	cout << endl << "The time of the operation is:    " << tau1 - tau0;
	return 0;
}