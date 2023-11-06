#include <iostream>
#include <vector>

#include "matrix_class.h"
#include "polynom_class.h"
#include "integration_methods.h"


using namespace std;
using namespace biv;

int main() {
	vector<double> nodes;
	double ans = computeIntegral<double>(makeSFGauss<double>(3, nodes), nodes);
	cout << ans << '\n';
	cout << "=================\n";

	double ans2 = computeIntegral<double>(makeISF<double>(3, nodes), nodes);
	cout << ans2 << "\n";
	
}