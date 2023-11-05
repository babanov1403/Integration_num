#include <vector>
#include "matrix_class.h"
#include <iostream>
#include <cmath>
//#include "burger.h"
namespace var_b {
	static const double alpha = 1.0 / 3;
	static const double a = 1.5;
	static const double b = 3.3;
}
/*
static const double alpha = 1;
static const double a = 0;
static const double b = 3;*/

double func(double x) {
	return (2.0 * cos(2.5 * x) * exp(x / 3.0) + 4 * sin(3.5 * x) * exp(-3 * x) + x);
}
double func_test(double x) {
	return x * x;
}

namespace biv {
	template <typename Number>
	Matrix<Number> makeISF(const vector<double>& nodes);
	template <typename Number>
	Matrix<Number> makeSFGauss(const vector<double>& nodes);
	double computeMoment(double k);
	template <typename Number>
	double computeIntegral(Matrix<Number> Mu, const vector<double>& nodes);
	template <typename Number>
	Matrix<Number> RhimannIntegrale(double, double, int, Matrix<Number>(*MatrixFunc)(double ksi));

	vector<double> makeNodes(int n, double a = var_b::a, double b = var_b::b) {
		vector<double> nodes(n);
		double h = (b - a) / n;
		nodes[0] = a;
		for (int i = 1; i < n; i++)
			nodes[i] = nodes[i - 1] + h;
		return nodes;
	}

	template <typename Number>
	Matrix<Number> makeISF(const vector<Number>& nodes) {
		Matrix<Number> Mu(nodes.size(), 1);
		for (int i = 0; i < nodes.size(); i++)
			Mu(i, 0) = computeMoment(static_cast<double>(i));
		Matrix<Number> A(nodes.size(), nodes.size());
		for (int i = 0; i < nodes.size(); i++)
			for (int j = 0; j < nodes.size(); j++)
				A(i, j) = pow(nodes[j] - var_b::a, i); //t = x-a   f(xj - a)
		Matrix<Number> Outp = GaussSlau(A, Mu);
		return Outp;
	}
	double computeMoment(double k) {
		k -= var_b::alpha;
		double ans;
		if (k == -1) ans = log(var_b::b - var_b::a);
		else {
			ans = pow((var_b::b - var_b::a), k + 1);
			ans /= k + 1;
		}
		//cout << ans << " " << k << '\n';
		return ans;
	}

	template <typename Number>
	double computeIntegral(Matrix<Number> Aj, const vector<double>& nodes) {
		int n = Aj.n;
		double integrale;
		vector<double> sl(n);
		for (int i = 0; i < n; i++) {
			sl[i] = Aj(i, 0) * func(nodes[i]);// xj = tj + a(a a + b / 2 b)
		}
		sort(sl.begin(), sl.end());
		for (auto i : sl) integrale += i;
		return integrale;
	}

	template <typename Number>
	Matrix<Number> makeSFGauss(int n) {
		Matrix<Number> Mu(2 * n, 1);
		for (int i = 0; i < 2 * n; i++) {
			Mu(i, 0) = computeMoment(i);
			cout << i << '\n';
		}
		cout << Mu;
		Matrix<Number> A(n, n);
		Matrix<Number> B(n, 1);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				A(i, j) = Mu(i + j, 0);
				B(i, 0) = -Mu(n + i, 0);
			}
		}
		cout << A;
		cout << B;
		Matrix<Number> A_ = GaussSlau(A, B);//a[i] of our poly		
	}
	/*template <typename Number>
	Matrix<Number> RhimannIntegrale(double x0, double x1, int N = 1e5, Matrix<Number> (*MatrixFunc)(double ksi)) {
		double h = (x1 - x0) / N;
		Matrix<Number> res = MatrixFunc(x0);
		Matrix<Number> curr = res;
		for (double ksi = x0+h; ksi <= x1; ksi += h) {
			curr = MatrixFunc(ksi);
			curr *= h;
			res += h;
		}
		return res;
	}*/
}

