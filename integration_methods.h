#include <vector>
#include <iostream>
#include <cmath>
#include <numeric>
#include <fstream>
#include <functional>
//#include "burger.h"
#include "matrix_class.h"

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
	double computeMoment(double k, double a, double b);
	template <typename Number>
	double computeIntegral(Matrix<Number> Mu, const vector<double>& nodes);
	//template <typename Number>
	//Matrix<Number> RhimannIntegrale(double, double, int, Matrix<Number>(*MatrixFunc)(double ksi));

	vector<double> makeNodes(int n, double a = var_b::a, double b = var_b::b) {
		vector<double> nodes(n);
		double h = (b - a) / (n-1);
		nodes[0] = a;
		for (int i = 1; i < n; i++)
			nodes[i] = nodes[i - 1] + h;
		return nodes;
	}
	
	template <typename Number>
	Matrix<Number> makeISF(int n, vector<double>& nodes, double a, double b) {
		nodes = makeNodes(n, a, b);
		//for (auto i : nodes) cout << i << ' ';
		//cout << '\n';
		Matrix<Number> Mu(nodes.size(), 1);
		for (int i = 0; i < nodes.size(); i++)
			Mu(i, 0) = computeMoment(i, a, b);
		Matrix<Number> A(nodes.size(), nodes.size());
		for (int i = 0; i < nodes.size(); i++)
			for (int j = 0; j < nodes.size(); j++)
				A(i, j) = pow(nodes[j] - var_b::a, i); //t = x-a   f(xj - a)
		Matrix<Number> Ashki = GaussSlau(A, Mu);
		bool flag = false;
		double Ashki_sum = 0, XD_A_sum = 0;
		for (int i = 0; i < n; i++) {
      		Ashki_sum += abs(Ashki(i, 0));
			XD_A_sum += Ashki(i, 0);
			if (Ashki(i, 0) < 0) flag = true;
		}
		cout << Ashki_sum << '\n';
		//fout << Ashki_sum << ',' << n << '|';
		//if (flag) cout << "\nthere is some negative A[i]!!!\n";
		return Ashki;
	}
	double computeMoment(double k, double a, double b) {
		k -= var_b::alpha;
		double ans;
		if (k == -1) {
			cout << "\nSHIT\n";
			ans = 1e9;
		}else {
			ans = pow((b - var_b::a), k + 1);
			ans /= k + 1;
			ans -= pow((a - var_b::a), k + 1)/(k+1);
		}
		return ans;
	}
	//important to note rhat nodes must lay between a...b
	template <typename Number>
	double computeIntegral(Matrix<Number> Aj, const vector<double>& nodes) {
		int n = Aj.n;
		double integrale = 0;
		vector<double> sl(n);
		for (int i = 0; i < n; i++) {
			sl[i] = Aj(i, 0) * func(nodes[i]);// xj = tj + a(a a + b / 2 b)
		}
		sort(sl.begin(), sl.end());
		for (auto i : sl) {
			integrale += i;
		}
		return integrale;
	}
	template <typename Number>
	Matrix<Number> makeSFGauss(int n, vector<double>& nodes, double a, double b) {
		Matrix<Number> Mu(2 * n, 1);
		for (int i = 0; i < 2 * n; i++) {
			Mu(i, 0) = computeMoment(i, a, b);
		}
		Matrix<Number> A(n, n);
		Matrix<Number> B(n, 1);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				A(i, j) = Mu(i + j, 0);
				B(i, 0) = -Mu(n + i, 0);
			}
		}
		Matrix<Number> A_ = GaussSlau(A, B).transpose();//a[i] of our poly	
		vector<Number> vec = A_[0];
		vec.push_back(1);
		Polynom<Number> p(vec);
		double tmp;
		nodes = p.findRoots(tmp);
		n = nodes.size();
		Matrix<Number> A_1(n, n);
		Matrix<Number> B_1(n, 1);
		for (int i = 0; i < n; i++) {
			B_1(i, 0) = Mu(i, 0);
			for (int j = 0; j < n; j++)
				A_1(i, j) = my_pow(nodes[j], i);
		}
		for (auto& q : nodes) q += var_b::a;
		Matrix<Number> Ashki = GaussSlau(A_1, B_1);
		bool flag = false;
		double Ashki_sum = 0;
		for (int i = 0; i < n; i++) {
			Ashki_sum += abs(Ashki(i, 0));
			if (Ashki(i, 0) < 0) flag = true;
		}
		//if (flag) cout << "\nthere is some negative A[i]!!!\n";
		return Ashki;
	}

	using vpairdd = vector<pair<double, double>>;
	// n - суммарное кло-во узлов, p - кол-во интервалов, rule - просто рофло-функция
	double makeCompoundSF(int n, int p, const function<vpairdd(int)>& rule) {
		vpairdd interv = rule(p);
		//только ньютон-котс
		vector<vector<double>> nodes(p);
		vector<Matrix<double>> Ashki(p);
		vector<int> node_counter(p, 0);
		if (interv.size() != p) {
			throw runtime_error("");
		}
		//медленный алгос, просто за квадрат суем куда надо чиселки
		while (n > 0) {
			for (int i = 0; i < p; i++) {
				if (n > 0) {
					node_counter[i]++;
					n--;
				}
				else break;
			}
		}
		for (int i = 0; i < p; i++) {
			if (i != p - 1) {
				Ashki[i] = makeISF<double>(node_counter[i], nodes[i], interv[i].first, interv[i].second);
			}
			else {
				Ashki[i] = makeISF<double>(node_counter[i], nodes[i], interv[i].first, interv[i].second);
			}
		}
		//make sure that we do it right
		double overall_curr = 0;
		for (int i = 0; i < p; i++) {
			double curr = 0;
			cout <<"\n========\n" << Ashki[i] << "========\n";
			cout << i << "\n";
			for (int j = 0; j < node_counter[i]; j++) {
				curr += Ashki[i](j, 0);
			}
			overall_curr += curr;
		}
		double ans = 0;
		cout << "COMPUTING...\n\n\n";
		for (int i = 0; i < p; i++) {
			double XD = computeIntegral<double>(Ashki[i], nodes[i]);
			ans += XD;
		}
		return ans;
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

