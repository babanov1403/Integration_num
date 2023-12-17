#include <iostream>
#include <vector>

#include "matrix_class.h"
#include "polynom_class.h"
#include "integration_methods.h"


using namespace std;
using namespace biv;
//vpairdd = vector<pair<double, double>>
//строит равномерное разбиение интервала интегрирование на p частей
vpairdd evenNet(int p) {
	vector<pair<double, double>> res(p);
	double h = (var_b::b - var_b::a) / p;
	double curr = var_b::a;
	for (int i = 0; i < p; i++) {
		res[i] = make_pair(curr, curr + h);
		curr += h;
	}
	for (int i = 0; i < p; i++) {
		cout << res[i].first << ' ' << res[i].second << "\n";
	}
	return res;
}

int computeOptimum(double s1, double s2, double epsilon, double L, int p, int m){
	double H = (var_b::b - var_b::a) / p;
	double H_0 = H * pow(epsilon * (1 - pow(L, -1.0 * m)) / abs(s1 - s2), 1.0 / m);
	int opt = (var_b::b - var_b::a) / H_0;
	return opt;
}
//при заданной точности оптимальное число разбиений такое:
int main() {
	//процесс Эйткена
	
	constexpr double l = 0.5;
	constexpr int p = 5;
	constexpr int N = 50;
	double Q = (makeCompoundSF(N, p, evenNet) - makeCompoundSF(N, l * p, evenNet)) / (makeCompoundSF(N, l * p, evenNet) - makeCompoundSF(N, l * l * p, evenNet));
	int m = abs(log(abs(Q)) / log(l));
	cout << m; // 9
	
	
	//правило Рунге
	constexpr double epsilon = 1e-2;
	constexpr int L = 2;
	double s1 = makeCompoundSF(N, p, evenNet);
	double s2 = makeCompoundSF(N, L*p, evenNet);
	cout << "=============\n";
	cout << abs(s1 - s2) << '\n';
	cout << s1 << " " << s2 << '\n';
	int p_nice = computeOptimum(s1, s2, epsilon, L, p, m);
	cout << makeCompoundSF(N, p_nice, evenNet) << '\n';
	cout << m << '\n';
	cout << p_nice;
}