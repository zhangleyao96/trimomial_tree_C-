//this file is to 
//1.calculate american option price with the average binomial tree
//2.calculate american option price with different N,by Trinomial Black¨CScholes,
//  and Trinomial Black¨CScholes with Richardson Extrapolation 
//3.conduct a comparison between these method in the aspects of delta,gamma and theta.




#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>

const double M_SQRT1_2 = sqrt(0.5);
double normalCFD(double value)
{
	return 0.5 * erfc(-value * M_SQRT1_2);
}


class Trinomial_Op
{
private:
	long N;
	double r, q, sigma,S,K,T;
	double Value, Delta, Gamma, Theta;
	double delta_t;
	double u;
	double d;
	double pu;
	double pm;
	double pd;
	double disc;

public:
	Trinomial_Op(int _N, double _r, double _q, double _sigma, double _S,double _K,double _T) :N(_N), r(_r), q(_q), sigma(_sigma), S(_S),K(_K),T(_T) {}
	double BlackScholes_Put(double _r, double _q, double _sigma, double _S, double _K, double _T);
	double calculate_America();
	double calculate_Euro();
	double calculate_Variance_Reduction();
};


double Trinomial_Op::BlackScholes_Put(double r, double q, double sigma, double S, double K, double T)
{
	double d1 = (std::log(S / K) + (r - q + sigma * sigma / 2)*T) / (sigma * sqrt(T));
	double d2 = (std::log(S / K) + (r - q - sigma * sigma / 2)*T) / (sigma * sqrt(T));
	return -S * exp(-q * T)*normalCFD(-d1) + K * exp(-r * T)*normalCFD(-d2);
}


double Trinomial_Op::calculate_America() {
	delta_t = T / static_cast<double>(N);
	u = exp(sigma*sqrt(3 * delta_t));
	d = exp(-sigma * sqrt(3 * delta_t));

	pu = 1.0 / 6.0 + (r - q - sigma * sigma / 2)*sqrt(delta_t / (12.0*sigma*sigma));
	pm = 2.0 / 3.0;
	pd = 1.0 / 6.0 - (r - q - sigma * sigma / 2)*sqrt(delta_t / (12.0*sigma*sigma));

	disc = exp(-r * delta_t);

	std::vector<double> v(2*N+1);

	for (long i = 0; i < 2*N+1; i++) {
		v[i] = std::max(K - S * std::pow(u, N - i),0.0);
	}
	for (long i = N - 1; i > 1; --i)
	{
		for (long j = 0; j <= 2*i; j++)
		{
		v[j] = std::max(disc * (pu*v[j] + pm * v[j + 1] + pd * v[j + 2]),std::max(K-S* pow(u, i - j),0.0));
		}
		//std::cout << i << std::endl;
	}
	double v20 = v[0];
	double v22 = v[2];
	double v24 = v[4];
	for (long j = 0; j <= 2; j++)
	{
		v[j] = std::max(disc * (pu*v[j] + pm * v[j + 1] + pd * v[j + 2]), std::max(K - S * pow(u, 1 - j),0.0));
	}
	double v10 = v[0];
	double v12 = v[2];
	double v11 = v[1];
	v[0] = std::max(disc * (pu*v[0] + pm * v[1] + pd * v[2]), std::max(K - S,0.0));
	
	Value = v[0];
	Delta = -(v10 - v12) / (S*d-S*u);
	Gamma = ((v20-v22)/(S*u*u-S)-(v22-v24)/(S-S*d*d))/(-S*d + S*u);
	Theta = (v11-Value) / (delta_t);


	return Value;
};


double Trinomial_Op::calculate_Euro() {
	delta_t = T / static_cast<double>(N);
	u = exp(sigma*sqrt(3 * delta_t));
	d = exp(-sigma * sqrt(3 * delta_t));

	pu = 1.0 / 6.0 + (r - q - sigma * sigma / 2)*sqrt(delta_t / (12.0*sigma*sigma));
	pm = 2.0 / 3.0;
	pd = 1.0 / 6.0 - (r - q - sigma * sigma / 2)*sqrt(delta_t / (12.0*sigma*sigma));

	disc = exp(-r * delta_t);

	std::vector<double> v(2 * N + 1);

	for (long i = 0; i < 2 * N + 1; i++) {
		v[i] = std::max(K - S * std::pow(u, N - i), 0.0);
	}
	for (long i = N - 1; i > 1; --i)
	{
		for (long j = 0; j <= 2 * i; j++)
		{
			v[j] = disc * (pu*v[j] + pm * v[j + 1] + pd * v[j + 2]);
			//v[j] = std::max(disc * (pu*v[j] + pm * v[j + 1] + pd * v[j + 2]),std::max(K-S* pow(u, i - j),0.0));
		}
		//std::cout << i << std::endl;
	}
	double v20 = v[0];
	double v22 = v[2];
	double v24 = v[4];
	for (long j = 0; j <= 2; j++)
	{
		v[j] = disc * (pu*v[j] + pm * v[j + 1] + pd * v[j + 2]);
		//v[j] = std::max(disc * (pu*v[j] + pm * v[j + 1] + pd * v[j + 2]), std::max(K - S * pow(u, 1 - j),0.0));
	}
	double v10 = v[0];
	double v12 = v[2];
	double v11 = v[1];
	v[0] = disc * (pu*v[0] + pm * v[1] + pd * v[2]);
	//v[0] = std::max(disc * (pu*v[0] + pm * v[1] + pd * v[2]), std::max(K - S,0.0));

	Value = v[0];
	Delta = -(v10 - v12) / (S*d - S * u);
	Gamma = ((v20 - v22) / (S*u*u - S) - (v22 - v24) / (S - S * d*d)) / (-S * d + S * u);
	Theta = (v11 - Value) / (delta_t);


	return Value;
};


double Trinomial_Op::calculate_Variance_Reduction() {
	delta_t = T / static_cast<double>(N);
	u = exp(sigma*sqrt(3 * delta_t));
	d = exp(-sigma * sqrt(3 * delta_t));

	pu = 1.0 / 6.0 + (r - q - sigma * sigma / 2)*sqrt(delta_t / (12.0*sigma*sigma));
	pm = 2.0 / 3.0;
	pd = 1.0 / 6.0 - (r - q - sigma * sigma / 2)*sqrt(delta_t / (12.0*sigma*sigma));

	disc = exp(-r * delta_t);

	std::vector<double> v(2 * N + 1);
	std::vector<double> v_euro(2 * N + 1);

	for (long i = 0; i < 2 * N + 1; i++) {
		v_euro[i] = std::max(K - S * std::pow(u, N - i), 0.0);
		v[i] = std::max(K - S * std::pow(u, N - i), 0.0);
	}
	for (long i = N - 1; i > 1; --i)
	{
		for (long j = 0; j <= 2 * i; j++)
		{
			v_euro[j] = disc * (pu*v_euro[j] + pm * v_euro[j + 1] + pd * v_euro[j + 2]);
			v[j] = std::max(disc * (pu*v[j] + pm * v[j + 1] + pd * v[j + 2]),std::max(K-S* pow(u, i - j),0.0));
		}
		//std::cout << i << std::endl;
	}
	double v20 = v[0]-v_euro[0]+ BlackScholes_Put(r, q, sigma, S*u*u, K, T*(N - 2) / N);
	double v22 = v[2]-v_euro[2]+ BlackScholes_Put(r, q, sigma, S, K, T*(N - 2) / N);
	double v24 = v[4]-v_euro[4]+ BlackScholes_Put(r, q, sigma, S*d*d, K, T*(N - 2) / N);

	for (long j = 0; j <= 2; j++)
	{
		v_euro[j] = disc * (pu*v_euro[j] + pm * v_euro[j + 1] + pd * v_euro[j + 2]);
		v[j] = std::max(disc * (pu*v[j] + pm * v[j + 1] + pd * v[j + 2]), std::max(K - S * pow(u, 1 - j),0.0));
	}

	double v10 = v[0] - v_euro[0]+ BlackScholes_Put(r, q, sigma, S*u, K, T*(N-1)/N);
	double v12 = v[2] - v_euro[2]+ BlackScholes_Put(r, q, sigma, S*d, K, T*(N - 1) / N);
	double v11 = v[1] - v_euro[1]+ BlackScholes_Put(r, q, sigma, S, K, T*(N - 1) / N);

	v_euro[0] = disc * (pu*v_euro[0] + pm * v_euro[1] + pd * v_euro[2]);
	v[0] = std::max(disc * (pu*v[0] + pm * v[1] + pd * v[2]), std::max(K - S,0.0));

	Value = v[0]-v_euro[0]+BlackScholes_Put(r, q, sigma, S, K, T);
	Delta = -(v10 - v12) / (S*d - S * u);
	Gamma = ((v20 - v22) / (S*u*u - S) - (v22 - v24) / (S - S * d*d)) / (-S * d + S * u);
	Theta = (v11 - Value) / (delta_t);
	return Value;
};






int main()
{
	
	Trinomial_Op Option1(1000, 0.03, 0.005, 0.285, 41.0,39.0,1);
	std::cout << std::setprecision(9)<< Option1.calculate_America() << std::endl;
	//std::cout << Option1.calculate_Euro() << std::endl;
	std::cout << Option1.calculate_Variance_Reduction() << std::endl;

	return 0;
}
