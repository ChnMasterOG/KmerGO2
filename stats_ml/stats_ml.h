/*
  Project Name	: KmerGO
  Version		: 2.0.0
  Author		: Qi Chen
  Date			: 2020-07-27
*/

#ifndef _STATS_ML_H
#define _STATS_ML_H

struct Logistic_Data {
	std::vector<float> features;
	int cls;
	Logistic_Data(std::vector<float> f, int c) :features(f), cls(c) {
	}
};

struct Logistic_Param {
	std::vector<double> w;
	double d;
	Logistic_Param(std::vector<double> w1, double d1) :w(w1), d(d1) {};
	Logistic_Param() :w(std::vector<double>()), d(0.0) {}
};

class Logistic {
public:
	Logistic(std::vector<Logistic_Data>& dataSet);
	void logisticRegression();
	bool samewb(const Logistic_Param& tparam, const Logistic_Param& param, double delta);
	void gradient(double lam);
	double logiFun(const Logistic_Param& p, const Logistic_Data& d);
	double Lw(Logistic_Param p);
	double innerWX(const Logistic_Param& p, const Logistic_Data& data);
	void predictClass(std::vector<Logistic_Data>& testDataSet);
private:
	std::vector<Logistic_Data> dataSet;
	Logistic_Param param;
};

void getrank(std::vector<float> group, std::vector<float>& grouprank);
float ranksums(std::vector<float>& group1, std::vector<float>& group2);
float spearman_correlation_coefficient(std::vector<float> group1, std::vector<float> group2);

#endif
