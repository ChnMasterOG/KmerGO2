/*
  Project Name	: KmerGO
  Version		: 2.0.0
  Author		: Qi Chen
  Date			: 2020-07-27
*/

#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>
#include<cmath>
#include"../stats_ml/stats_ml.h"

Logistic::Logistic(std::vector<Logistic_Data>& dataSet) {
	std::vector<double> pw(dataSet[0].features.size(), 0.0);
	Logistic_Param pt(pw, 0.0);
	param = pt;

};

void Logistic::logisticRegression() {
	double lamda = 0.1;		//step size
	double delta = 0.001;	//threshold of breaking iteration
	double objLw = Lw(param);
	Logistic_Param tpa(param.w, param.d);
	gradient(lamda);
	double newObjLw = Lw(param);
	int iter = 0;
	while (fabs(newObjLw - objLw) > delta || !samewb(tpa, param, delta)) {
		objLw = newObjLw;
		tpa = Logistic_Param(param.w, param.d);
		gradient(lamda);
		newObjLw = Lw(param);
		++iter;
	}
}

bool Logistic::samewb(const Logistic_Param& tparam, const Logistic_Param& param, double delta) {
	for (uint32_t i = 0; i < tparam.w.size(); i++) {
		if (fabs(tparam.w[i] - param.w[i]) > delta) {
			return false;
		}
	}
	if (fabs(tparam.d - param.d) > delta) {
		return false;
	}
	return true;
}

void Logistic::gradient(double lam) {
	for (uint32_t i = 0; i < param.w.size(); i++) {
		double tmp = 0.0L;
		for (uint32_t j = 0; j < dataSet.size(); j++) {
			tmp += (dataSet[j].cls - logiFun(param, dataSet[j])) * dataSet[j].features[i] * lam;
		}
		param.w[i] += (tmp);
	}
	double tmp = 0.0L;
	for (uint32_t j = 0; j < dataSet.size(); j++) {
		tmp += (dataSet[j].cls - logiFun(param, dataSet[j])) * lam;
	}
	param.d += tmp;

}
double Logistic::logiFun(const Logistic_Param& p, const Logistic_Data& d) {
	double inner = innerWX(p, d);
	double le = exp(inner) / (1 + exp(inner));
	return le;
}

double Logistic::Lw(Logistic_Param p) {
	double l = 0.0L;
	for (uint32_t i = 0; i < dataSet.size(); i++) {
		double inner = innerWX(p, dataSet[i]);
		l += (dataSet[i].cls * inner - (log10(1 + exp(inner))));
	}
	return l;
}

double Logistic::innerWX(const Logistic_Param& p, const Logistic_Data& data) {
	if (p.w.size() != data.features.size()) {
		exit(0);
	}
	double innerP = 0.0L;
	for (uint32_t i = 0; i < p.w.size(); i++) {
		innerP += (p.w[i] * data.features[i]);
	}
	innerP += p.d;
	return innerP;
}

void Logistic::predictClass(std::vector<Logistic_Data>& testDataSet) {
	for (uint32_t i = 0; i < testDataSet.size(); i++) {
		double py1 = 0.0L;
		double py0 = 0.0L;
		double inner = innerWX(param, testDataSet[i]);
		py1 = exp(inner) / (1 + exp(inner));
		py0 = 1 - py1;
		if (py1 >= py0) testDataSet[i].cls = 1;	// class 1
		else testDataSet[i].cls = 0;	// class 0
	}
}
