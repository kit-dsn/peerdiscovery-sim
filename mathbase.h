#pragma once

#include <vector>
#include <assert.h>
#include <math.h>

using namespace std;


double calcSquareError(const vector<double>& v1, const vector<double>& v2) {
	size_t minsize = min(v1.size(), v2.size());
	double errorsum = 0;
	for (size_t i = 0; i < minsize; ++i) {
		errorsum += (v1[i]-v2[i]) * (v1[i]-v2[i]);
	}
	return errorsum;
}



void convolve(vector<double>& conv, const vector<double>& x, const vector<double>& y) {
	size_t n = x.size();
	size_t m = y.size();

	conv.resize(n+m-1);
	fill(conv.begin(), conv.end(), 0);

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < m; ++j) {
			conv[i+j] += x[i] * y[j];
		}
	}

	return;
}

vector<double> convolve(const vector<double>& x, const vector<double>& y) {

	vector<double> r;
	convolve(r, x, y);
	return r;
}

vector<double> convolve(const vector<double>& x, int convpower = 2) {
	if (convpower < 2) {
		return x;
	}

	auto convdist = convolve(x, x);
	convpower -= 2;
	while (convpower > 0) {	// TODO this could be done logarithmically if larger powers are to be computed...
		convdist = convolve(x, convdist);
		--convpower;
	}
	return convdist;
}

template<typename T>
void getCdf(const vector<T>& pdf, vector<T>& cdf) {
	cdf.resize(pdf.size());
	if (pdf.size() == 0)
		return;

	cdf[0] = pdf[0];
	for (size_t i = 1; i < cdf.size(); ++i) {
		cdf[i] = cdf[i-1] + pdf[i];
	}

	return;
}

template<typename T>
vector<T> getCdf(const vector<T>& pdf) {
	vector<T> cdf;
	getCdf(pdf, cdf);
	return cdf;
}

void multVector(vector<double>& pdf, const double mult) {
	for (double& d : pdf) {
	    d *= mult;
	}
}

template<typename T>
vector<double> multVector(const vector<T>& pdf, const double mult) {
	vector<double> newpdf(pdf.size());
	for (size_t i = 0; i < newpdf.size(); ++i) {
	    newpdf[i] = pdf[i] * mult;
	}
	return newpdf;
}
void scalePdf(vector<double>& pdf) {
	double sum = 0;
	for (const double d : pdf) {
	    sum += d;
	}

	multVector(pdf, 1.0/sum);
}

template<typename T>
vector<double> scalePdf(const vector<T>& pdf) {
	vector<double> scaledpdf(pdf.size());
	double sum = accumulate(pdf.begin(), pdf.end(), 0.0);

	for (size_t i = 0; i < pdf.size(); ++i) {
		scaledpdf[i] = pdf[i] / sum;
	}
	return scaledpdf;
}

template<typename T>
vector<double> scaleCdf(const vector<T>& cdf) {
	vector<double> scaledcdf(cdf.size());
	double sum = cdf.back();

	for (size_t i = 0; i < cdf.size(); ++i) {
		scaledcdf[i] = cdf[i] / sum;
	}
	return scaledcdf;
}

void scaleCdf(vector<double>& cdf) {
	double sum = cdf.back();

	for (double& d : cdf) {
	    d /= sum;
	}
}

void sumpdf(vector<double>& sum, const vector<double>& add) {
	if (add.size() > sum.size()) {
		sum.resize(add.size());
	}
	for (size_t i = 0; i < sum.size(); ++i) {
		sum[i] += add[i];
	}
	// don't scale here as multiple additions then favor later added pdfs
}

/*
 * extends distributions to a specified length by appending 'v' values
 * does NOT modify empty distributions
 */
void extendDist (vector<double>& dist, size_t size, double v = 0) {
	if (dist.size() > 0 && dist.size() < size) {
		dist.resize(size, v);
	}
}

void equalizeDistLength (vector<double>& distA, vector<double>& distB, double v = 0) {
	if (distA.size() < distB.size()) {
		extendDist(distA, distB.size(), v);
	} else if (distA.size() > distB.size()) {
		extendDist(distB, distA.size(), v);
	}
	assert(distA.size() == distB.size());
}



/*
 * generates a histogram out of a vector of values.
 */
template<typename T>
void genHist(vector<double>& hist, const vector<T>& values, const T min, const T max,
		const size_t nbins, bool relative = true, int strechby = 1) {

	hist.resize(nbins*strechby);

	if (values.size() == 0)
		return;

	T stepsize = (max - min) / nbins;
	for (auto v : values) {
		if (v >= min && v <= max) {
			size_t pos = (v-min)/stepsize*strechby;
			if (pos >= hist.size()) {
				pos = hist.size() - 1;
			}
			assert(pos >= 0 && pos < hist.size());
			hist[pos]++;
		}
	}

	if (relative) {
		for (double& v : hist) {
			v /= values.size()*strechby;
		}
	}
	if (strechby == 1) {
		return;
	} else {
		double v = -999;
		for(size_t i = 0; i < hist.size(); ++i) {
			if (i % strechby == 0) {
				v = hist[i];
			} else {
				hist[i] = v;
			}

		}
	}
}


template<typename T>
vector<double> genHist(const vector<T>& values, const T min, const T max, const size_t nbins, bool relative = true, int strechby = 1) {
	vector<double> hist;
	genHist<T>(hist, values, min, max, nbins, relative, strechby);
	return hist;
}

template<typename T>
double mean(const vector<T>& values) {
	double offset = 0;
	return accumulate(values.begin(), values.end(), offset) / values.size();
}

template<typename T>
double expected(const vector<T>& values) {
	double avg = 0;
	for (size_t i = 0; i < values.size(); ++i) {
		avg += values[i] * i;
	}
	return avg;
}

template<typename T>
vector<pair<int, pair<double, double>>> genLogHist(const vector<T>& values, double logbase, T max = 0) {
	vector<pair<int, pair<double, double>>> hist;

	if (max == 0) {
		for (const auto v : values) {
			if (v > max)
				max = v;
		}
	}

	size_t nbins = (size_t) (log(max)/log(logbase) + 1);
	hist.resize(nbins);
	for (const auto v : values) {
		size_t pos = (size_t) (log(v)/log(logbase));
		if (pos >= hist.size()) {
			cerr << "value too large " << v << " -> " << log(v)/log(logbase) << " -> " << pos << endl;
			continue;
		}
		hist[pos].first++;
	}

	for (size_t i = 0; i < nbins; ++i) {
		hist[i].second.first = pow(logbase, i); // this value is part of the bin [
		hist[i].second.second = pow(logbase, i+1); // this value is NOT part of the bin )
	}

	return hist;
}



template<typename T>
double variance(const vector<T>& values) {
	T m = mean<T>(values);
	T sum = 0;
	for (long double v : values) {
		sum += (v - m) * (v - m);
	}
	return sum / (double)values.size();
}

template<typename T>
double stdDev(const vector<T>& values) {
	return sqrt(variance(values));
}

template<typename T>
T median(vector<T> values) {
	sort(values.begin(), values.end());
	return values[values.size()/2];
}

template<typename T>
double confidence(const vector<T>& values) {
	return 1.96 * stdDev<T>(values) / sqrt(values.size());
}

/*
 * takes two histograms as arguments that represent PDFs
 * calculates the probability that when selecting a single value from the distribution of A
 * and a single value from the distribution of b, the value taken from a is smaller than the value
 * of b.
 *
 * \sum_t P(A=t) * (P(B > t) + 0.5*P(B = t))
 *
 * see: https://stats.stackexchange.com/questions/102778/correlations-between-continuous-and-categorical-nominal-variables
 */
template<typename T>
double calcRankCorrelation(const vector<T>& a, const vector<T>& b) {
	auto pdfA = scalePdf<T>(a);
	auto pdfB = scalePdf<T>(b);
	auto cdfB = getCdf<double>(pdfB);


	double sum = 0;
	for (size_t t = 0; t < a.size(); ++t) {
		sum += pdfA[t] * ((1-cdfB[t]) + 0.5*pdfB[t]);

	}

	return sum;
}

