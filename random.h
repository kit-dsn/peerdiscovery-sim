#pragma once

#include <random>
#include <stdlib.h>
#include <cmath>
#include <assert.h>
#include <climits>
#include <iostream>

#include <vector>
#include <algorithm>

namespace {

	std::mt19937 trand_gen(1);

	void SetSeed(int seed) {
		trand_gen.seed(seed);
	}

	// returns values SMALLER than max! (as done in bitcoin)
	bool GetRandBool() {
		return trand_gen() % 2;
	};
	int GetRandInt(int max) {
		if (max == 0) return 0;
		return trand_gen() % max;
	};
	int GetRandInt(int min, int max) {
		return min + (trand_gen() % (max-min));
	};
	int GetRandIntMax() {
		return trand_gen();
	};

	double GetRandProb() {	// value between 0 and 1
		return trand_gen() / (double)trand_gen.max();
	};

	double GetRandDouble(double min, double max) {	// value between min and max
		return min + (trand_gen() / (double)trand_gen.max())*(max-min);
	};

	int GetUniformInt(int min, int max) {
		return GetRandInt(max-min)+min;
	};

	double GetExpDist(double mean, double bound = 0)
	{
	  while (1)
	    {
	      double r = -mean*std::log (GetRandProb());
	      if (bound == 0 || r < bound) {
	          return r;
	      }
	      // otherwise, try again
	    }
	}

	double GetPowerDist(double alpha, double xmin) {
		return xmin*pow(1.0-GetRandProb(), 1.0/(1.0-alpha));
	}



	struct tr_normdist_t {
		bool tr_init = false;

		std::mt19937 gen;	// mersenne twister
		std::normal_distribution<> normdist;

		int lowerlimit;
		int upperlimit;

		int getVal() {
			int r;
			do {
				r = (int)normdist(gen);
			} while (r < lowerlimit || r > upperlimit);

			return r;
		}

		tr_normdist_t(double mean, double sd, int llimit = 0, int ulimit = INT_MAX, int seed = 0)
		: gen(seed), normdist(mean, sd), lowerlimit(llimit), upperlimit(ulimit) {};
	};
	struct tr_normdist_t *tr_normdist = NULL;

	void InitNormDist(int mean, int sd, int seed) {
		assert(tr_normdist == NULL);
		tr_normdist = new tr_normdist_t(mean, sd, 0, INT_MAX, seed);
	}

	int GetNormDist() {
		return tr_normdist->getVal();
	}


	struct tr_possiondist_t {
		bool tr_init = false;

		std::mt19937 gen;	// mersenne twister
		std::poisson_distribution<int> poissondist;

		int upperlimit;

		int getVal() {
			int r;
			do {
				r = poissondist(gen);
			} while (r > upperlimit);

			return r;
		}

		tr_possiondist_t(double mean, int ulimit = INT_MAX, int seed = 0)
		: gen(seed), poissondist(mean), upperlimit(ulimit) {};
	};
	struct tr_possiondist_t *tr_poissondist = NULL;

	void InitPoissonDist(int mean, int seed) {
		if(tr_poissondist != NULL) {
			delete tr_poissondist;
		}
		tr_poissondist = new tr_possiondist_t(mean, INT_MAX, seed);
	}

	int GetPoissonDist() {
		return tr_poissondist->getVal();
	}


	int getInverseElement(const std::vector<double>& cdf) {
		double p = GetRandProb();

		auto low = lower_bound(cdf.begin(), cdf.end(), p);
		int d = low - cdf.begin();

		if (d >= cdf.size()) {
			d = cdf.size() - 1;
		}
		return d;
	}

}
