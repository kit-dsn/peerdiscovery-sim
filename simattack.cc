#include <fstream>
#include <iostream>

#include "mathbase.h"
#include "random.h"
#include "gp.h"
#include "json.hpp"


using json = nlohmann::json;
using namespace std;


struct simresult {
	vector<int> tp;
	vector<int> tn;
	vector<int> fp;
	vector<int> fn;

	simresult(const size_t n) : tp(n), tn(n), fp(n), fn(n) {}

};

simresult simAttacker(vector<double>& pdfA, vector<double>& pdfB, const double aprioriA,
		const double pSendIPconnected, const double pSendIPnotconnected,
		const size_t nObservations, const size_t nTrials) {

	equalizeDistLength(pdfA, pdfB);

	simresult sr(nObservations);

	auto cdfA = getCdf(pdfA);
	auto cdfB = getCdf(pdfB);



	/*
	 * one trial is one possible, randomly selected edge
	 */
	for (size_t i = 0; i < nTrials; ++i) {

		// that edge either exists or not, according to the apriori probability
		// the existence of the edge remains throughout all observations
		bool connected = i < nTrials * aprioriA;
		vector<double>& selectedCdf = connected ? cdfA : cdfB;



		/*
		 * the attacker starts with zero observations and hence only the apriori probabilities
		 * each observation modifies the likelihoods
		 */
		double likelihoodA = aprioriA;
		double likelihoodB = 1 - aprioriA;

		// one observation is the reception of one addr message from one of the
		// two ends of the edge.
		// an observation can be either "no ip" in case the other ip was not sent
		// or an age of the ip
		for (size_t ob = 0; ob < nObservations; ++ob) {

			bool sendip;

			if (connected) {
				sendip = GetRandProb() < pSendIPconnected;
			} else {
				sendip = GetRandProb() < pSendIPnotconnected;
			}

			if (sendip) {
				likelihoodA *= pSendIPconnected;	// P(sendIP | connected)
				likelihoodB *= pSendIPnotconnected;	// P(sendIP | !connected)

				size_t age;
				age = getInverseElement(selectedCdf);
				assert(age < pdfA.size());
				likelihoodA *= pdfA[age];
				likelihoodB *= pdfB[age];

			} else {
				likelihoodA *= 1 - pSendIPconnected;	// P(!sendIP | connected)
				likelihoodB *= 1 - pSendIPnotconnected;	// P(!sendIP | !connected)
			}


			bool detectEdge = (likelihoodA > likelihoodB);

			if (connected && detectEdge) ++sr.tp[ob];
			if (connected && !detectEdge) ++sr.fn[ob];
			if (!connected && detectEdge) ++sr.fp[ob];
			if (!connected && !detectEdge) ++sr.tn[ob];
		}

	}
	return sr;

}



int main(int argc, char** argv) {


	if (argc < 6) {
		cout << "USAGE: ./simattack <simresult.json> <ts> <nObservations> <nTrials> <seed>" << endl;
		return 0;
	}

	std::ifstream ifs(argv[1]);
	int ts = atoi(argv[2]);
	size_t nObservations = atoi(argv[3]);
	int nTrials = atoi(argv[4]);
	SetSeed(atoi(argv[5]));

	json j;
	ifs >> j;

	for (const auto& e : j) {
		if (e["ts"] != ts) {
			continue;
		}

		vector<double> statAgeConnected = e["statAgeConnected"];
		vector<double> statAgeNotConnected = e["statAgeNotConnected"];
		double Pconn = e["Pconn"];
		double pSendIPconnected = e["pSendIPconnected"];
		double pSendIPnotconnected = e["pSendIPnotconnected"];
		scalePdf(statAgeConnected);
		scalePdf(statAgeNotConnected);




		simresult sr = simAttacker(statAgeConnected, statAgeNotConnected, Pconn,
				pSendIPconnected, pSendIPnotconnected, nObservations, nTrials);

		for (size_t ob = 0; ob < nObservations; ++ob) {
			cout << ob+1 << "\t" << e["maxFakeAgeConnection"] << "\t" << e["maxAddrSend"]
					<< "\t" << sr.tp[ob]/(double)(sr.tp[ob]+sr.fp[ob])
					<< "\t" << sr.tp[ob]/(double)(sr.tp[ob]+sr.fn[ob])<< endl;
		}

	}

}
