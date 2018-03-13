#include "json.hpp"

#include "random.h"
#include "mathbase.h"
#include "gp.h"

#include <sstream>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <map>
#include <algorithm>
#include <iterator>
#include <unistd.h>

using namespace std;
using json = nlohmann::json;

//#define LOG_DEBUG_ENABLED
#ifdef LOG_DEBUG_ENABLED

#define LOG(msg) \
  do {                      \
    cout << msg;    \
  } while (false)

#else

#define LOG(msg) \
  do {                      \
  } while (false)
#endif

typedef size_t peerid_t;
typedef uint32_t simtime_t;

enum etype {
	ADD, DEL, SENDADDR, ADDRSTAT
};

struct event {
	peerid_t receiver;
	etype type;
	peerid_t data;
};

struct peer {
	bool online;
	set<peerid_t> incoming;
	set<peerid_t> outgoing;
	unordered_map<peerid_t,simtime_t> addrlist;

	peer() : online (true) {

	}
};


const simtime_t oneSecond = 1;	// resolution: second
const simtime_t oneMinute = 60 * oneSecond;
const simtime_t oneHour = 60 * oneMinute;
const simtime_t oneDay = 24 * oneHour;




simtime_t maxFakeAgeConnection = oneHour;
simtime_t maxAddrTotalAge = 24 * oneHour;
simtime_t meanAddrSendInterval = oneHour;
size_t maxAddrSend = 10;
int seed = 42;

const simtime_t meanSessionLength = 3*oneHour;
const char SLdist = 'u';

simtime_t maxSimTime = 7 * oneDay;
const simtime_t startSimTSDefault = oneHour;
const simtime_t statInterval = 3 * oneHour;
const simtime_t warmupTime = 24 * oneHour;

const size_t nObservations = 20;
const size_t nTrials = 10000;

const size_t nOutgoingConnections = 8;
const size_t maxIncomingConnections = 125 - nOutgoingConnections;
const size_t numPeersDefault = 1000;


multimap<peerid_t,event> eventMap;
vector<peer> peers;

// peer --> sorted list of join/leave events with timestamp
// the size_t pair is the index of the event that is next executed
vector<pair<vector<pair<simtime_t, etype>>, size_t>> churnmap;


vector<size_t> statAgeConnected;
vector<size_t> statAgeNotConnected;

size_t statAddrConnectionNew = 0;
size_t statAddrConnectionExisting = 0;
size_t statAddrMsgNew = 0;
size_t statAddrMsgUpdateNewer = 0;
size_t statAddrMsgUpdateOlder = 0;
size_t statAddrUpdate = 0;
size_t statAddrDelete = 0;
size_t statAddrOffline = 0;

size_t statEventJoin = 0;
size_t statEventLeave = 0;
size_t statEventAddr = 0;
size_t statEventAddrScheduled = 0;
size_t statEventAddrCancelled = 0;

const simtime_t statCompress = oneSecond;

string jsonfname;
vector<json> jsonvec;

//http://stackoverflow.com/questions/372484/how-do-i-programmatically-check-memory-use-in-a-fairly-portable-way-c-c
size_t memory_used () {
	size_t size = 0;
	FILE *file = fopen("/proc/self/statm", "r");
	if (file) {
		unsigned int vm = 0;
		int rc = fscanf (file, "%ul", &vm);  // Just need the first num: vm size
		if (rc > 0) {
			size = (size_t)vm * getpagesize();
		}
		fclose (file);
	}
	return size;
}


simtime_t getSLDist(const double mean) {
	double camelShare = 0.99;
	double camelLow = mean/2.0;

	switch (SLdist) {
	case 'f':
		return mean;
	case 'u':
		return GetRandInt(mean*2);
	case 'e':
		return (int32_t)GetExpDist(mean, INT_MAX);
	case 'p': {
		// this just sets the MEDIAN and not the mean.
		// Hence, it's not linear anymore, i.e. getSLDist(a*b) != a * getSLDist(b)
		double alpha = 2.0;
		double x0 = mean / (pow(2.0, -1.0/(1.0-alpha)));
		//cout << x0 << endl;
		return (int32_t)GetPowerDist(alpha, x0);
	}
	case 'c':	// camel distribution
		return (int32_t)(GetRandProb() < camelShare ? camelLow : (mean-camelShare*camelLow)/(1.0-camelShare));
	}

	return 0;
}

simtime_t getSLDiff(const double mean) {
	double camelShare = 0.99;
	double camelLow = mean/2.0;

	switch (SLdist) {
	case 'f':
		return GetRandInt(mean);
	case 'u':
		return abs(GetRandInt(mean*2)-GetRandInt(mean*2));
	case 'e':
		return abs((int32_t)(GetExpDist(mean, INT_MAX)-GetExpDist(mean, INT_MAX)));
	case 'p': {
		double alpha = 2.0;
		double x0 = mean / (pow(2, -1.0/(1.0-alpha)));
		return abs((int32_t)(GetPowerDist(alpha, x0)-GetPowerDist(alpha, x0)));
	}
	case 'c':
		return GetRandProb() < camelShare ? GetRandInt((int)camelLow)
				: GetRandInt((int)((mean-camelShare*camelLow)/(1.0-camelShare)));
	}

	return 0;
}


simtime_t getRandomAddrSendInterval() {
	return GetRandInt(meanAddrSendInterval);
	//return meanAddrSendInterval-1;
}

simtime_t getAddrSendInterval() {
	return meanAddrSendInterval;
}

simtime_t getAddrFakeAge() {
	/*return 0;

	int age;
	do {
		age = GetNormDist();
	} while (age < 0 || age >= maxFakeAgeConnection);
	return age;
	*/
	return GetRandInt(maxFakeAgeConnection);
}





/*
 * checks whether peer p will be online (without interruption) from tsnow until tsscheduled
 */
bool checkPeerOnline(const peerid_t p, const simtime_t tsnow, const simtime_t tsscheduled) {
	assert(p < churnmap.size());
	assert(peers[p].online);

	if (churnmap[p].second < churnmap[p].first.size()) {
		if (churnmap[p].first[churnmap[p].second].first > tsscheduled) {
			// next scheduled event is later than event to be scheduled
			return true;
		} else {
			return false;
		}
	} else {
		// there is no more event
		return true;
	}

}

void scheduleSendaddrEvent(const simtime_t ts, const peerid_t from, const peerid_t to, const bool randomized) {
	simtime_t newts = randomized ? ts + getRandomAddrSendInterval() : ts + getAddrSendInterval();
	if (!(checkPeerOnline(from, ts, newts) && checkPeerOnline(to, ts, newts))) {
		return;
	}
	event e;
	e.type = SENDADDR;
	e.receiver = from;
	e.data = to;
	eventMap.insert(make_pair(newts, e));
	++statEventAddrScheduled;
}


/*
 * handles both incoming and outgoing connections
 */
void handleNewConnection(const simtime_t ts, const peerid_t client, const peerid_t server, bool directed = false) {
	if (!directed) {
		handleNewConnection(ts, client, server, true);
		handleNewConnection(ts, server, client, true);
	} else {
		assert(client < peers.size());
		// add connection to addrlist
		auto it = peers[client].addrlist.find(server);
		if (it == peers[client].addrlist.end()) {
			++statAddrConnectionNew;
			peers[client].addrlist[server] = ts - getAddrFakeAge();
		} else {
			++statAddrConnectionExisting;
			it->second = ts - getAddrFakeAge();
		}

		// schedule sendaddr event
		scheduleSendaddrEvent(ts, client, server, true);
	}
}



void establishOutgoing(const simtime_t ts, const peerid_t client, bool force = false) {
	assert(client < peers.size());
	while (peers[client].outgoing.size() < nOutgoingConnections) {
		const peerid_t server = GetRandInt(peers.size());
		if (server == client) {	// connection to itself
			continue;
		}
		if (peers[server].online == false) {	// remote offline
			continue;
		}
		if (peers[client].incoming.count(server) > 0 || peers[client].outgoing.count(server) > 0) {	// already connected
			continue;
		}
		if (!force && peers[server].incoming.size() >= maxIncomingConnections) {	// no incoming connection slots at remote
			continue;
		}

		peers[client].outgoing.insert(server);
		peers[server].incoming.insert(client);
		//cout << "connection " << client << "\t" << server << endl;
		handleNewConnection(ts, client, server);
	}
}


void genAdj(const simtime_t startSimTS) {
	cout << "genAdj: startSimTS = " << startSimTS << " peers.size() = " << peers.size() << endl;
	for (peerid_t i = 0; i < peers.size() ; ++i) {
		establishOutgoing(startSimTS, i);
	}
}

void createChurn(const simtime_t startSimTS) {
	// create DEL events for existing peers
	peerid_t nPeers = peers.size();
	for (peerid_t i = 0; i < nPeers; ++i) {
		simtime_t ts = startSimTS + getSLDiff(meanSessionLength);
		event enew;
		enew.receiver = i;
		enew.type = DEL;

		eventMap.insert(make_pair(ts, enew));

	    if (churnmap.size() <= i) {
	    	churnmap.resize(i+1);
	    }
	    churnmap[i].first.emplace_back(ts, DEL);
	}

	simtime_t tjoin = startSimTS;
	while(tjoin < maxSimTime) {
		tjoin += getSLDist(meanSessionLength)/numPeersDefault;

		event e;
		e.type = ADD;
		e.receiver = nPeers;
		eventMap.insert(make_pair(tjoin, e));

		e.type = DEL;
		e.receiver = nPeers;
		simtime_t tleave = tjoin + getSLDist(meanSessionLength);
		eventMap.insert(make_pair(tleave, e));

	    if (churnmap.size() <= nPeers) {
	    	churnmap.resize(nPeers+1);
	    }
	    churnmap[nPeers].first.emplace_back(tjoin, ADD);
	    churnmap[nPeers].first.emplace_back(tleave, DEL);

		++nPeers;
	}

}



void handleReceiveAddr(const simtime_t ts, const peerid_t sender, const peerid_t receiver,
		const pair<const peerid_t,simtime_t>& addr) {

	auto it = peers[receiver].addrlist.find(addr.first);
	if (it == peers[receiver].addrlist.end()) {
		++statAddrMsgNew;
		peers[receiver].addrlist[addr.first] = addr.second;
	} else {
		if (it->second < addr.second) {	// received timestamp is newer
			it->second = addr.second;
			++statAddrMsgUpdateNewer;
		} else {
			++statAddrMsgUpdateOlder;
		}
	}

	if (ts > warmupTime) {
		size_t pos = (ts - addr.second)/statCompress;
		assert(pos < statAgeConnected.size());
		if (peers[sender].incoming.count(addr.first) > 0 || peers[sender].outgoing.count(addr.first) > 0) {
			++statAgeConnected[pos];	// online and connected
		} else if (peers[addr.first].online == true) {
			++statAgeNotConnected[pos];	// online but not connected
		} else {
			// offline --> ignore (attacker knows offline peers)
		}

	}

}


void updateConnectionTimestamps(peer& p, const simtime_t curts) {
	for (const auto& neighborlist : {p.incoming, p.outgoing}) {
		for (const auto& neighbor : neighborlist) {
			auto it = p.addrlist.find(neighbor);
			if (it == p.addrlist.end()) {
				cout << neighbor << endl;
				assert(false);
			}

			if (it->second < curts - maxFakeAgeConnection) {
				it->second = curts - getAddrFakeAge();
				++statAddrUpdate;
			}
		}
	}

	auto it = p.addrlist.begin();
	while (it != p.addrlist.end()) {
		if (maxAddrTotalAge < curts && it->second < curts - maxAddrTotalAge) {
			//cout << "erase: " << it->first << "\t" <<  maxAddrTotalAge  << "\t" << curts << "\t" <<  it->second  << endl;
			it = p.addrlist.erase(it);
			++statAddrDelete;
		} else {
			++it;
		}
	}

}

simtime_t readChurtrace(const char* churntracefile) {
	ifstream ifs(churntracefile);
	simtime_t ts;
	simtime_t t0 = 0;
	event e;
	size_t nPeers = 0;
	for(string line; getline(ifs, line ); ) {
	    istringstream iss(line);
	    vector<string> strs{istream_iterator<string>{iss}, istream_iterator<string>{}};

	    ts = stoi(strs[0]);
	    size_t peerid = stoi(strs[1]);
	    bool join = (strs[2] == "c");
		//cout << ts << "\t" << peerid << "\t" << join << endl;



	    // simulation starts when first peer leaves. all peers before are initialized
	    if (t0 == 0) {
	    	// init phase
	    	if (join == false) {
	    		cout << "nPeers: " << nPeers << " t0:" << ts << endl;
	    		t0 = ts;
	    		peers.resize(nPeers);
	    	}
	    	nPeers = peerid + 1;
	    }

	    if (t0 > 0) {	// not in else because t0 can be set above in if
		    if (churnmap.size() <= peerid) {
		    	churnmap.resize(peerid+1);
		    }
		    churnmap[peerid].first.emplace_back(ts, join ? ADD : DEL);

	    	e.receiver = peerid;
	    	e.type = join ? ADD : DEL;
	    	eventMap.insert(make_pair(ts, e));
	    }

	}

	/*
	 * genAdj calls establishOutgoing which schedules ADDR events which requires the churnmap to be populated completely
	 * therefore, call it down here and not before
	 */
	genAdj(t0);

	/*for (size_t i = 0; i < churnmap.size(); ++i) {
		cout << i << endl;
		for (const auto& p : churnmap[i]) {
			cout << p.first << " " << p.second << endl;
		}
	}*/

	cout << "lastTS: " << ts << " duration: " << ts -t0 << endl;

	return t0;
}

simtime_t initSim(const char* churntracefile) {
	assert(maxFakeAgeConnection <= maxAddrTotalAge);
	SetSeed(seed);
	simtime_t startTS;
	if (churntracefile != NULL) {
		startTS = readChurtrace(churntracefile);
	} else {
		startTS = startSimTSDefault;
		InitNormDist(maxFakeAgeConnection / 2, 20*oneMinute, seed);
		peers.resize(numPeersDefault);
		createChurn(startTS);
		genAdj(startTS);
	}

	event e;
	e.type = ADDRSTAT;
	eventMap.insert(make_pair(startTS, e));

	statAgeConnected.resize(maxAddrTotalAge/statCompress+1);
	statAgeNotConnected.resize(maxAddrTotalAge/statCompress+1);

	return startTS;
}


void doAddrstat(const simtime_t ts, const simtime_t startTS) {
	// nOnline: the total number of currently only peers
	int nOnline = 0;
	int nConnections = 0;
	vector<int> onlinecount;
	vector<int> addrsize;
	vector<int> ageOnline(maxAddrTotalAge/oneSecond+1);
	vector<int> ageOffline(maxAddrTotalAge/oneSecond+1);

	for (const auto& p : peers) {
		if (p.online) {
			++nOnline;
			nConnections += p.outgoing.size();
			int peersonline = 0;
			int peersoffline = 0;
			for (const auto& addr : p.addrlist) {
				size_t pos = (ts-addr.second)/oneSecond;
				if (pos >= ageOnline.size()) {
					continue;
				}
				if (peers[addr.first].online) {
					++peersonline;
					++ageOnline[(ts-addr.second)/oneSecond];
				} else {
					++peersoffline;
					++ageOffline[(ts-addr.second)/oneSecond];
				}
			}
			onlinecount.push_back(peersonline);
			addrsize.push_back(p.addrlist.size());
		}
	}

	// the average number of peers in the addrlist that are online (neighbor+foreign peers)
	double avgonlinecount = accumulate(onlinecount.begin(), onlinecount.end(), 0) / (double)onlinecount.size();

	// avgaddrsize: the average size of the complete list including online+offline, neighbors+foreign peers
	double avgaddrsize = accumulate(addrsize.begin(), addrsize.end(), 0) / (double)addrsize.size();

	// avgonlineshare: the probability of a random ip in the addrlist to be online
	double avgonlineshare = avgonlinecount / avgaddrsize;

	// Pconn: the probability of a random edge to exist
	double Pconn = nConnections/((nOnline*(nOnline-1))/2.0);

	// the average number of connections per peer
	double avgConnectionCount = nConnections*2.0 / (double)nOnline;

	// the probability that a random _online_ and _foreign_ IP is in a random addrlist
	double PHasForeignIP = (avgonlinecount - avgConnectionCount) / (nOnline - 1 - avgConnectionCount);	// -1 = own peer

	// the probability that a random peer from addrlist is being sent in one timestep
	double Psent = maxAddrSend < avgaddrsize ? maxAddrSend / avgaddrsize : 1;


	double pSendIPconnected = Psent;	// P(sendIP | connected)
	double pSendIPnotconnected = PHasForeignIP * Psent;	// P(sendIP | !connected)


	json j = {
			{"maxFakeAgeConnection", maxFakeAgeConnection},
			{"maxAddrTotalAge", maxAddrTotalAge},
			{"meanAddrSendInterval", meanAddrSendInterval},
			{"maxAddrSend", maxAddrSend},
			{"seed", seed},
			{"ts", ts - startTS},

			{"avgonlineshare", avgonlineshare},
			{"avgonlinecount", avgonlinecount},
			{"avgaddrsize", avgaddrsize},
			{"nConnections", nConnections},
			{"nOnline", nOnline},
			{"correlation", calcRankCorrelation(ageOnline, ageOffline)},
			{"avgConnectionCount", avgConnectionCount},
			{"PHasForeignIP", PHasForeignIP},

			{"Pconn", Pconn},
			{"pSendIPconnected", pSendIPconnected},
			{"pSendIPnotconnected", pSendIPnotconnected},

			{"statAgeConnected", statAgeConnected},
			{"statAgeNotConnected", statAgeNotConnected},

			{"statEventJoin",  statEventJoin},
			{"statEventLeave",  statEventLeave},
			{"statEventAddr",  statEventAddr},
			{"statEventAddrScheduled",  statEventAddrScheduled},
			{"statEventAddrCancelled",  statEventAddrCancelled},

			{"statAddrConnectionNew",  statAddrConnectionNew},
			{"statAddrConnectionExisting",  statAddrConnectionExisting},
			{"statAddrMsgNew",  statAddrMsgNew},
			{"statAddrMsgUpdateNewer",  statAddrMsgUpdateNewer},
			{"statAddrMsgUpdateOlder",  statAddrMsgUpdateOlder},
			{"statAddrUpdate",  statAddrUpdate},
			{"statAddrDelete",  statAddrDelete},
			{"statAddrOffline",  statAddrOffline}
	};

	jsonvec.push_back(j);
	json jvec(jsonvec);

	ofstream ofs(jsonfname);
	ofs.width(2);
	ofs << jvec;
}

void printEventMap() {
	while (!eventMap.empty()) {
		simtime_t ts = eventMap.begin()->first;
		event &e = eventMap.begin()->second;
		cout << ts << "\t" << e.type << "\t" << e.receiver << "\t" << e.data << endl;
		eventMap.erase(eventMap.begin());
	}
}

void sim(const simtime_t startTS) {
	cout << "start sim" << endl;
	simtime_t ts = startTS;
	size_t eventcounter = 0;
	while (!eventMap.empty() && ts < startTS + maxSimTime) {
		ts = eventMap.begin()->first;
		event &e = eventMap.begin()->second;
		//cout << ts << "\t" << e.type << "\t" << e.receiver << "\t" << e.data << endl;
		if (++eventcounter % 100000 == 0) {
			cout << ts - startTS << "\t" << eventcounter << " mem: " << memory_used()/1024/1024 << " mb" << endl;
		}

		switch (e.type) {
		case ADD:
			LOG(ts << " add peer " << e.receiver << endl);
			++statEventJoin;
			if (e.receiver == peers.size()) {
				// new peer added
				peers.resize(e.receiver+1);
			} else if (e.receiver < peers.size()) {
				// peer rejoined, incoming, outgoing, addrlist were cleared on exit,
				// just set online
				assert(peers[e.receiver].online == false);
				peers[e.receiver].online = true;
			} else {
				cerr << "e.receiver > peers.size()" << endl;
				assert(false);
			}

			assert(e.receiver < churnmap.size());
			assert(churnmap[e.receiver].first[churnmap[e.receiver].second].first == ts);
			assert(churnmap[e.receiver].first[churnmap[e.receiver].second].second == e.type);
			++churnmap[e.receiver].second;

			establishOutgoing(ts, e.receiver);

			break;

		case DEL:
		{
			LOG(ts << " del peer " << e.receiver << endl);
			++statEventLeave;
			assert(e.receiver < peers.size());
			assert(peers[e.receiver].online == true);

			assert(e.receiver < churnmap.size());
			assert(churnmap[e.receiver].first[churnmap[e.receiver].second].first == ts);
			assert(churnmap[e.receiver].first[churnmap[e.receiver].second].second == e.type);
			++churnmap[e.receiver].second;

			peer &p = peers[e.receiver];
			p.online = false;

			for (const peerid_t r : p.outgoing) {
				peers[r].incoming.erase(e.receiver);
			}
			for (const peerid_t r : p.incoming) {
				peers[r].outgoing.erase(e.receiver);
				establishOutgoing(ts, r);	// let remote create new connection for lost outgoing connection
			}

			p.incoming.clear();
			p.outgoing.clear();
			statAddrOffline += p.addrlist.size();
			p.addrlist.clear();

			break;
		}
		case SENDADDR:
		{
			assert(e.receiver < peers.size() && e.data < peers.size());
			peer &p = peers[e.receiver];
			if (p.online == false) {
				++statEventAddrCancelled;
				//LOG("peer went offline" << endl);
				break;
			}
			if (peers[e.data].online == false) {
				++statEventAddrCancelled;
				//LOG("neighbor went offline" << endl);
				break;
			}
			/*
			 * this alone does not prevent peers from going offline and rejoining
			 * and the events from the old session are still executed
			 * --> prevent events that will be executed after a peer left the network
			 * from being scheduled (see scheduleSendaddrEvent())
			 */
			++statEventAddr;
			LOG(ts << " sendaddr " << e.receiver << " -> " << e.data << endl);
			//cout << ts << " sendaddr " << e.receiver << " -> " << e.data << endl;

			// check whether connected IPs need updated timestamps
			updateConnectionTimestamps(p, ts);

			// select subset of addrlist
			set<peerid_t> tosend;
			if (maxAddrSend >= p.addrlist.size()) {
				// send whole addrlist
				for(const auto& addr : p.addrlist) {
					handleReceiveAddr(ts, e.receiver, e.data, addr);
				}
			} else {
				while (tosend.size() < maxAddrSend) {
					tosend.insert(GetRandInt(p.addrlist.size()));	// random_shuffle might be faster
				}

				auto it = p.addrlist.begin();
				peerid_t previousid = 0;
				for (peerid_t entryid : tosend) {
					//cout << "advance to " << entryid << " total size: " << p.addrlist.size() << endl;
					advance(it, entryid - previousid);
					handleReceiveAddr(ts, e.receiver, e.data, *it);
					previousid = entryid;
				}
			}


			// schedule next sendaddr event
			scheduleSendaddrEvent(ts, e.receiver, e.data, false);

			break;
		}
		case ADDRSTAT:
			doAddrstat(ts, startTS);
			eventMap.insert(make_pair(ts + statInterval, e));

			break;

		}	// switch end

		eventMap.erase(eventMap.begin());
	}
	doAddrstat(ts, startTS);


}

void reset() {
	eventMap.clear();
	peers.clear();

	statAgeConnected.clear();
	statAgeNotConnected.clear();
}

void usage() {
	cout << "Simulates a simple IP address propagation mechanism" << endl;
	cout << "USAGE: ./sim [churntrace tracefile.txt][d maxFakeAgeConnection] [dx maxAddrTotalAge] "
			<< "[delta meanAddrSendInterval] [n maxAddrSend] [maxsimtime maxSimTime] [seed seed]" << endl;
}
int main(int argc, char** argv) {

	stringstream ss;
	char* churntracefile = NULL;

	for (int i = 1; i < argc - 1; i += 2) {	// increment by TWO
		if (argv[i] == string("churntrace")) {
			churntracefile = argv[i+1];
			cout << "using churntracefile " << churntracefile << endl;
		} else if (argv[i] == string("d")) {
			maxFakeAgeConnection = atoi(argv[i+1]);
			cout << "d = " << maxFakeAgeConnection << endl;
		} else if (argv[i] == string("dx")) {
			maxAddrTotalAge = atoi(argv[i+1]);
			cout << "dx = " << maxAddrTotalAge << endl;
		} else if (argv[i] == string("delta")) {
			meanAddrSendInterval = atoi(argv[i+1]);
			cout << "delta = " << meanAddrSendInterval << endl;
		} else if (argv[i] == string("n")) {
			maxAddrSend = atoi(argv[i+1]);
			cout << "n = " << maxAddrSend << endl;
		} else if (argv[i] == string("maxsimtime")) {
			maxSimTime = atoi(argv[i+1]);
			cout << "maxSimTime = " << maxSimTime << endl;
		} else if (argv[i] == string("seed")) {
			seed = atoi(argv[i+1]);
			cout << "seed = " << seed << endl;
		} else {
			cout << "unrecognized option: " << argv[i] << endl;
			usage();
			return -1;
		}
	}
	if (argc % 2 == 0) {
		cout << "strange argument count" << endl;
		usage();
		return -1;
	}

	stringstream fname;
	fname << "simresult_d" << maxFakeAgeConnection << "_dx" << maxAddrTotalAge
			<< "_delta" << meanAddrSendInterval << "_n" << maxAddrSend << ".json";
	jsonfname = fname.str();




	simtime_t startTS = initSim(churntracefile);
	sim(startTS);
}
