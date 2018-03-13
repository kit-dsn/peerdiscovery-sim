#pragma once

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;


template<typename T>
void writePlotData(const vector<T>& v, ofstream& ofs, bool withe = true) {
	for (auto it = v.begin(); it != v.end(); ++it) {
		if (std::isnan(*it) == false) {
			ofs << *it << endl;
		} else {
			ofs << 0 << "\n";
		}
	}
	if (withe) {
		ofs << "e" << "\n";
	}
}

template<typename T>
void gpplot(const vector<T>& v, const vector<T>& w, const char* xrange = NULL, const char* yrange = NULL,
		const char* desc = NULL, bool store = false, bool png = false) {

	ofstream gpfile;
	ofstream gpdfile;
	char fname[100] = "__cpp.gp";
	char fnamegpd[100] = "__cpp.gpd";
	char cmd[100];

	if (store) {
		if (desc != NULL) {
			snprintf(fnamegpd, sizeof(fnamegpd), "%s.gpd", desc);
			snprintf(fname, sizeof(fname), "%s.gp", desc);
		}
		gpdfile.open(fnamegpd);
	}

	gpfile.open(fname);



	if (png == true) {
		gpfile << "set terminal png" << endl;
		//assert (desc != NULL);
		gpfile << "set output '" << desc << ".png'" << endl;
	}
	gpfile << "set grid lw 0.5" << endl;

	if (xrange != NULL) {
		gpfile << "set xrange " << xrange << endl;
	}
	if (yrange != NULL) {
		gpfile << "set yrange " << yrange << endl;
	}

	if (desc != NULL) {
		gpfile << "set title '" << desc << "'" << endl;
	}
	gpfile << "plot '-' title '1', '-' title '2'" << endl;

	writePlotData<T>(v, gpfile);
	if (store) {
		writePlotData<T>(v, gpdfile);
	}

	writePlotData<T>(w, gpfile);
	if (store) {
		writePlotData<T>(w, gpdfile);
	}

	if (png == false) {
		gpfile << "pause -1" << endl;
	}
	gpfile.close();
	
	if (store)
		gpdfile.close();

	if (!store) {
		snprintf(cmd, sizeof(cmd), "gnuplot %s", fname);
		int rc = system(cmd);
		if (rc == -1) {
			cerr << "gnuplot call to system failed" << endl;
		}
	}

}

template<typename T>
void gpplot(vector<T> v, const char* xrange = NULL, const char* yrange = NULL,
		const char* desc = NULL, bool store = false, bool png = false, bool execute = true) {

	ofstream gpfile;
	ofstream gpdfile;
	char fname[100] = "__cpp.gp";
	char fnamegpd[100] = "__cpp.gpd";
	char cmd[100];

	if (desc != NULL) {
		snprintf(fnamegpd, sizeof(fnamegpd), "%s.gpd", desc);
		snprintf(fname, sizeof(fname), "%s.gp", desc);
	}

	if (store) {
		gpdfile.open(fnamegpd);
	}
	gpfile.open(fname);

	if (png == true) {
		gpfile << "set terminal png" << endl;
		//assert (desc != NULL);
		gpfile << "set output '" << desc << ".png'" << endl;
	}
	gpfile << "set grid lw 0.5" << endl;

	if (xrange != NULL) {
		gpfile << "set xrange " << xrange << endl;
	}
	if (yrange != NULL) {
		gpfile << "set yrange " << yrange << endl;
	}

	if (desc != NULL) {
		gpfile << "set title '" << desc << "'" << endl;
	}
	gpfile << "plot '-' " << endl;

	writePlotData<T>(v, gpfile);
	if (store) {
		writePlotData<T>(v, gpdfile, false);
	}

	if (png == false) {
		gpfile << "pause -1" << endl;
	}
	gpfile.close();

	if (store)
		gpdfile.close();

	if (execute) {
		snprintf(cmd, sizeof(cmd), "gnuplot %s", fname);
		int rc = system(cmd);
		if (rc == -1) {
			cerr << "gnuplot call to system failed" << endl;
		}
	}


}
