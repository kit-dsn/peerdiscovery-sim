all:
	g++ -Wall -Wno-unused-function -g --std=c++11 -O3 sim.cc -o sim
	g++ -Wall -Wno-unused-function -g --std=c++11 -O3 simattack.cc -o simattack
