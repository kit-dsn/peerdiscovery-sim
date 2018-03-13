# Overview

Simulates a simple IP address propagation mechanism and an attacker that wants to infer the network topology by observing address messages.


# Requirements & Build
Requires C++11, no additional dependencies.

Run ``make`` to build all required executables.


# Usage

Run ./sim to simulate the address propagation mechanism of a network. Parameters to ./sim are:
```
./sim [churntrace tracefile.txt][d maxFakeAgeConnection] [dx maxAddrTotalAge] [delta meanAddrSendInterval] [n maxAddrSend] [maxsimtime maxSimTime] [seed seed]
```

  * [churntrace churntrace.txt]: use given tracefile for churn of peers. If option is not given, churn is created based on some hard-coded normal distribution.
  * [d maxFakeAgeConnection]: The maximum age of connected IP addresses ($\delta_d$).
  * [dx maxAddrTotalAge]: Not connected addresses older dx are removed from a client's address list ($\delta_x$).
  * [delta meanAddrSendInterval]: How often clients send address messages to their peers ($\delta_s$).
  * [n maxAddrSend]: Number of addresses to send per address message ($n$).
  * [maxsimtime maxSimTime]: Maximum simulation time in seconds.
  * [seed seed]: Seed for PRNG.

./sim produces a json file that contains the results of the simulation run. This file can be used to simulate a topology inference attack using ./simattack.
```
./simattack <simresult.json> <ts> <nObservations> <nTrials> <seed>
```

  * simresult.json: The file created by ./sim
  * ts: The timestamp from the simulation that should be used. ```grep ts *.json``` to get available timestamps from json file.
  * nObservations: Maximum number of observations (i.e., observed address messages) by the simulated adversary.
  * nTrials: Number of random trials for each observation. (Use at least 100,000 to get any results)
  * seed: Seed for PRNG.
  
./simattack prints a table with the columns nObservations, maxFakeAgeConnection, maxAddrSend, precision, recall (in that order).