### SDL Simple example
1. To setup, run `source setup.sh`
2. `make` twice - once inside the SDL directory, and once outside
3. Download http://uaf-10.t2.ucsd.edu/~bsathian/SDL/test_10_events.root
4. Run the following command
`./bin/doAnalysis -i test_10_events.root -t trackingNtuple/tree -n 1 -d -v 1`
5. Tweak the argument of -n to change the number of events (up to 10)
