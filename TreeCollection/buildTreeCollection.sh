g++ -DDEBUG -Wall -O0 -g3 -I/usr/local/include/eigen3 ProblemParser.cc TreeCollection.cc MinSqTree.cc newick.cc -o TreeCollection.debug
g++ -Wall -O3 -march=corei7 -g3 -I/usr/local/include/eigen3 ProblemParser.cc TreeCollection.cc MinSqTree.cc newick.cc -o TreeCollection
