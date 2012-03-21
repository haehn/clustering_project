#include <fstream>
#include <iostream>

#include "MinSqTree.h"
#include "ProblemParser.h"

using namespace std;

int main(int argc, char *argv[]) {

   std::vector<MinSquareTreeCollection::DblMatrix> pmatrices;
   MinSquareTreeCollection::IntMatrix pmapping;
   std::vector<std::string> plabels;
   PhyTree *ptree = NULL;

   MinSquareTreeCollection *mstc = NULL;

   ifstream matrices_data(argv[1]);          // command line args
   ifstream mapping_data(argv[2]);           // used to direct input
   ifstream labels_data(argv[3]);            // files to program
   ifstream tree_data(argv[4]);              // order is: distvar, genomemap, labels, tree
   ofstream outfile(argv[5]);                // fifth argument is output newick filename
   string tree;

   try {
      pmatrices = ProblemParser::parse_matrices(matrices_data);
      pmapping = ProblemParser::parse_mapping(mapping_data);
      plabels = ProblemParser::parse_labels(labels_data);
      ptree = ProblemParser::parse_tree(tree_data);

      mstc = new MinSquareTreeCollection(pmatrices,pmapping,plabels,*ptree);
      delete ptree;
      ptree = NULL;

      mstc->compute(false);
   }
   catch(std::exception &e) {
      if(ptree) {
         delete ptree;
      }
      if(mstc) {
         delete mstc;
      }
      cerr << e.what() << endl;
   }

   tree = mstc->getTree();                   // output tree as a string
   cout << tree << endl;                     // write to stdout
   cout << mstc->getScore() << endl;
   outfile << tree << endl;                  // write to file
   outfile.close();
   
   delete mstc;
   return 0;
}
