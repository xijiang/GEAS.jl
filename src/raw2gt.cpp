#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <vector>

/**
 * Remove the header line.
 * Remove the first 6 columns of a plink.raw file.
 * Remove all the spaces between the genotypes.
 * Write genotypes to stdout
 * Calculate the allele frequencies into argv[1].
 */
using namespace std;

int main(int argc, char *argv[])
{
  ios_base::sync_with_stdio(false);
  if(argc !=2){
    cerr << "Usage: cat plink.raw | "<<argv[0]<<" freq.txt >genotypes"<<endl;
    return 1;
  }
  string line, t;
  int nlc{-6}, div{0};
  getline(cin, line);
  {
    stringstream ss(line);
    string t;
    while(ss>>t) ++nlc;
  }
  double frq[nlc]{0.};
  while(getline(cin, line)){
    stringstream ss(line);
    for(auto i{0}; i<6; ++i) ss>>t;
    getline(ss, line);
    line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
    cout<<line<<'\n';
    for(auto i=0; i<nlc; ++i) frq[i] += line[i] - '0';
    div += 2;
  }
  
  ofstream foo(argv[1]);
  foo.precision(17);
  for(auto&x:frq) foo<<x/div<<'\n';
  return 0;
}
