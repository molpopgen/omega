/* omega.cc
	Copyright 2008 Kevin Thornton, UC Irvine
		krthornt@uci.edu
		
		This code is released under the terms of the GNU public license.
		See http://www.gnu.org/copyleft/gpl.html for license details
		
		To compile:
			c++ -o omega omega.cc -O2 -lsequence

			c++ -o omega omega.cc -O2 -lsequence -lz (with libsequence >= 1.7.8)
			
		To run:
			ms [params] | ./omega
			
		Note: singletons (based on minor allele frequency) are excluded from the analysis
*/

#include <Sequence/SimData.hpp>
#include <Sequence/PolySNP.hpp>
#include <Sequence/SimParams.hpp>
#include <limits>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
using namespace std;
using namespace Sequence;

std::pair<double,double> omega_max(const SimData * data)
{
  //vector<string> temp(1, std::string(data->numsites(),'0'));
  //copy(data->begin(),data->end(),back_inserter(temp));
  
  //SimData d2;
  //  d2.assign(&*(data->pbegin()),data->numsites(),&temp[0],temp.size());
  SimData d2(*data);
  d2.ApplyFreqFilter(2);
  if(d2.empty()) { return make_pair(strtod("NAN",NULL),-1); }
  //cerr << d2 << endl;
  //  exit(1);
  //PolySNP adata(&d2,true);
  PolySNP adata(&d2);
  vector<vector<double> > ld = adata.Disequilibrium(2);

  double omega_max = numeric_limits<double>::min();
  double snp=-1;
  unsigned S = d2.numsites();
  for(unsigned pos = 1 ; pos < S-1 ; ++pos)
    {
      unsigned l=pos+1;
      double position = d2.position(pos);
      double srsqL=0.,srsqR=0.,sumrsqLR=0.;
      for(unsigned i=0;i<ld.size();++i)
	{
	  //cerr << position << ' ' << l << ' ' << S <<' ' << ld[i][0] << ' ' << ld[i][1] 
	  //<<' ' << ld[i][3] << '\n';
	  //assert(ld[i][0] != ld[i][1]);
	  if( ld[i][1] <= position )
	    {
	      srsqL += ld[i][2];
	    }
	  else if (ld[i][1] > position)
	    {
	      srsqR += ld[i][2];
	    }
	  if( ld[i][0]<=position && ld[i][1] > position)
	    {
	      sumrsqLR += ld[i][2];
	    }
	}
      //cerr << srsqL << ' ' << srsqR << ' ' << sumrsqLR << '\n';
      double numerator = (1./(double(l*(l-1))/2.+ double( (S-l)*(S-l-1) )/2.))*(srsqL+srsqR);
      double denominator = (1./(double(l*(S-l))))*sumrsqLR;
      double omega = numerator/denominator;
      if(finite(omega) && omega > omega_max)
	{
	  omega_max=omega;
	  snp = position;
	}
    }
  //cerr << *data << '\n';
  return make_pair(omega_max,snp);
}

int main(int argc, char **argv)
{
  SimData d;
  int rv;
	string temp;
getline(cin,temp);
cerr << temp << '\n';
  while( (rv=d.fromfile(stdin)) != EOF)
    {
      pair<double,double> o = omega_max(&d);
      cout << o.first << '\t' << o.second << endl;
    }
}
