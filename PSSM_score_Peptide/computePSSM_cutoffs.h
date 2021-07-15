#pragma once
#ifndef COMPUTE_PSSM_CUTOFFS_____H
#define COMPUTE_PSSM_CUTOFFS_____H

#include <vector>
#include <map>
#include <string>
#include <cmath>  // because of unix
using namespace std;

#include "PSSM.h"
#include "randomPeptides.h"
#include "alphabet.h"

class computePSSM_cutoffs{
public:
	computePSSM_cutoffs(vector<PSSM> & PSSM_array,
						size_t TotalNumberOfRandoSeq,
						alphabet & alph,
						const string & CutofsPerPSSM_FileName,
						int totalMemes,
						double PercentOfRandomHitsPerPSSM);

private:
	void generateRandomPeptides();
	void computecCutoffsBasedOnRandomPeptides();
	

	vector<PSSM>  & _PSSM_array;
	size_t _totalNumberOfRandoSeq;
	//vector<double> const & _aaFreq;
	//vector<char> const & _correspondingCharacters;
	map<string, randomPeptides> _randomPeptideDataSet; // map between seqType and object containing the random peptides dataset
	string const & _CutofsPerPSSM_FileName;
	alphabet& _alph;
	int _totalMemes;
	double _PercentOfRandomHitsPerPSSM;
















};

#endif
