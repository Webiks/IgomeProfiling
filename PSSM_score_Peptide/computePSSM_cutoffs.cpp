#include "computePSSM_cutoffs.h"

#include <fstream>
using namespace std;
#include "SEQ.h"

void PSSM_scoresFromSeqVector(const PSSM& PSSM1, const vector<vector <size_t>> & vseq1, vector<double> & scores);

computePSSM_cutoffs::computePSSM_cutoffs(vector<PSSM> & PSSM_array,
	size_t TotalNumberOfRandoSeq,
	alphabet& alph,
	const string & CutofsPerPSSM_FileName,
	int totalMemes,
	double PercentOfRandomHitsPerPSSM) : 
	_PSSM_array(PSSM_array), 
	_totalNumberOfRandoSeq(TotalNumberOfRandoSeq), 
	_alph(alph), 
	_CutofsPerPSSM_FileName(CutofsPerPSSM_FileName),
	_totalMemes(totalMemes),
	_PercentOfRandomHitsPerPSSM(PercentOfRandomHitsPerPSSM)
{
	generateRandomPeptides();
	computecCutoffsBasedOnRandomPeptides();
}

string SizetToString(size_t sz) {
	stringstream ss;
	ss << sz;
	return ss.str();
}

void computePSSM_cutoffs::generateRandomPeptides() {
	size_t NumberOfRandoSeq = _totalNumberOfRandoSeq;  
	srand(931); // Set srand for generating random pepties // TODO set srand from input argument
	map<string, randomPeptides>::iterator it = _randomPeptideDataSet.begin(); // use iteration and insert to add values to map
	for (size_t length = 5; length <= 12; length++)
	{
		randomPeptides tmp(_alph._aaFreq, NumberOfRandoSeq, length);
		tmp.generateRandomSequences();
		string seqType = SizetToString(length);
		_randomPeptideDataSet.insert(it, pair<string, randomPeptides>(seqType, tmp));
		
		cout << "Finish to create " << NumberOfRandoSeq << " random sequences of type " << seqType << endl;
		//randomPeptideDataSet.insert([seqType]=tmp;
		randomPeptides tmpCys(_alph._aaFreq, NumberOfRandoSeq, length, 1);
		tmpCys.generateRandomSequences();
		seqType = "C" + SizetToString(length) + "C";
		// randomPeptideDataSet[seqType]=tmpCys;
		_randomPeptideDataSet.insert(it, pair<string, randomPeptides>(seqType, tmpCys));
		cout << "Finish to create " << NumberOfRandoSeq << " random sequences of type " << seqType << endl;
	}
}

double get_percentile_cutoff(const vector<double> & scores, double TopPercent)
{
	vector <double> sorted_vector = scores;
	sort(sorted_vector.begin(), sorted_vector.end());
	double TopPercentSize = TopPercent*scores.size();
	if (floor(TopPercentSize) == 0) {
		cout << "ERROR: too small random datast for determining PSSM score cutoff at level of " << TopPercent << " (TopPercentSize is: " << TopPercentSize << ")" << endl;
		exit(3);
	}
	int TopPercentIndex = int(scores.size() - floor(TopPercentSize));
	double cutoff = sorted_vector[TopPercentIndex];
	return (cutoff);
}

void computePSSM_cutoffs::computecCutoffsBasedOnRandomPeptides() {
	
	for (size_t i = 0; i<_PSSM_array.size(); ++i)
	{
		for (auto SeqDataset = _randomPeptideDataSet.begin(); SeqDataset != _randomPeptideDataSet.end(); ++SeqDataset)
		{
			vector<double> tmp;
			PSSM_scoresFromSeqVector(_PSSM_array[i], SeqDataset->second._simulatedSequences, tmp);
			_PSSM_array[i].dist_of_random_peptides_scores[SeqDataset->first] = tmp;
		}
		cout << "Score PSSM " << i << " for random dataset" << endl;
	}
	
	// determin cutoff for each seq type and print to file
	ofstream PSSM_Scores_Cutoff;
	PSSM_Scores_Cutoff.open(_CutofsPerPSSM_FileName);
	double PercentOfAcceptedPeptidesPerType = _PercentOfRandomHitsPerPSSM / _totalMemes; // for each seq type the cutoff will be the percent of hits accepted for all PSSMs devided by the number of PSSMs
	for (size_t i = 0; i<_PSSM_array.size(); ++i)
	{
		PSSM_Scores_Cutoff << "###\t" << _PSSM_array[i].PSSM_name << "\t";
		SEQ tmp(_PSSM_array[i]._PSSM_consensus_seq,"",1,_alph);
		PSSM_Scores_Cutoff << tmp.getStringOfSeq() << "\t";
		PSSM_Scores_Cutoff << _PSSM_array[i].PSSM_MaxScore << endl;
		// go over all seq types
		for (auto SeqTypeRandomPeptidesScores = _PSSM_array[i].dist_of_random_peptides_scores.begin(); SeqTypeRandomPeptidesScores != _PSSM_array[i].dist_of_random_peptides_scores.end(); ++SeqTypeRandomPeptidesScores)
		{
			vector <double> scores = SeqTypeRandomPeptidesScores->second;
			string seqType = SeqTypeRandomPeptidesScores->first;
			double cutoff = get_percentile_cutoff(scores, PercentOfAcceptedPeptidesPerType);
			_PSSM_array[i].set_SeqTypCutoff(seqType, cutoff);
			PSSM_Scores_Cutoff << seqType << "\t" << cutoff << endl;
		}
		cout << "Define cutoff for PSSM #" << i << endl;
	}
	PSSM_Scores_Cutoff.close();
}