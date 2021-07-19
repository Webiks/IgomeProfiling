#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>
#include <cfloat> // eli because of unix
#include <cmath>  // because of unix
#include <random>
using namespace std;

#include "PSSM.h"
#include "SEQ.h"
#include "HIT.h"
#include "randomPeptides.h"
#include "computePSSM_cutoffs.h"
#include "read_PSSMs_from_MAST_file.h"
#include "alphabet.h"


string PadSeq (string seq, size_t PaddingLength);
void fromFileToVectorOfString(const string fileName,vector<string> & allLines);
void readFileToPSSM_array(const string fileName, vector <PSSM> & PSSM_array, map<string,int> & PSSM_Name_To_ArrayIndex);
void readFileToSeq_array (const string fileName, alphabet& alph, vector <SEQ> &Seq_array, bool useCountDuplicate);
void get_top_hits (const vector<SEQ> & sorted_seq,double fraction, vector <SEQ> & top_hits);
void fill_Seq_Hits_PSSM_map(const vector<SEQ> & seq_Hits_vector,map<string,vector<string>> &PSSM_Hits, string PSSM_name);
vector<SEQ> get_sort_seq_vector_by_scores (vector<SEQ> & seq_vector, vector <double> & scores_vector);
// void Find_PSSM_Hits(const PSSM& PSSM1, const vector<SEQ> & vseq1, vector<SEQ> & hits);
void Find_PSSM_Hits(const PSSM& PSSM1, const vector<SEQ> & vseq1, vector<HIT> & hits, const size_t verbose);
void get_aaFreq_and_CorrChar (string const PSSM_FileName, vector <double> & aaFreq, vector <char> & correspondingCharacters);
void CalclateCutoffsForSetOfPSSMs (vector <PSSM> & PSSM_array,string RandomPeptides_FileName,size_t TotalNumberOfRandoSeq, vector <double> const & aaFreq, vector<char> & correspondingCharacters, const string & RandomHits_FileName, const string & CutofsPerPSSM_FileName);
void get_PSSM_cutoff_scores_from_file (const string fileName,map<string,int> & PSSM_Name_To_ArrayIndex,vector <PSSM> & PSSM_array);
// void print_PSSM_Hits (const PSSM & PSSM,const vector <SEQ> & hits, ostream& out);
void print_PSSM_Hits (const PSSM & PSSM,const vector <HIT> & hits, ostream& out);
string SizetToString(const size_t sz);

// The program can be used for two reasons
// 1. To compute cutoffs for an array of PSSM (from which score, the match to a peptide is statistically significant). This is main1().
// 2. Given a set of peptides in a Fasta file and a set of PSSMs, asign for each PSSM the peptides that are significantly associated with it. This is main2().


int computeCutoffsOfPssmMain(int argc, char *argv[]);
int associatePeptidesWithPSSM(int argc, char *argv[]);
int assignPvalueToPSSMaRRAY(int argc, char *argv[]);
size_t get_running_mode(int argc, char *argv[]);
int main(int argc, char *argv[]) {
	size_t running_mode=get_running_mode(argc, argv);
	cout<<"MODE: "<<running_mode<<endl;
	if (running_mode == 3){
		cout<<"MODE: 3, call: associatePeptidesWithPSSM()"<<endl;
		return associatePeptidesWithPSSM(argc, argv);
	}
	else if (running_mode == 2){
		cout<<"MODE: 2, call: assignPvalueToPSSMaRRAY();"<<endl;
		return assignPvalueToPSSMaRRAY(argc, argv);
	}
	else if (running_mode == 1){
		cout<<"MODE 1, call computeCutoffsOfPssmMain()"<<endl;
		return computeCutoffsOfPssmMain(argc, argv);
	}
}

// The function computeCutoffsOfPssmMain gets two parameters in the command line. That is argc = 4.
// The first is -pssm <name_of_pssm_in_Mast_format>. 
// The second is- pssm_cutoffs <filename_for PSSM_cutoffs>
size_t get_running_mode(int argc, char *argv[]){
	int max_required_params = 8;
	int min_required_params = 2;
	if ((argc > (max_required_params * 2) + 2) || (argc < (min_required_params * 2) + 2)) {// each with its flag and mode flag, check the value of argc. If not enough parameters have been passed, inform user and exit.
		cout << "Usage is in one of few modes: <<endl";
		cout << "[1: CalcPSSM_Cutoff] " << argv[0] << " -pssm <PSSMs_in_MAST_Format> -pssm_cutoffs <filename_for PSSM_cutoffs> -CalcPSSM_Cutoff -total_memes <total memes, used if input is splitted otherwise 0>" << endl; // Inform the user of how to use the program
		cout << "[2: CalcPSSM_Pval] "<<argv[0]<<" -pssm <PSSMs_in_MAST_Format> -pssm_cutoffs <filename_for PSSM_cutoffs> -seq <input_seq_FASTA> -out <out> -NrandPSSM <number_of_random_PSSMs> -userFactor <is use factor> -CalcPSSM_Pval" << endl; // Inform the user of how to use the program
		cout << "[3: CalcPSSM_Hits] " << argv[0] << " -pssm <PSSMs_in_MAST_Format> -pssm_cutoffs <filename_for PSSM_cutoffs> -seq <input_seq_FASTA> -out <out> -CalcPSSM_Hits " << endl; // Inform the user of how to use the program
		cout << "\n\nThe number of provided arguments is "<<argc<<endl;
		exit(12);
	}
	size_t mode=0;
	for (int i = 1; i < argc; ++i) // iterate over flags
	{
		if (string(argv[i]) == "-CalcPSSM_Hits") 
			mode=3;
		else if (string(argv[i]) == "-CalcPSSM_Pval")
			mode=2;
		else if (string(argv[i]) == "-CalcPSSM_Cutoff")
			mode=1;
	}
	if ((mode!=1.0) && (mode!=2.0) && (mode!=3.0)) 
	{
		cout<<"ERROR: could not recognize running mode ["<<mode<<"]"<<endl;
		cout<<"The command was: "<<endl;
		for (int i=0;i<argc;++i)
			cout<<argv[i]<<" ";
		cout<<endl;    
		exit (24);
	}
	return mode;
}

void getFileNamesFromArgv(int argc, char *argv[], string & PSSM_FileName, string & CutofsPerPSSM_FileName, int & totalMemes, double & PercentOfRandomHitsPerPSSM) {
	// parse ARGV arguments
	size_t num_required_params = 4;
	if (argc != (num_required_params * 2)+2) {// each with its flag and mode_flag, check the value of argc. If not enough parameters have been passed, inform user and exit.
		cout << "Usage is -pssm <PSSMs_in_MAST_Format> -pssm_cutoffs <filename_for PSSM_cutoffs> -total_memes <total memes, used if input is splitted otherwise 0>\n"; // Inform the user of how to use the program
		exit(17);
	}
	cout << argv[0] <<" ";
	for (int i = 1; i < argc; ++i) // iterate over flags
	{
		/* We will iterate over argv[] to get the parameters stored inside.
		* Note that we're starting on 1 because we don't need to know the
		* name of the program, which is stored in argv[0] */
		if (string(argv[i]) == "-pssm") PSSM_FileName = string(argv[i + 1]);
		else if (string(argv[i]) == "-pssm_cutoffs") CutofsPerPSSM_FileName = string(argv[i + 1]);
		else if (string(argv[i]) == "-total_memes") totalMemes = stoi(string(argv[i + 1]));
		else if (string(argv[i]) == "-cutoff_random_peptitdes_percentile") PercentOfRandomHitsPerPSSM = stod(string(argv[i + 1]));
	}
}

void getFileNamesFromArgv(int argc, char *argv[], string & PSSM_FileName, string & CutofsPerPSSM_FileName, string & Seq_FASTA_FileName, string & Hits_Out_FileName) {
	// parse ARGV arguments
	size_t num_required_params = 4;
	if (argc != (num_required_params * 2)+2) {// each with its flag and mode_flag, check the value of argc. If not enough parameters have been passed, inform user and exit.
		cout << "Usage is -pssm <PSSMs_in_MAST_Format> -pssm_cutoffs <filename_for PSSM_cutoffs> -seq <input_seq_FASTA> -out <out>\n"; // Inform the user of how to use the program
		exit(12);
	}
	cout << argv[0] <<" ";
	for (int i = 1; i < argc; ++i) // iterate over flags
	{
		/* We will iterate over argv[] to get the parameters stored inside.
		* Note that we're starting on 1 because we don't need to know the
		* name of the program, which is stored in argv[0] */
		if (string(argv[i]) == "-pssm") PSSM_FileName = string(argv[i + 1]);
		else if (string(argv[i]) == "-pssm_cutoffs") CutofsPerPSSM_FileName = string(argv[i + 1]);
		else if (string(argv[i]) == "-seq") Seq_FASTA_FileName = string(argv[i + 1]);
		else if (string(argv[i]) == "-out") Hits_Out_FileName = string(argv[i + 1]);
	}
}

void getFileNamesFromArgv(int argc, char *argv[], string & PSSM_FileName, string & CutofsPerPSSM_FileName, string & Seq_FASTA_FileName, string & Hits_Out_FileName, size_t & numberOfRandomPSSM) {
	// parse ARGV arguments
	size_t num_required_params = 5;
	if (argc != (num_required_params * 2)+2) {// each with its flag and mode_flag, check the value of argc. If not enough parameters have been passed, inform user and exit.
		cout << "Usage is -pssm <PSSMs_in_MAST_Format> -pssm_cutoffs <filename_for PSSM_cutoffs> -seq <input_seq_FASTA> -out <out>\n"; // Inform the user of how to use the program
		exit(12);
	}
	cout << argv[0] <<" ";
	for (int i = 1; i < argc; ++i) // iterate over flags
	{
		/* We will iterate over argv[] to get the parameters stored inside.
		* Note that we're starting on 1 because we don't need to know the
		* name of the program, which is stored in argv[0] */
		if (string(argv[i]) == "-pssm") PSSM_FileName = string(argv[i + 1]);
		else if (string(argv[i]) == "-pssm_cutoffs") CutofsPerPSSM_FileName = string(argv[i + 1]);
		else if (string(argv[i]) == "-seq") Seq_FASTA_FileName = string(argv[i + 1]);
		else if (string(argv[i]) == "-out") Hits_Out_FileName = string(argv[i + 1]);
		else if (string(argv[i]) == "-NrandPSSM") numberOfRandomPSSM=size_t(atoi(argv[i + 1]));
		cout<<argv[i]<<" ";
	}
	if (numberOfRandomPSSM == 0) {
		cout<<"ERROR: could not parse the number of random PSSMs for Pvalue computation"<<endl;
		exit (99);
	}
	cout<<endl;
}

void getFileNamesFromArgv(int argc, char *argv[], string & PSSM_FileName, string & CutofsPerPSSM_FileName, string & Seq_FASTA_FileName, string & Hits_Out_FileName,
 						  size_t & numberOfRandomPSSM, bool &useFactor, bool & isOutputSequences, string & sequenceHitMotifPath) {
	// parse ARGV arguments
	size_t num_required_params = 8;
	//The first 5 params only - > 12 argc ,
	//The first 5 params with useFactor - > 13 argc
	//The first 5 params with isOutputSequences and sequenceHitMotifPath -> 15 argc
	//All params together -> 16 
	// each with its flag and mode_flag, check the value of argc. If not enough parameters have been passed, inform user and exit.
	if ((argc != (num_required_params + 4)) & (argc != (num_required_params + 5)) & (argc != (num_required_params * 2)-1) & (argc != (num_required_params * 2))) {
		cout << "Usage is -pssm <PSSMs_in_MAST_Format> -pssm_cutoffs <filename_for PSSM_cutoffs> -seq <input_seq_FASTA> -out <out>\n"; // Inform the user of how to use the program
		exit(12);
	}
	cout << argv[0] <<" ";
	for (int i = 1; i < argc; ++i) // iterate over flags
	{
		/* We will iterate over argv[] to get the parameters stored inside.
		* Note that we're starting on 1 because we don't need to know the
		* name of the program, which is stored in argv[0] */
		if (string(argv[i]) == "-pssm") PSSM_FileName = string(argv[i + 1]);
		else if (string(argv[i]) == "-pssm_cutoffs") CutofsPerPSSM_FileName = string(argv[i + 1]);
		else if (string(argv[i]) == "-seq") Seq_FASTA_FileName = string(argv[i + 1]);
		else if (string(argv[i]) == "-out") Hits_Out_FileName = string(argv[i + 1]);
		else if (string(argv[i]) == "-NrandPSSM") numberOfRandomPSSM = size_t(atoi(argv[i + 1]));
		else if (string (argv[i]) == "-useFactor") useFactor = true;
		else if (string (argv[i]) == "-outputSequences") isOutputSequences=true;
		else if (string (argv[i]) == "-sequenceHitMotifPath") sequenceHitMotifPath = string(argv[i + 1]);
		cout<<argv[i]<<" ";
	}
	if (numberOfRandomPSSM == 0) {
		cout<<"ERROR: could not parse the number of random PSSMs for Pvalue computation"<<endl;
		exit (99);
	}
	cout<<endl;
}

int computeCutoffsOfPssmMain(int argc, char *argv[])
{
	size_t TotalNumberOfRandoSeq=100000;
	string PSSM_FileName = "";
	string CutofsPerPSSM_FileName = "";
	int totalMemes = 0;
	double PercentOfRandomHitsPerPSSM = 0.05;
	getFileNamesFromArgv(argc, argv, PSSM_FileName, CutofsPerPSSM_FileName, totalMemes, PercentOfRandomHitsPerPSSM);
	readPSSM_info_from_file rpif(PSSM_FileName);
	if (totalMemes == 0) 
		totalMemes = rpif._PSSM_array.size();
	computePSSM_cutoffs cpc1(rpif._PSSM_array, TotalNumberOfRandoSeq, rpif._alph, CutofsPerPSSM_FileName, totalMemes, PercentOfRandomHitsPerPSSM);
	return 0;
}

int associatePeptidesWithPSSM(int argc, char *argv[])
{
	// -pssm <PSSMs_in_MAST_Format> -seq <input_seq_FASTA> -out <out> -calc_cutoff TRUE / FALSE -pssm_cutoffs <filename_for PSSM_cutoffs>
	
	// the function get 4 arguments from the command line. That is argc = 8.
	// 1. -pssm <PSSMs_in_MAST_Format>
	// 2. -pssm_cutoffs <filename_with_PSSM_cutoffs>
	
	string PSSM_FileName = "";
	string Seq_FASTA_FileName = "";
	string CutofsPerPSSM_FileName = "";
	string Hits_Out_FileName = "";

	getFileNamesFromArgv(argc, argv, PSSM_FileName, CutofsPerPSSM_FileName, Seq_FASTA_FileName, Hits_Out_FileName);
	readPSSM_info_from_file rpif(PSSM_FileName);
	rpif.update_PSSM_cutoff_from_file(CutofsPerPSSM_FileName);


	vector<SEQ> Seq_array;
	readFileToSeq_array(Seq_FASTA_FileName, rpif._alph, Seq_array, false);

	ofstream HitsReport_OutFileHandle;
	HitsReport_OutFileHandle.open(Hits_Out_FileName);
	for (size_t i = 0; i<rpif._PSSM_array.size(); ++i)
	{
		// vector <SEQ> hits;
		vector <HIT> hits;
		Find_PSSM_Hits(rpif._PSSM_array[i], Seq_array, hits, 1);
		// do something with hits...
		cout << "Finished to score peptides for PSSM #" << i << endl;
		print_PSSM_Hits(rpif._PSSM_array[i], hits, HitsReport_OutFileHandle);
		hits.clear();
	}

	HitsReport_OutFileHandle.close();
	return 0;
}

void print_PSSM_Hits(const PSSM & PSSM,const vector<HIT> & hits, ostream& out)
{
	SEQ tmp(PSSM._PSSM_consensus_seq, "", 1, PSSM._alph);
	out<<"### Total "<<hits.size()<<" hits of PSSM '"<<PSSM.PSSM_name<<"' with consensus sequence: "<< tmp.getStringOfSeq()<<" and max score: "<<PSSM.PSSM_MaxScore<<" are"<<endl;
	for (size_t i=0;i<hits.size();i++)
	{
		size_t PSSM_length = PSSM.PSSM_length();
		SEQ tmp = hits[i]._seq;
		string seqToPrint = tmp.getStringOfSeq();
		out << seqToPrint << "\t";
		out << hits[i]._seq._CopyNumber << "\t";
		out << hits[i]._seq._Seq_Name << "\t";
		out << hits[i]._match_score << "\t";
		out << hits[i]._match_pos << "\t";
		// build the match string
		string match_seq = "";
		if (hits[i]._match_pos > 0) // only part of the seq is matched -> the PSSM is shorther than the seq
		{
			match_seq = tmp.getStringOfSeq(hits[i]._match_pos);
		}
		else
		{
			if (hits[i]._match_pos < 0) // gaps to the left of match
			{
				size_t num_of_gaps_to_the_left = 0 - hits[i]._match_pos;
				match_seq.append(num_of_gaps_to_the_left, '-');
			}
			match_seq.append(seqToPrint); // the match itself
			if (match_seq.length() < PSSM_length) // gaps to the right of match
			{
				size_t num_of_gaps_to_the_right = PSSM_length - match_seq.length();
				match_seq.append(num_of_gaps_to_the_right, '-'); // add to end of sequence
			}
		}	
		// out << tmp.getStringOfSeq(hits[i]._match_pos)<<endl;
		// tmp.getStringOfSeq(hits[i]._match_pos) // TO DO: if only part of the seq => in case the PSSM is shorther than the seq
		out << match_seq << endl;
	}
}

// TO DO...
// void Find_PSSM_Hits(const PSSM& PSSM1, const vector<SEQ> & vseq1, vector<SEQ> & hits)

void Find_PSSM_Hits(const PSSM& PSSM1, const vector<SEQ> & vseq1, vector<HIT> & hits,const size_t verbose)
{
// get seq dataset
// scan with PSSMs and consider hits only the sequences getting score above cutoff
	size_t skipped_seq_counter = 0;
	for (size_t j =0; j < vseq1.size(); ++j) {
		//string best_match_string;
		int best_match_start=0;
		double score = PSSM1.computeScore(vseq1[j]._Seq,best_match_start);
		// cout << "score = " << score << endl;
		auto it = PSSM1.seq_type_cutoff.find(vseq1[j]._Seq_Type); // get an iterator pointing to the specific key
		double cutoff=DBL_MAX;
		if (it == PSSM1.seq_type_cutoff.end()) // key not found
		{
      if (verbose>0)
      {
      			cout << "WARNING: skipping sequence " << vseq1[j].getStringOfSeq() << " - can't find the appropriate cutoff for sequence type: "<< vseq1[j]._Seq_Type<<endl;
      }
			++skipped_seq_counter;
		}
		else
		{
			cutoff=it->second;  // get the actual sequence for the specific dataset by it index
		}
		if ((cutoff!=DBL_MAX)&&(score>=cutoff))
		{
			HIT hit(vseq1[j]);
			hit.set_hit_info(score,best_match_start);
			hits.push_back(hit);
			// hits.push_back(vseq1[j]);
			// hits[hits.size()-1].setHitScore(score);
		}
	}
	if (skipped_seq_counter > 0) {
		cout << "== Total skipped sequences due to not appropriate cutoff: " << skipped_seq_counter << endl;
		//cin.get();
	}
}

void PSSM_scoresFromSeqVector(const PSSM& PSSM1, const vector<vector<size_t> > & vseq1, vector<double> & scores) {
	for (size_t j =0; j < vseq1.size(); ++j) {
		int best_match_start=0;
		double score = PSSM1.computeScore(vseq1[j],best_match_start);
		scores.push_back(score);
	}
}

vector<size_t> PadSeq(const vector<size_t>& seq, size_t PaddingLength)
{
		vector<size_t> tempSeqPadded(PaddingLength,GAP);
		tempSeqPadded.insert(tempSeqPadded.end(), seq.begin(), seq.end());
		for (size_t i=0;i<PaddingLength;++i)
		{
			tempSeqPadded.push_back(GAP); // padding with '-' will get same score as unobserved charcter
		}
		return tempSeqPadded;
}

void readFileToSeq_array (const string fileName, alphabet& alph, vector<SEQ> &Seq_array, bool useCountDuplicate) {
	vector<string> allLines;
	size_t total_seq=0;
	fromFileToVectorOfString(fileName,allLines);
	for (size_t i=0; i<allLines.size(); ++i) {
		double count = 1;
		string currLine = allLines[i];
		string FirstChar=currLine.substr(0,1);
//		FirstChar=">";
		int compare=FirstChar.compare(">");
		if (currLine.substr(0,1).compare(">")==0) // header
		{
			string Seq;
			string name = currLine.substr(1);
			if (useCountDuplicate) {
				auto lastIndex = currLine.find_last_of('_');
            	count = stod(currLine.substr(lastIndex + 1));
			}
			++i;
			currLine = allLines[i];
			while (currLine.substr(0,1).compare(">")!=0)
			{
				Seq.append(currLine);
				++i;
				if (i<allLines.size())
				{
					currLine = allLines[i]; 
				}
				else
				{
					break;
				}
			}
			SEQ currSeq(Seq, name, count, alph);
			total_seq++;
			Seq_array.push_back(currSeq);
			i--; // got to new seq
		}
	}
	cout<<"Total Seq: "<<total_seq<<endl;
}

void fromFileToVectorOfString(const string fileName,vector<string> & allLines) {
	ifstream myfile;
	myfile.open(fileName.c_str());
	if (myfile.is_open())
	{
		string tmp;
		while (getline(myfile, tmp)) {
			allLines.push_back(tmp);
		}
		myfile.close();
	}
	else
	{
		cout << "Error: could not open the file " << fileName << endl;
exit(15);
	}

}

vector<SEQ> get_sort_seq_vector_by_scores(vector<SEQ> & seq_vector, vector<double> & scores_vector)
{

	// create corresponding index data
	vector<int> index(seq_vector.size(), 0);
	for (int i = 0; i != index.size(); i++) {
		index[i] = i;
	}

	// sort indexes by the scores value
	sort(index.begin(), index.end(),
		[&](const int& a, const int& b) {
		return (scores_vector[a] < scores_vector[b]);
	}
	);

	// get the data according to the sorted indices
	vector<SEQ> tmp_sorted_seq;
	vector<double> tmp_sorted_scores;

	for (int i = 0; i != index.size(); i++) {
		tmp_sorted_seq.push_back(seq_vector[index[i]]);
		tmp_sorted_scores.push_back(scores_vector[index[i]]);
		// cout << j << endl;
	}

	scores_vector = tmp_sorted_scores;
	return tmp_sorted_seq;
}

void fill_Seq_Hits_PSSM_map(const vector<SEQ> & seq_Hits_vector, map<string, vector<string>> &PSSM_Hits, string PSSM_name) {
	// get vector of hit sequences for specific PSSM and update the table mapping hit sequence to PSSM name
	for (size_t i = 0; i < seq_Hits_vector.size(); ++i) {
		string hit_sequence = seq_Hits_vector[i].getStringOfSeq();
		PSSM_Hits[hit_sequence].push_back(PSSM_name);
	}
}

void get_top_hits(const vector<SEQ> & sorted_seq, double fraction, vector <SEQ> & top_hits)
{
	double TopPercentSize = fraction*sorted_seq.size();
	int TopPercentIndex = int(sorted_seq.size() - floor(TopPercentSize));
	for (size_t i = TopPercentIndex; i < sorted_seq.size(); ++i)
	{
		top_hits.push_back(sorted_seq[i]);
	}
}

void writeSequenceHits(const PSSM& pssm1, const vector <HIT> & hits, const string sequenceHitsPath){
	ofstream fileSequenceHit;
	fileSequenceHit.open(sequenceHitsPath);
	fileSequenceHit << "MOTIF " << pssm1.PSSM_name << endl;
	cout<<hits[1]._seq.getStringOfSeq() <<endl;
    for (size_t k = 0; k < hits.size(); ++k) {
        fileSequenceHit << hits[k]._seq.getStringOfSeq() << " " << hits[k]._seq._CopyNumber << endl;
    }
	fileSequenceHit.close();
}


double numberOfTotalHitsPerPSSM(const PSSM& pssm1, const vector<SEQ> & Seq_array, const string sequenceHitsPath, const bool isOutputSequences, const size_t verbose=0) {
	vector <HIT> hits;
	Find_PSSM_Hits(pssm1, Seq_array, hits, verbose); //2 compute how many peptides are significant for this shuffled pssm
	double sum = 0;
	//print the sequence that had hit with the motif 
	if (isOutputSequences){
		writeSequenceHits(pssm1, hits, sequenceHitsPath);
	}
	for (size_t k = 0; k < hits.size(); ++k) {
		sum += hits[k]._seq._CopyNumber;
	}
	return sum;
}


int assignPvalueToPSSMaRRAY(int argc, char *argv[])
{
	// -pssm <PSSMs_in_MAST_Format> -seq <input_seq_FASTA> -out <out> -calc_cutoff TRUE / FALSE -pssm_cutoffs <filename_for PSSM_cutoffs>

	// the function get 4 arguments from the command line. That is argc = 8.
	// 1. -pssm <PSSMs_in_MAST_Format>
	// 2. -pssm_cutoffs <filename_with_PSSM_cutoffs>

	string PSSM_FileName = "";
	string Seq_FASTA_FileName = "";
	string CutofsPerPSSM_FileName = "";
	string Hits_Out_FileName = "";
	size_t numberOfRandomPSSM = 0;
	bool useFactor = false;
	bool isOutputSequences = false;
	string sequenceHitMotifPath = "";
	getFileNamesFromArgv(argc, argv, PSSM_FileName, CutofsPerPSSM_FileName, Seq_FASTA_FileName, Hits_Out_FileName, numberOfRandomPSSM,
						 useFactor, isOutputSequences, sequenceHitMotifPath);
	cout<<"Number of random PSSMs to calculate Pval: "<< numberOfRandomPSSM <<endl;
	readPSSM_info_from_file rpif(PSSM_FileName);
	rpif.update_PSSM_cutoff_from_file(CutofsPerPSSM_FileName);

	vector<SEQ> Seq_array;
	readFileToSeq_array(Seq_FASTA_FileName, rpif._alph, Seq_array, true);

	ofstream listOfPvaluesFile;
	listOfPvaluesFile.open(Hits_Out_FileName);
	listOfPvaluesFile << "## PSSM_name\tp_Value\tTrue_Hits: num_of_hits" <<endl;
	for (size_t i = 0; i < rpif._PSSM_array.size(); ++i) {
	//for (size_t i = 0; i < 1; ++i) {
		double numberOfHitsInRealPSSM = numberOfTotalHitsPerPSSM(rpif._PSSM_array[i], Seq_array, sequenceHitMotifPath, isOutputSequences, 1);
		vector<double> numSigPeptides;
		default_random_engine gen(483); // TODO seed should be fro input
		for (size_t j = 0; j < numberOfRandomPSSM; ++j) {
			PSSM randomPSSM = rpif._PSSM_array[i].randomize(gen); //1 generate a random PPSM.
			double sum = numberOfTotalHitsPerPSSM(randomPSSM, Seq_array, sequenceHitMotifPath, false, 0);
			numSigPeptides.push_back(sum);// store the number
		}
		sort(numSigPeptides.begin(), numSigPeptides.end());
		// for (size_t p = 0; p < numSigPeptides.size(); ++p) {
		//	listOfPvaluesFile << numSigPeptides[p] << " ";
		// }
		// listOfPvaluesFile << endl;

		int place = numberOfRandomPSSM-1;
		while (place>=0) {
			if (numberOfHitsInRealPSSM > numSigPeptides[place]) break;
			place--;
		}
		if (place == -1) place = 0; //so that we get p value = 1 in this case.
		//cout << "place = " << place << endl;
		double p_Value = (numberOfRandomPSSM - place +0.0) / numberOfRandomPSSM;
		
		if (useFactor) {
			int numberOfSeq = 0;
			for (size_t k = 0; k < Seq_array.size(); ++k) {
				numberOfSeq += Seq_array[k]._CopyNumber;
			}
			float factor = float(1000000) / float(numberOfSeq);
			numberOfHitsInRealPSSM *= factor;
		}

		listOfPvaluesFile << rpif._PSSM_array[i].PSSM_name << "\t" << p_Value << "\tTrue_Hits: " << numberOfHitsInRealPSSM <<endl; // << " total true hits " << numberOfHitsInRealPSSM << endl;
		cout << "finished with PSSM " << i << endl;
	}
	listOfPvaluesFile.close();
	return 0;
}
