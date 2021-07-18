#include "SEQ.h"
/// NOTE: to suppurt regex, gcc>=4.9 is needed!!!
#define min_seq_length 6
#define max_seq_length 14

void SEQ::setCopyNmber()
{
	_CopyNumber = 1.0;
	smatch match;
	regex rgx(".*Repeats_([\\d.]+)_Type");
	if (regex_search(_Seq_Name, match, rgx))
	{
		string CopyNumber_str = match[1].str();
		_CopyNumber = atof(CopyNumber_str.c_str());
	}
}
void SEQ::setSeq_Type()
{
	int seq_length = _Seq.size();
	stringstream ss;
	if (seq_length < 1)
	{
		cout << "Invalid sequence " << _Seq_Name << " ignored..." << endl;
	}
	char first = _Seq[0];
	char last = _Seq[seq_length - 1];
	//if ((first == 'C') && (last == 'C'))
	if ((first == _alph._alphabetMap['C']) && (last == _alph._alphabetMap['C']) && (seq_length > min_seq_length)) // if not long enough don't account for CysLoop - e.g. C4C=6
	{
		seq_length = seq_length - 2;
		ss << "C" << seq_length << "C";
	}
	else
	{
		ss << seq_length;
	}
	_Seq_Type = ss.str();
}

string SEQ::getStringOfSeq() const {
	string res = "";
	for (size_t i = 0; i < _Seq.size(); ++i) {
		res.append(1,_alph._correspondingCharacters[_Seq[i]]);
	}
	return res;
}

string SEQ::getStringOfSeq(size_t pos) const {
	string res = "";
	for (size_t i = pos; i < _Seq.size(); ++i) {
		res.append(1, _alph._correspondingCharacters[_Seq[i]]);
	}
	return res;
}

SEQ::SEQ(string & SeqString, const string & SeqName, const double CopyNumber, alphabet& alph) :  _Seq_Name(SeqName), _CopyNumber(CopyNumber), _alph(alph) {
	for (size_t k = 0; k < SeqString.length(); ++k) {
		char i = SeqString.at(k);
		_Seq.push_back(alph._alphabetMap[i]);
	}
	setSeq_Type();
	//setCopyNmber();
	
};
