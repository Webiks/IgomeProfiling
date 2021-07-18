#include <iostream>
#include <string>
#include <fstream>
#include <regex>
#include <sstream>
#include <vector>

#include "types.hpp"
#include "trim.hpp"

using namespace std;
SequencesMap loadSequences(string faaPath, int& numSequences, SequencesRpmMap& sequncesRpm, bool verbose) {
    SequencesMap sequences;
    ifstream file(faaPath);
    string line;
    list<string>* sequencesByType;
    
    regex pattern("^>.+_(.+)$");
    smatch matches;
    auto end = sequences.end();
    int countAllSequrnces = 0;
    int count = 0;
    while (getline(file, line)) {
        if (line[0] == '>') {
            replace(line.begin(), line.end(), '_', ' ');  // replace '_' by ' '
            istringstream iss(line);
            vector<string> result;
            for(string s;iss>>s;){
                result.push_back(s);
            }
            auto seqType = result[3];
            int unique_rpm = stoi(result[7]);
            rtrim(seqType);
        
            auto iter = sequences.find(seqType);
            if (iter == end) {
                sequencesByType = new list<string>();
                sequences[seqType] = sequencesByType;
            } else {
                sequencesByType = iter->second;
            }
            
            getline(file, line);
            sequencesByType->push_back(line);
            sequncesRpm[line] = unique_rpm;
            count++ ;
            countAllSequrnces += unique_rpm;
        }
    }
    numSequences = countAllSequrnces;
    cout << "total sequences: " << count << endl;
    return sequences;
}
