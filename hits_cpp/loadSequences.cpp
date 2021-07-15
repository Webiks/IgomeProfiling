#include <iostream>
#include <string>
#include <fstream>
#include <regex>

#include "types.hpp"
#include "trim.hpp"

using namespace std;
SequencesMap loadSequences(string faaPath, int& numSequences, bool verbose) {
    SequencesMap sequences;
    ifstream file(faaPath);
    string line;
    list<string>* sequencesByType;
    
    regex pattern("^>.+_(.+)$");
    smatch matches;
    auto end = sequences.end();

    int count = 0;
    while (getline(file, line)) {
        if (line[0] == '>') {
            auto lastIndex = line.find_last_of('_');
            auto seqType = line.substr(lastIndex + 1);
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
            count++;
        }
    }
    numSequences = count;
    cout << "total sequences: " << count << endl;
    return sequences;
}
