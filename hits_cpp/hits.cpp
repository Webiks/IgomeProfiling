#include <string>
#include <iostream>
#include <math.h>
#include <fstream>
#include <regex>
#include <chrono>
#include <iomanip>

#include "cxxopts.hpp"
#include "types.hpp"
#include "trim.hpp"
#include "meme.hpp"
#include "memes.hpp"
#include "shufflesGenerator.hpp"

Memes loadMemes(string memePath, int limit, bool verbose);
void loadCutoffs(string cutoffsPath, Memes& memes, int limit, bool verbose);
SequencesMap loadSequences(string faaPath, int& numSequences, SequencesRpmMap& seeuncesRpm, bool verbose);

// TODO support Repeats_
// TODO move createShuffles, isHit and getHits?
MemeShufflesMap createShuffles(Memes& memes, int shuffles) {
    MemeShufflesMap memesShuffles;

    if (shuffles) {
        // Get keys (before mutation)
        vector<string> keys;
        auto memesIter = memes.getMemes().begin();
        auto memesEnd = memes.getMemes().end();
        while (memesIter != memesEnd) {
            keys.push_back(memesIter->first);
            memesIter++;
        }

        // Add shuffles
        auto keysIter = keys.begin();
        auto keysEnd = keys.end();
        while (keysIter != keysEnd) {
            auto shuffler = getShuffler(memes.getMemes()[*keysIter], shuffles);
            while (shuffler.next()) {
                memesShuffles[*keysIter].push_back(shuffler.generate());
            }
            keysIter++;
        }
    }

    return memesShuffles;
}

bool isHit(Meme& meme, AlphabetMap& alphabet, string seqType, string& seq, bool verbose) {
    auto iter = meme.getCuttofs().find(seqType);
    if (iter == meme.getCuttofs().end()) {
        return false;
    }
    auto cutoffValue = iter->second;
    auto seqLen = seq.length();
    auto memeLen = meme.getRows().size();

    int start = -memeLen + 1;
    int end = min(memeLen, seqLen);
    double totalScore = 0;
    double score = 0;
    char c;
    for (int i = start; i < end; i++) {
        totalScore = 1;
        for (int j = 0; j < memeLen; j++) {
            c = ((i + j < 0) || (i + j) >= seqLen) ? '-' : seq[i + j];
            score = meme.getRows()[j][alphabet[c]];
            // TODO check if equal to -inf, set total_score and break
            totalScore *= score;
        }
        if (log(totalScore) > cutoffValue) {
            return true;
        }
    }
    return false;
}

void memeHits(Meme& meme, AlphabetMap& alphabet, SequencesMap& sequences, int& hits, 
    int printInterval, bool isOutputSequences, SequencesRpmMap& sequncesRpm, bool isUseUniqueRpm, bool verbose) {
    auto sequencesTypesIter = sequences.begin();
    auto sequencesTypesEnd = sequences.end();
    int counter = 0;

    while (sequencesTypesIter != sequencesTypesEnd) {
        auto sequencesIter = sequencesTypesIter->second->begin();
        auto sequencesEnd = sequencesTypesIter->second->end();
        while (sequencesIter != sequencesEnd) {
            if (counter && (counter % printInterval) == 0) {
                cout << "seq: " << counter << ", overall hits: " << hits << endl;
            }
            counter++;
            if (isHit(meme, alphabet, sequencesTypesIter->first, *sequencesIter, verbose)) {
                if (isUseUniqueRpm){
                    meme.addHitSequence(*sequencesIter, isOutputSequences, sequncesRpm.find(*sequencesIter)->second);
                } else {
                    meme.addHitSequence(*sequencesIter, isOutputSequences);
                }
                hits += sequncesRpm.find(*sequencesIter)->second;
            }
            sequencesIter++;
        }
        sequencesTypesIter++;
    }
    cout << "meme hits: " << meme.getHitCount() << endl;
}

void writeSequenceHits(SequencesCount& hitSequences,  SequencesRpmMap& sequncesRpm, string motif, string sequenceHitMotifPath){
    ofstream fileSequenceHit;
    fileSequenceHit.open(sequenceHitMotifPath);
    fileSequenceHit << "MOTIF " << motif << endl;
    auto sequencesTypesIter = hitSequences.begin();
    auto sequencesTypesEnd = hitSequences.end();
    while (sequencesTypesIter != sequencesTypesEnd){
        fileSequenceHit << sequencesTypesIter->first << " " << sequncesRpm.find(sequencesTypesIter->first)->second << endl;
        sequencesTypesIter++;
    }
    fileSequenceHit.close();
}

int getHits(Memes& memes, SequencesMap& sequences, MemeShufflesMap& shuffles, bool isOutputSequences, string sequenceHitMotifPath, SequencesRpmMap& sequncesRpm, bool verbose) {
    if (verbose) {
        cout << "GET HITS" << endl;
    }
    auto alphabet = memes.getAlphabet();
    auto memesIter = memes.getMemes().begin();
    auto memesEnd = memes.getMemes().end();
    int hits = 0;
    int shuffleHits = 0;
    int counter = 0;
    int printInterval = 100000;
    if (verbose) {
        printInterval = 10000;
    }

    while (memesIter != memesEnd) {
        if (verbose) {
            cout << "Calculating hits for " << memesIter->first << endl;
        }
        memeHits(memesIter->second, alphabet, sequences, hits, 
            printInterval, isOutputSequences, sequncesRpm, true, verbose);
        //print all the sequences that have hit with the motif
        if (isOutputSequences){
            SequencesCount hitSequences = (memesIter->second).getHitSequences();
            writeSequenceHits(hitSequences, sequncesRpm, memesIter->first, sequenceHitMotifPath);
        }
        auto memeShuffles = &shuffles[memesIter->first];
        if (memeShuffles->size()) {
            counter = 0;
            auto shufflesIter = memeShuffles->begin();
            auto shufflesEnd = memeShuffles->end();
            while (shufflesIter != shufflesEnd) {
                if (verbose) {
                    cout << "Calculating hits for shuffle " << ++counter << endl;
                }
                memeHits(*shufflesIter, alphabet, sequences, shuffleHits,
                    printInterval, isOutputSequences, sequncesRpm, false, verbose);
                shufflesIter++;
            }            
        }

        memesIter++;
    }
    cout << "total hits: " << hits << endl;
    return hits;
}

MemeRatingMap getRatings(Memes& memes, MemeShufflesMap& shuffles, bool verbose, float shufflesPercent) {
    MemeRatingMap ratings;

    auto memeIter = shuffles.begin();
    auto memeEnd = shuffles.end();

    while (memeIter != memeEnd) {
        auto hits = memes.getMemes()[memeIter->first].getHitCount();

        auto shuffleIter = memeIter->second.begin();
        auto shuffleEnd = memeIter->second.end();
 
        int add = memeIter->second.size();
        if (shufflesPercent != 0.0) {
            add++;
        }
        int shuffleArray[add];
        int count = 0;
        int max = 0;
        //create list of the suffle hits
        while (shuffleIter != shuffleEnd) {
                shuffleArray[count] = shuffleIter->getHitCount();
                if (shuffleIter->getHitCount() > max){
                    max=shuffleIter->getHitCount();
                }
                shuffleIter++;
                count++;
        }
        if (add != memeIter->second.size()) {
            shuffleArray[memeIter->second.size()] = (int)((float)max * shufflesPercent) + max;
        }
        //sort the list 
        int n = sizeof(shuffleArray) / sizeof(shuffleArray[0]);
        sort(shuffleArray, shuffleArray + n);
        
        int shuffleHitStart = shuffleArray[0];
        int shuffleHitEnd = shuffleArray[n-1];
        float score;
        if (hits == 0) {
            score = 0.0;
        }
        else if (hits > shuffleHitEnd) {
            score = 1.0;
        }
        else if (hits < shuffleHitStart) {
            score = 0.0;
        }else{
            for (int place = 0; place < n ;place++) {
                if (hits <=shuffleArray[place]) {
                    float min_1 = (float)shuffleArray[place - 1];
                    float max_1 = (float)shuffleArray[place];
                    float z_maen = ((float)hits - min_1) / (max_1 - min_1);
                    float p = (float)place / (float)n;
                    score = (z_maen / (float)n) + p;
                    break;
                }
			} 
        }
        if (verbose) {
            cout << "Motif " << memeIter->first << ", Total shuffles: " << memeIter->second.size() <<
                ", Score: " << score << endl;
        }
        ratings[memeIter->first] = score;
        memeIter++;
    }

    return ratings;
}

float getRpmFactor(int numSequences) {
    float factor;
    if (numSequences != 0) {
        factor = 1000000 / numSequences;
    }
    return factor;
}

void factorHits(Memes& memes, float factor) {
    auto memesIter = memes.getMemes().begin();
    auto memesEnd = memes.getMemes().end();
    while (memesIter != memesEnd) {
        memesIter->second.factorHitCount(factor);
        memesIter++;
    }
}


void writeResults(Memes& memes, MemeRatingMap& ratings, MemeShufflesMap& shuffles, string& outputPath, bool verbose, int shufflesDigits) {
    auto memesIter = memes.getMemes().begin();
    auto memesEnd = memes.getMemes().end();
    auto ratingEnd = ratings.end();
    ofstream file(outputPath);

    while (memesIter != memesEnd) {
        file << "MOTIF " << memesIter->second.getMotif() << endl;
        file << "HITS " << memesIter->second.getHitCount() << endl;
        auto ratingIter = ratings.find(memesIter->first);
        if (ratingIter != ratingEnd) {
            file << "SHUFFLES " << shuffles[memesIter->first].size() << endl;
            file << "RANK " << std::fixed << std::setprecision(shufflesDigits) <<ratingIter->second << endl;
        }
        memesIter++;
    }
}

int main(int argc, char *argv[])
{
    cxxopts::Options options("Hits", "Calculate hits given mime, cutoffs and sequences");
    options.add_options()
        ("m,memes", "Path to memes file", cxxopts::value<string>())
        ("c,cutoffs", "Path to cutoffs file", cxxopts::value<string>())
        ("s,sequences", "Path to sequences file", cxxopts::value<string>())
        ("o,output", "Path to results file", cxxopts::value<string>())
        ("maxMemes", "Limit number of memes to process (0 = all)", cxxopts::value<int>()->default_value("0"))
        ("outputSequences", "Write matched sequences (not memory efficient)", cxxopts::value<bool>()->default_value("false"))
        ("shuffles", "Create shuffles and rate memes by them (0 = disable)", cxxopts::value<int>()->default_value("0"))
        ("shufflesPercent", "Percent from shuffle with greatest number of hits (0-1)", cxxopts::value<float>()->default_value("0.2"))
        ("shufflesDigits", "Number of digits after the point to print in scanning files", cxxopts::value<int>()->default_value("2"))
        ("useFactor", "To multiply by factor hits for normalization", cxxopts::value<bool>()->default_value("false"))
        ("sequenceHitMotifPath", "Path for results of sequence that had hit with motif", cxxopts::value<string>()->default_value(""))
        ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"));
    auto result = options.parse(argc, argv);

    auto memesPath = result["memes"].as<string>();
    auto cutoffsPath = result["cutoffs"].as<string>();
    auto sequencesPath = result["sequences"].as<string>();
    auto outputPath = result["output"].as<string>();
    auto isOutputSequences = result["outputSequences"].as<bool>();
    auto maxMemes = result["maxMemes"].as<int>();
    auto shuffles = result["shuffles"].as<int>();
    auto shufflesPercent = result["shufflesPercent"].as<float>();
    auto shufflesDigits = result["shufflesDigits"].as<int>();
    auto useFactor = result["useFactor"].as<bool>();
    auto sequenceHitMotifPath = result["sequenceHitMotifPath"].as<string>();
    auto isVerbose = result["verbose"].as<bool>();

    auto begin = chrono::steady_clock::now();
    int numSequences;
    SequencesRpmMap sequncesRpm;
    auto memes = loadMemes(memesPath, maxMemes, isVerbose);
    loadCutoffs(cutoffsPath, memes, maxMemes, isVerbose);
    SequencesMap sequences = loadSequences(sequencesPath, numSequences, sequncesRpm, isVerbose);
    auto memesShuffles = createShuffles(memes, shuffles);
    getHits(memes, sequences, memesShuffles, isOutputSequences, sequenceHitMotifPath, sequncesRpm ,isVerbose);
    MemeRatingMap memesRating;
    if (shuffles) {
        memesRating = getRatings(memes, memesShuffles, isVerbose, shufflesPercent);
    }
    if (useFactor) {
        float factor = getRpmFactor(numSequences);
        factorHits(memes, factor);
    }
    writeResults(memes, memesRating, memesShuffles, outputPath, isVerbose, shufflesDigits);

    auto end = chrono::steady_clock::now();
    cout << chrono::duration_cast<chrono::seconds>(end - begin).count() << "[s]" << endl;
}
