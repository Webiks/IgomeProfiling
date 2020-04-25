#include <string>
#include <iostream>
#include <math.h>
#include <fstream>
#include <regex>
#include <chrono>

#include "cxxopts.hpp"
#include "types.hpp"
#include "trim.hpp"
#include "meme.hpp"
#include "memes.hpp"
#include "shufflesGenerator.hpp"

Memes loadMemes(string memePath, int limit, bool verbose);
void loadCutoffs(string cutoffsPath, Memes& memes, int limit, bool verbose);
SequencesMap loadSequences(string faaPath, bool verbose);

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
    int printInterval, bool isOutputSequences, bool verbose) {
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
                meme.addHitSequence(*sequencesIter, isOutputSequences);
                hits++;
            }
            sequencesIter++;
        }
        sequencesTypesIter++;
    }
    cout << "meme hits: " << meme.getHitCount() << endl;
}

int getHits(Memes& memes, SequencesMap& sequences, MemeShufflesMap& shuffles, bool isOutputSequences, bool verbose) {
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
            printInterval, isOutputSequences, verbose);
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
                    printInterval, isOutputSequences, verbose);
                shufflesIter++;
            }            
        }

        memesIter++;
    }
    cout << "total hits: " << hits << endl;
    return hits;
}

MemeRatingMap getRatings(Memes& memes, MemeShufflesMap& shuffles, bool isOutputSequences, bool verbose) {
    MemeRatingMap ratings;

    auto memeIter = shuffles.begin();
    auto memeEnd = shuffles.end();

    while (memeIter != memeEnd) {
        auto hits = memes.getMemes()[memeIter->first].getHitCount();
        int rank = 0;

        auto shuffleIter = memeIter->second.begin();
        auto shuffleEnd = memeIter->second.end();
        while (shuffleIter != shuffleEnd) {
            if (hits > shuffleIter->getHitCount()) { // TODO should be >= ? 
                rank++;
            }
            shuffleIter++;
        }
        auto score = (float)rank / memeIter->second.size();
        if (verbose) {
            cout << "Motif " << memeIter->first << ", Total shuffles: " << memeIter->second.size() <<
                ", Rank: " << rank << ", Score: " << score << endl;
        }
        ratings[memeIter->first] = score;
        memeIter++;
    }

    return ratings;
}

void writeResults(Memes& memes, MemeRatingMap& ratings, string& outputPath, bool verbose) {
    auto memesIter = memes.getMemes().begin();
    auto memesEnd = memes.getMemes().end();
    auto ratingEnd = ratings.end();
    ofstream file(outputPath);

    while (memesIter != memesEnd) {
        file << "MOTIF " << memesIter->second.getMotif() << endl;
        file << "HITS " << memesIter->second.getHitCount() << endl;
        auto ratingIter = ratings.find(memesIter->first);
        if (ratingIter != ratingEnd) {
            file << "RANK " << ratingIter->second << endl;
        }
        auto sequencesIter = memesIter->second.getHitSequences().begin();
        auto sequencesEnd = memesIter->second.getHitSequences().end();
        while (sequencesIter != sequencesEnd) {
            file << sequencesIter->first << " " << sequencesIter->second << endl;
            sequencesIter++;
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
        //shufflesintersections
        ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"));
    auto result = options.parse(argc, argv);

    auto memesPath = result["memes"].as<string>();
    auto cutoffsPath = result["cutoffs"].as<string>();
    auto sequencesPath = result["sequences"].as<string>();
    auto outputPath = result["output"].as<string>();
    auto isOutputSequences = result["outputSequences"].as<bool>();
    auto maxMemes = result["maxMemes"].as<int>();
    auto shuffles = result["shuffles"].as<int>();
    auto isVerbose = result["verbose"].as<bool>();

    auto begin = chrono::steady_clock::now();

    auto memes = loadMemes(memesPath, maxMemes, isVerbose);
    loadCutoffs(cutoffsPath, memes, maxMemes, isVerbose);
    SequencesMap sequences = loadSequences(sequencesPath, isVerbose);
    auto memesShuffles = createShuffles(memes, shuffles);
    getHits(memes, sequences, memesShuffles, isOutputSequences, isVerbose);
    auto memesRating = getRatings(memes, memesShuffles, isOutputSequences, isVerbose);
    writeResults(memes, memesRating, outputPath, isVerbose);

    auto end = chrono::steady_clock::now();
    cout << chrono::duration_cast<chrono::seconds>(end - begin).count() << "[s]" << endl;
}
