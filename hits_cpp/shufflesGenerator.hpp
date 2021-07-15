#pragma once

#include "meme.hpp"
#include <algorithm>
using namespace std;

class MemeShuffler {
public:
    MemeShuffler(Meme& source, int maxShuffles, ShufflePatterns& patterns) :
        _source(source), _patterns(patterns), _currentShuffle(-1), _lenPatterns(patterns.size()), _numPattern(-1) {
        auto len = patterns.size();
        if (len < maxShuffles) {
            maxShuffles = len;
        }
        this->_maxShuffles = maxShuffles;
    }

    bool next() {
        _currentShuffle++;
        _numPattern++;
        return hasNext();
    }

    Meme generate() {
        auto meme = Meme(this->_source);
        meme.getRows().clear();
        auto pattern = this->_patterns[this->_numPattern];
        auto iter = pattern.begin();
        auto end = pattern.end();
        while (iter != end) {
            meme.getRows().push_back(this->_source.getRows()[*iter]);
            iter++;
        }
        //check that we don't change place of the same consensus:
        for (int i = 0; i < this->_source.getRows().size(); i++){
            auto max_source = max_element(this->_source.getRows()[i].begin(), this->_source.getRows()[i].end()) - this->_source.getRows()[i].begin();
            auto max_shuffle = max_element(meme.getRows()[i].begin(), meme.getRows()[i].end()) - meme.getRows()[i].begin();
            if (max_source == max_shuffle && (this->_lenPatterns - this->_numPattern) > (this->_maxShuffles - this->_currentShuffle)){
                _numPattern++;
                return this->generate();
            }
        }  
        return meme;
    }

    bool hasNext() {
        return _currentShuffle < _maxShuffles;
    }

    int getMaxShuffles() {
        return this->_maxShuffles;
    }
private:
    Meme _source; //TODO pointer?
    ShufflePatterns _patterns; //TODO pointer?
    int _maxShuffles;
    int _currentShuffle;
    int _lenPatterns;
    int _numPattern;
};

MemeShuffler getShuffler(Meme& meme, int maxPatterns);
