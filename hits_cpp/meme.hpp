#pragma once
#include "types.hpp"

class Meme {
public:
    Meme() : _hitCount(0) {

    }

    Meme(const Meme& other) :
        _motif(other._motif),
        _alength(other._alength),
        _nsites(other._nsites),
        _hitCount(other._hitCount),
        _rows(other._rows),
        _cutoffs(other._cutoffs),
        _hitSequences(other._hitSequences) {
    }

    string getMotif() {
        return this->_motif;
    }

    void setMotif(string motif) {
        this->_motif = motif;
    }

    int getALength() {
        return this->_alength;
    }

    void setALength(int alength) {
        this->_alength = alength;
    }

    int getNSites() {
        return this->_nsites;
    }

    void setNSites(int nsites) {
        this->_nsites = nsites;
    }

    MemeRows& getRows() {
        return this->_rows;
    }

    CutoffsMap& getCuttofs() {
        return this->_cutoffs;
    }

    SequencesCount& getHitSequences() {
        return this->_hitSequences;
    }

    int getHitCount() {
        return this->_hitCount;
    }

    void factorHitCount(float factor) {
        this->_hitCount *=factor;
    }

    void addHitSequence(string& sequence, bool isStoreSequences = true, int count = 1) {
        this->_hitCount += count;
        if (isStoreSequences) {
            if (this->_hitSequences.find(sequence) == this->_hitSequences.end()) {
                this->_hitSequences[sequence] = count; 
            } else {
                this->_hitSequences[sequence] += count;
            }
        }
    }

    void normalize() {
        if (this->_nsites == 0) {
            return;
        }
        auto rows = this->_rows.size();
        auto columns = this->_alength + 1;
        
        auto iter = this->_rows.begin();
        auto end = this->_rows.end();

        for (auto rowsIter = this->_rows.begin(); rowsIter != this->_rows.end(); ++rowsIter) {
            for (auto columnsIter = rowsIter->begin(); columnsIter != rowsIter->end(); ++columnsIter) {
                *columnsIter = ((*columnsIter * this->_nsites) + 1) / (this->_nsites + columns);
            }
        }
    }
private:
    string _motif;
    int _alength;
    int _nsites;
    int _hitCount;
    MemeRows _rows;
    CutoffsMap _cutoffs;
    SequencesCount _hitSequences;
};
