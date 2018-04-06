#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include "LinearProbing.h"
#include "AminoAcid.h"

using namespace std;

class Sequence
{
  friend class LinearHashTable;

  string seq;     // nucleotide sequence of given strand from 5' to 3'
  string seqMRNA; // full mRNA generated from seq, including non-coding regions
  string mRNA;    // mRNA of coding region with polyATail
  int seqLength;  // length of full nucelotide sequence

  int compareThreeNums(int num1, int num2, int num3) const;
  int findFirstStopCodon(const string &s, size_t findStop, int offset) const;
  int findCodon(const string &s, const string &codon, size_t found, int offset) const;
  char genRandomBase() const;

public:
  Sequence();
  int calcNumMutations(const string &orig, const string &mut) const;
  void findRepeats() const;
  void findPairRepeats(vector<string> &repeats, vector<int> &counts, 
                                 vector<int> &starts) const;
  void findSingleRepeats(vector<string> &repeats, vector<int> &counts, 
                                 vector<int> &starts) const;
  void genComplementSeq(const Sequence &s);
  void genMRNA(int aLength, int startPos, int stopPos);
  string genMutatedAminoAcidSeq(const LinearHashTable &AALst);
  void genMutatedSeq(const LinearHashTable &AALst, const Sequence &s,
                     int percent, int start, int stop);
  void genReadingFrameSeqs(Sequence &s2, Sequence &s3) const;
  void genSeqMRNA();
  string getCDNA(int start, int stop) const { return seq.substr(start, stop); }
  string getMRNA() const { return mRNA; }
  string getSeq() const { return seq; }
  int getSeqLength() const { return seqLength; }
  string getSeqMRNA() const { return seqMRNA; }
  friend istream& operator>> (istream &inf, Sequence &s);
  void printRepeats(vector<string> repeats, vector<int> counts, 
                    vector<int> starts) const;
  void searchForLongestChain(const string &s, int &start, int &stop) const;
  void setSequenceMembers();

}; // class Sequence

#endif // SEQUENCE_H_