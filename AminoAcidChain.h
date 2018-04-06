#ifndef AMINOACIDCHAIN_H_
#define AMINOACIDCHAIN_H_

#include <cstdlib>
#include <cmath>

#include "LinearProbing.h"
#include "AminoAcid.h"
#include "Structure.h"
#include "Sequence.h"

using namespace std;

class AminoAcidChain
{
  friend class LinearHashTable;

public:
  Sequence sequence;    // contains DNA sequence + mRNA sequences
  int metPos;           // position in sequence at start codon
  int stopPos;          // position in sequence at stop codon
  AminoAcid** AASeq;    // amino acid sequence of given
  string aaSeq;         // string version of amino acid sequence
  int aaLength;         // length of amino acid sequence
  Structure* structure; // structural formula of amino acid chain

  AminoAcidChain();
  int getaaLength() const { return aaLength; }
  AminoAcid** getAASeq() const { return AASeq; }
  void genAminoAcidSeq(const LinearHashTable &AALst);
  void genComplement(const AminoAcidChain &a);
  void genReadingFrames(AminoAcidChain &c2, AminoAcidChain &c3) const;
  void genStructures();
  friend istream& operator>> (istream &inf, AminoAcidChain &a);
  void printAA() const;
  void printAll() const;
  void printCDNA() const;
  void printmRNA() const;
  void printSeq() const;
  void printStructure() const;
  void setAAChainMembers();

}; // class AminoAcidChain

#endif // AMINOACIDCHAIN_H_