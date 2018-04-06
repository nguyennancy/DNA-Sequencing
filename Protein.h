#ifndef PROTEIN_H_
#define PROTEIN_H_

#include "AminoAcidChain.h"


class Protein
{
  friend class AminoAcidChain; // to access AAList in main

public:
  AminoAcidChain chain1;  // first reading frame
  AminoAcidChain chain1c; // complementary strand first reading frame
  AminoAcidChain chain2;  // second reading frame; first letter will be truncated
  AminoAcidChain chain2c; // complementary strand second reading frame
  AminoAcidChain chain3;  // third reading frame; first 2 letters will be truncated
  AminoAcidChain chain3c; // complementary strand third reading frame

  void genAminoAcidSequences(const LinearHashTable &AAList);
  void genOtherStrands();
  friend istream& operator>> (istream &inf, Protein &p);
  void printAAs() const;
  void printCDNAs() const;
  void printmRNAs() const;
  void printSeqs() const;
  void printStructures() const;
  AminoAcidChain& returnLongestChain();
  int returnLongestChainNum() const;
  
}; // class Protein

#endif // PROTEIN_H_