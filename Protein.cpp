#include <cstdlib> // for srand and rand
#include <ctime>   // for time in srand
#include "Protein.h"


// generates the amino acid sequences/chains + structures
void Protein::genAminoAcidSequences(const LinearHashTable &AAList)
{
  chain1.genAminoAcidSeq(AAList); chain1.genStructures();
  chain1c.genAminoAcidSeq(AAList); chain1c.genStructures();

  chain2.genAminoAcidSeq(AAList); chain2.genStructures();
  chain2c.genAminoAcidSeq(AAList); chain2c.genStructures();

  chain3.genAminoAcidSeq(AAList); chain3.genStructures();
  chain3c.genAminoAcidSeq(AAList); chain3c.genStructures();
} // genAminoAcidSequences()


// generates the other reading frames + complement nucelotide sequences and
// stores them from 5' to 3' + set values for some data members
void Protein::genOtherStrands()
{
  chain1c.genComplement(chain1); // Reading frame 1 complementary
  chain1c.setAAChainMembers(); 

  // Reading frames 2 and 3 orig
  chain1.genReadingFrames(chain2, chain3);
  chain1c.genReadingFrames(chain2c, chain3c);
} // genOtherStrands


istream& operator>> (istream &inf, Protein &p)
{
  inf >> p.chain1;

  p.genOtherStrands();
  return inf;
} // operator>> ()


void Protein::printAAs() const
{
  cout << "Reading frame 1" << endl;
  cout << "Original amino acid strand:" << endl << endl;
  chain1.printAA(); 
  
  cout << "Complementary amino acid strand:" << endl << endl;
  chain1c.printAA();

  cout << endl << endl << "Reading frame 2" << endl;
  cout << "Original amino acid strand:" << endl << endl;
  chain2.printAA(); 

  cout << "Complementary amino acid strand:" << endl << endl;
  chain2c.printAA();

  cout << endl << endl << "Reading frame 3" << endl;
  cout << "Original amino acid strand:" << endl << endl;
  chain3.printAA(); 

  cout << "Complementary amino acid strand:" << endl << endl;
  chain3c.printAA();
} // printAAs()


void Protein::printCDNAs() const
{
  cout << "Reading frame 1" << endl;
  cout << "Original cDNA strand:" << endl;
  chain1.printCDNA(); 

  cout << endl << "Complementary cDNA strand:" << endl;
  chain1c.printCDNA(); 

  cout << endl << endl << "Reading frame 2" << endl;
  cout << "Original cDNA strand:" << endl;
  chain2.printCDNA(); 

  cout << endl << "Complementary cDNA strand:" << endl;
  chain2c.printCDNA(); 

  cout << endl << endl << "Reading frame 3" << endl;
  cout << "Original cDNA strand:" << endl;
  chain3.printCDNA(); 

  cout << endl << "Complementary cDNA strand:" << endl;
  chain3c.printCDNA(); 

} // printCDNAs()


void Protein::printmRNAs() const
{
  cout << "Reading frame 1" << endl;
  cout << "Original mRNA strand:" << endl;
  chain1.printmRNA(); 

  cout << "Complementary mRNA strand:" << endl;
  chain1c.printmRNA();

  cout << endl << endl << "Reading frame 2" << endl;
  cout << "Original mRNA strand:" << endl;
  chain2.printmRNA(); 

  cout << "Complementary mRNA strand:" << endl;
  chain2c.printmRNA();

  cout << endl << endl << "Reading frame 3" << endl;
  cout << "Original mRNA strand:" << endl;
  chain3.printmRNA(); 

  cout << "Complementary mRNA strand:" << endl;
  chain3c.printmRNA();
} // printmRNAs()


void Protein::printSeqs() const
{
  cout << "Reading frame 1" << endl;
  cout << "Original strand:" << endl;
  chain1.printSeq(); 

  cout << "Complementary strand:" << endl;
  chain1c.printSeq();

  cout << endl << endl << "Reading frame 2" << endl;
  cout << "Original strand:" << endl;
  chain2.printSeq(); 

  cout << "Complementary strand:" << endl;
  chain2c.printSeq();

  cout << endl << endl << "Reading frame 3" << endl;
  cout << "Original strand:" << endl;
  chain3.printSeq(); 

  cout << "Complementary strand:" << endl;
  chain3c.printSeq();
} // printSeqs()


void Protein::printStructures() const
{
  cout << "Reading frame 1" << endl;
  cout << "Original amino acid strand:" << endl; 
  chain1.printStructure();

  cout << endl << "Complementary amino acid strand:" << endl;
  chain1c.printStructure();

  cout << endl << endl << "Reading frame 2" << endl;
  cout << "Original amino acid strand:" << endl; 
  chain2.printStructure();

  cout << endl << "Complementary amino acid strand:" << endl;
  chain2c.printStructure();

  cout << endl << endl << "Reading frame 3" << endl;
  cout << "Original amino acid strand:" << endl; 
  chain3.printStructure();

  cout << endl << "Complementary amino acid strand:" << endl;
  chain3c.printStructure();
} // printStructures()


AminoAcidChain& Protein::returnLongestChain()
{
  int num = returnLongestChainNum();

  switch(num)
  {
    case 0: return chain1;
    case 1: return chain1c;
    case 2: return chain2;
    case 3: return chain2c;
    case 4: return chain3;
    case 5: return chain3c;
  } // switch

  return chain1; // to quiet warnings

} // returnLongestChain()


// returns the index that will indicate which AminoAcidChain has the longest
// chain
int Protein::returnLongestChainNum() const
{
  int array[6];

  // assigning the chains with a specific number to distinguish them
  array[0] = chain1.aaLength; array[1] = chain1c.aaLength;
  array[2] = chain2.aaLength; array[3] = chain2c.aaLength;
  array[4] = chain3.aaLength; array[5] = chain3c.aaLength;

  return max_element(array, array + 6) - array;
} // returnLongestChainNum()
