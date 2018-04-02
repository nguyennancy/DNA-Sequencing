// Sources: https://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
//          January 31st, 2018 - used for its.eof() in getline()

#include "AminoAcidChain.h"


AminoAcidChain::AminoAcidChain(): metPos(0), stopPos(0), AASeq(NULL), aaLength(0)
{} // AminoAcidChain()


// generates the amino acid sequence/chain based on the mRNA sequence
void AminoAcidChain::genAminoAcidSeq(const LinearHashTable &AALst)
{
  char symbl = '!';
  int start = 0; // position of the start codon
  int stop = 0; // position of the stop codon
  string seqMRNA = sequence.getSeqMRNA(); // copy the sequence

  sequence.searchForLongestChain(seqMRNA, start, stop);
  metPos = start;
  stopPos = stop; 

  for(int pos = start; pos < stop;)
  {
    string codon;

    for(int i = 0; i < 3; i++, pos++)
      codon += seqMRNA[pos];

    AASeq[aaLength] = new AminoAcid();
    AASeq[aaLength]->findAminoAcid(codon, symbl);
    AASeq[aaLength]->updateAA(AALst, symbl);
    aaSeq += symbl;
    aaLength++;

  } // for each nucleotide

  // generate mRNA after knowing start and stop pos, and aaLength
  sequence.genMRNA(aaLength, start, stop); 

} // genAminoAcidSeq2()


// generates complementary strand of AminoAcidChain a
void AminoAcidChain::genComplement(const AminoAcidChain &a)
{
  sequence.genComplementSeq(a.sequence);
} // genComplement()


// shifts the reading frame to the sequence's respective frame
void AminoAcidChain::genReadingFrames(AminoAcidChain &c2, AminoAcidChain &c3) const
{
  sequence.genReadingFrameSeqs(c2.sequence, c3.sequence);
  c2.setAAChainMembers();
  c3.setAAChainMembers();
  
} // genReadingFrames()


// creates the structural form of amino acid @ default pH 7
void AminoAcidChain::genStructures()
{
  structure = new Structure(238, 35); // 238 is approx width of the full screen
  structure->createAAStructure(AASeq, 7.0, aaLength);
} // genStructures()


// read in file
istream& operator>> (istream &inf, AminoAcidChain &a)
{
  inf >> a.sequence;

  a.setAAChainMembers();

  return inf;
} // operator>> ()


void AminoAcidChain::printAA() const
{
  for(int i = 0; i < aaLength - 1; i++)
    cout << aaSeq[i] << " - ";
  cout << aaSeq[aaLength - 1] << endl << endl;
  cout << "Amino Acid Sequence length: " << aaLength << endl << endl << endl;
} // printAA()


// used for printing information about the longest possible gene from said sequence
void AminoAcidChain::printAll() const
{
  cout << "mRNA sequence:" << endl;
  printmRNA();
  cout << endl;

  cout << "Amino Acid sequence:" << endl;
  printAA();
  cout << endl;

  cout << "Amino Acid sequence in structural form:" << endl;
  printStructure();
} // printAll()


void AminoAcidChain::printCDNA() const
{
  cout << "5' "<< sequence.getCDNA(metPos, stopPos) << " 3'" << endl << endl;
} // print()


void AminoAcidChain::printmRNA() const
{
  cout << "5' "<< sequence.getMRNA() << " 3'"<< endl << endl;
} // printmRNA()


void AminoAcidChain::printSeq() const
{
  cout << "5' " << sequence.getSeq() << " 3'" << endl << endl;
} // printSeq()


void AminoAcidChain::printStructure() const
{
  structure->print();
  cout << "Amino Acid Sequence length: " << aaLength << endl << endl;
} // printOrigStructure()


// allocates space for the amino acid sequence and sets corresponding members
// note: more space is allocated than may be needed but always enough
void AminoAcidChain::setAAChainMembers()
{
  sequence.setSequenceMembers();
  AASeq = new AminoAcid*[sequence.getSeqLength() / 3];
} // setSequenceMembers()
