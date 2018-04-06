#include <vector>
#include <iomanip> // for setw
#include "Sequence.h"


Sequence::Sequence(): seqLength(0)
{} // Sequence()


// compares 3 numbers and returns the smallest positive number
// used to determine which stop codon is reached first to end gene
int Sequence::compareThreeNums(int num1, int num2, int num3) const
{
  if(num1 == -1) // first stop codon isn't found
  {
    if(num2 == -1) // second stop codon isn't found, return third stop codon
      return num3;

    else if(num3 == -1) // third stop codon isn't found, return sec stop codon
      return num2;

    else // second and third codons are found, return the smaller position
      return (num2 <= num3) ? num2 : num3;
  } // if num1 is negative

  else if(num2 == -1) // second stop codon isn't found, but first codon is found
  {
    if(num3 == -1) // third codon isn't found, return first codon
      return num1;

    else // first and third stop codons are found, return the smaller position
      return (num1 <= num3) ? num1 : num3;
  } // else if num2 is negative 

  else if(num3 == -1) // third stop codon isn't found, but first and second are
  {
    // first and second stop codons are found, return the smaller position
    return (num1 <= num2) ? num1 : num2;  
  } // if num3 is negative

  else // all stop codons have been found
  {
    int min1or2 = (num1 <= num2) ? num1 : num2;
    return (min1or2 <= num3) ? min1or2 : num3;
  } // else all nums are positive

  return -1; // to quiet warnings
} // compareThreeNums


// searches for a particular codon in a given s (mRNA)
// used to find the start and stop codons in a mRNA sequence
int Sequence::findCodon(const string &s, const string &codon, size_t found, int offset) const
{
  int pos = 0;

  found = s.find(codon, found + offset);

  while(found != string::npos)
  {
    pos = static_cast<int>(found);

    if(pos % 3 == 0 || pos == s.length() - 1)
      return pos;

    found = s.find(codon, found + 1);
  } // while codon is still found in s

  return -1;

} // findCodon()


// returns the index of the first stop codon in the sequence s
int Sequence::findFirstStopCodon(const string &s, size_t findStop, int offset) const
{
  int stop1 = -1, stop2 = -1, stop3 = -1; // initialize to -1 to distinguish
  size_t foundStop1 = findStop; 
  size_t foundStop2 = findStop; 
  size_t foundStop3 = findStop;

  stop1 = findCodon(s, "UAA", foundStop1, offset);
  stop2 = findCodon(s, "UAG", foundStop2, offset);
  stop3 = findCodon(s, "UGA", foundStop3, offset);

  return compareThreeNums(stop1, stop2, stop3);

} // findFirstStopCodon()


// generates a random base for mutation
char Sequence::genRandomBase() const
{
  int num = 0;
  num = rand() % 4;

  switch(num)
  {
    case 0: return 'A'; break;
    case 1: return 'C'; break;
    case 2: return 'G'; break;
    case 3: return 'T'; break;
  } // switch

  return '!'; // to quiet warnings
} // genRandomBase()


// calculates the number of mutations (differences between two amino acid sequences)
// also takes into consideration the change in length of the amino acid sequence
// after mutation
int Sequence::calcNumMutations(const string &orig, const string &mut) const
{
  int numMutations = 0;
  int origLen = orig.length();
  int mutLen = mut.length();

  for(int pos = 0; (pos < origLen) && (pos < mutLen); pos++)
  {
    if(orig[pos] == mut[pos])
      continue;

    else
      numMutations++;
  } // for()

  if(origLen != mutLen)
    numMutations += abs(origLen - mutLen);

  return numMutations;

} // calcNumMutations


// find repeated fragments of at least size 5
void Sequence::findRepeats() const
{
  
  vector<string> repeatFrags; // contains all the repeats
  vector<int> fragCount;      // will keep count of each valid repeat frag
  vector<int> startIndex;     // will keep count of the start location of each valid repeat
 
  cout << "Full sequence (Reading frame 1 original): " << seq << endl << endl;
  findSingleRepeats(repeatFrags, fragCount, startIndex);
  findPairRepeats(repeatFrags, fragCount, startIndex);
  printRepeats(repeatFrags, fragCount, startIndex);

} // findRepeats()


// looks for alternating pairs of bases in a sequence
void Sequence::findPairRepeats(vector<string> &repeats, vector<int> &counts, 
                                 vector<int> &starts) const
{
  int count = 0;
  int prevJ = 0;   // holds j's prev index in the case that end of the sequence has been reached
  string pairFrags;
  string prevPair; // holds prev pairFrag in the case that end of the sequence has been reached
  
  pairFrags.append("ATTAACCAAGGACGGCCTTCGTTG");
  prevPair.append(pairFrags, 0, 2);

  for(int i = 0; i < pairFrags.length() - 1; i += 2)
  {
    for(int j = 0; j < seqLength - 1; j++)
    {
      if(pairFrags[i] == seq[j] && pairFrags[i + 1] == seq[j + 1])
      {
        count++;
        j++; // to jump to next pair
      } // if strings are equal

      else // not equal
      {
        if(count >= 5)
        {
          string c;
          c.append(prevPair, 0, 2);
          repeats.push_back(c);
          counts.push_back(count);
          
          if(j != 0)
            prevJ = j;

          if(pairFrags[i] != prevPair[0] && pairFrags[i+1] != prevPair[1])
            prevJ += 1;

          starts.push_back( prevJ - (count *  2) );
        } // if letter repeats 5 times or more consecutively

        count = 0; // reset
      } // else

      prevJ = j;
      prevPair[0] = pairFrags[i];
      prevPair[1] = pairFrags[i + 1];

    } // for each pair in sequence

  } // for each pair
} // findPairRepeats()


// looks for a single repeating base in sequence
void Sequence::findSingleRepeats(vector<string> &repeats, vector<int> &counts, 
                                 vector<int> &starts) const
{
  int count = 0; // keeps count of current repeatFrag being searched
  string singleFrags;
  singleFrags.append("ACGT");

  for(int i = 0; i < singleFrags.length(); i++)
  {
    for(int j = 0; j < seqLength; j++)
    {
      if(singleFrags[i] == seq[j])
        count++;

      else // not equal
      {
        if(count >= 5)
        {
          string c;
          c.append(1, singleFrags[i]);
          repeats.push_back(c);
          counts.push_back(count);
          starts.push_back(j - count);
        } // if letter repeats 5 times or more consecutively

        count = 0; // reset
      } // else

    } // for each letter in sequence

  } // for each letter
} // findSingleRepeats()


// starts at the 3' end of the sequence and generates complement
// in 5' to 3'
void Sequence::genComplementSeq(const Sequence &s)
{
  for(int i = s.seqLength - 1; i > -1; i--)
  {
    switch(s.seq[i]) 
    {
      case 'A':
        seq.append("T");
        break;
      case 'T':
        seq.append("A");
        break;
      case 'C':
        seq.append("G");
        break;
      case 'G':
        seq.append("C");
        break;
    } //switch()

  } // for each letter in seq
} // genComplementSeq()


// generates the actual mRNA sequence + poly A tail attached
void Sequence::genMRNA(int aLength, int startPos, int stopPos)
{
  if(aLength == 0)
    return;

  int num = rand() % 200;

  for(int i = startPos; i < stopPos + 3; i++)
    mRNA += seqMRNA[i];

  mRNA.append("AAUAAA");
  
  for(int i = 0; i < num; i++) // random number of A's are appended at the end
    mRNA.append("A"); 
} // genMRNA()


// generates a mutated amino acid chain; similar to the non-mutated version of this
// function but only updates the string containing the amino acid sequence to save
// time and space
string Sequence::genMutatedAminoAcidSeq(const LinearHashTable &AALst)
{
  string aaSeq;   // amino acid sequence of to-be-mutated seq
  char symbl = '!';
  AminoAcid temp; // used to just access functions, poor practice though
  int start = 0; int stop = 0;

  searchForLongestChain(seqMRNA, start, stop);

  for(int pos = start; pos < stop; )
  {
    string codon;

    for(int i = 0; i < 3; i++, pos++)
      codon += seqMRNA[pos];
    
    temp.findAminoAcid(codon, symbl);
    aaSeq += symbl;
  } // for each position in cDNA nucleotide sequence

  return aaSeq;
} // genMutatedAminoAcidSeq()


// generates a mutated sequence based on the percentage the user inputs
void Sequence::genMutatedSeq(const LinearHashTable &AALst, const Sequence &s,
                             int percent, int start, int stop)
{
  int randomPos = 0;
  // num of how many bases to mutate
  int numBasesMutated = ( (s.seqLength - 3) * percent) / 100; 
  char mutatedBase = '!';
  bool visited[s.seqLength]; // to prevent each base from being mutated multiple times

  for(int i = 0; i < s.seqLength; i++) 
    visited[i] = false;

  seq.append(s.seq); // copy original seq into seq to be mutated

  for(int count = 0; count < numBasesMutated; count++)
  {
    randomPos = rand() % stop; // generates random number 0 <= x < a.stopPos

    // keep generating a new number if random number is < beginning position 
    // of cDNA (after M)
    while(randomPos < start + 3 || visited[randomPos] == true) 
      randomPos = rand() % stop;  

    mutatedBase = genRandomBase();
    seq[randomPos] = mutatedBase;
    visited[randomPos] = true;
  } // while there are still bases left to mutate

  setSequenceMembers();
} // genMutatedSeq()


// shifts the reading frame for respective sequences' reading frames
void Sequence::genReadingFrameSeqs(Sequence &s2, Sequence &s3) const
{
  for(int i = 1; i < seqLength - 1; i++)
  {
    s2.seq += seq[i];
    s3.seq += seq[i + 1];
  } // for each letter in chain1

  s2.seq += seq[seqLength - 1];
} // genReadingFrameSeqs()


// searches for the longest amino acid sequence that starts with a start codon and
// ends in a stop codon; continues to search in current sequence even if a valid
// gene with a start and stop codon have been found
void Sequence::searchForLongestChain(const string &s, int &start, int &stop) const
{
  int tempLength = 0; int tempLength2 = 0; // length between start and stop codons
  int tempStart = 0; int tempStart2 = 0;   // position of start codons
  int tempStop = 0; int tempStop2 = 0;     // position of stop codons
  size_t foundStart = 0, foundStop = 0;

  tempStart = findCodon(s, "AUG", foundStart, 0);
  tempStop = findFirstStopCodon(s, foundStop, 0);

  if(tempStart == -1 || tempStop == -1 || tempStop < tempStart)
    return; // if no start or stop codon is found, return 

  tempLength = tempStop - tempStart;
  tempStart2 = tempStart;
  tempStop2 = tempStop;

  for(int seqsLeft = s.length() - tempStop; seqsLeft > tempLength; 
      seqsLeft = s.length() - tempStop2 - 3) // - 3 for the stop codon
  {
    // start search for start codon after previous stop codon
    tempStart2 = findCodon(s, "AUG", foundStart, 1 + tempStart2);
    tempStop2 = findFirstStopCodon(s, foundStop, 1 + tempStop2);

    if(tempStart2 == -1 || tempStop2 == -1 || tempStop < tempStart)
      break; // exit loop if another start/stop codon cannot be found

    tempLength2 = tempStop2 - tempStart2;

    if(tempLength > tempLength2)
      continue; // if length of first chain is longer, continue searching

    else // tempLength < tempLength2
    {
      tempLength = tempLength2;
      tempStart = tempStart2;
      tempStop = tempStop2;
    } // else

  } // while it is still possible to generate a sequence longer than initially

  start = tempStart;
  stop = tempStop;

} // searchForLongestChain()


// changes all T's in nucelotide sequence to U's; not the true method of 
// transcription but produces same output
void Sequence::genSeqMRNA()
{
  for(int i = 0; i < seqLength; i++)
  {
    switch(seq[i]) 
    {
      case 'A':
        seqMRNA.append("A");
        break;
      case 'T':
        seqMRNA.append("U");
        break;
      case 'C':
        seqMRNA.append("C");
        break;
      case 'G':
        seqMRNA.append("G");
        break;
    } //switch()

  } // for each nucleotide in seq

} // genMRA()


// reads in the sequence provided in the file and stores it in sequence
istream& operator>> (istream &inf, Sequence &s)
{
  string buffer;

  getline(inf, buffer); // reads in the first line of the file

  while(getline(inf, buffer)) // read in AminoAcidChain until the end of file
    s.seq.append(buffer); 

  return inf;

} // operator>> ()


// prints table of the repeat fragments found
void Sequence::printRepeats(vector<string> repeats, vector<int> counts, 
                            vector<int> starts) const
{
  int size = repeats.size();

  if(size == 0)
  {
    cout << "No repeats of at least length 5 were found in sequence." 
      << endl << endl;
    return;
  } // if

  for(int i = 0; i < size; i++)
  {
    cout << setw(20) << left << "Repeat fragment" << setw(20) << "Num of repeats"
      << setw(20) << "Starting position of repeat" << endl;
    cout << setw(20) << left << repeats[i] << setw(20) << counts[i] << setw(20)
      << starts[i] << endl << endl;
  } // for each item in repeats

} // printRepeats()


// sets the length and generates the seqMRNA
void Sequence::setSequenceMembers()
{
  seqLength = seq.length();
  genSeqMRNA();
} // setMembers()