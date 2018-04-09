#include "AminoAcid.h"
#include "AminoAcidChain.h"
#include "LinearProbing.h"
#include "Protein.h"
#include "Sequence.h"

using namespace std;

void aminoAcidSeq(const Protein &protein);
void calcMutationProb(const LinearHashTable &l, AminoAcidChain* chain);
void findRepeatFrags(const Protein &protein);
void generateSeqs();
int getChoice(int max);
void mRNASeq(const Protein &protein);
void nucleotideSeq(const Protein &protein);
void printStructuralFormMenu();
void printStructuralForms1(const Protein &protein);
void printStructuralForms2(const Protein &protein, float pH);
void returnLongestSequence(const Protein &protein);
void structuralForm(Protein &protein);

static string buffer; // used to clear cin's buffer

int main(int argc, char** argv)
{
  int choice = 0;
  int maxChoice = 8; // maximum number choice user can enter/choose in menu

  Protein protein; // contains all the possible reading frames given a sequence
  LinearHashTable AAList; // contains all amino acid info; used as a reference
  AminoAcidChain* longestChain = NULL; // longest viable gene found
  srand(time(NULL)); // for generating random numbers based on the time

  if(argv[1])
  {
      ifstream inf(argv[1]);
      inf >> protein;
      inf.close();
  } // if file is valid

  AAList.create();
  protein.genAminoAcidSequences(AAList);
  longestChain = &(protein.returnLongestChain());

  // Main Menu
  do 
  {
    cout << endl << "DNA Sequence Manipulation and Analysis" << endl;
    cout << "0. Quit" << endl;
    cout << "1. View longest amino acid sequence obtained" << endl;
    cout << "2. View nucleotide sequences" << endl;
    cout << "3. View mRNA sequences" << endl;
    cout << "4. View amino acid sequences" << endl;
    cout << "5. View structural formula of amino acid sequence" << endl;
    cout << "6. Calculate mutation probability" << endl;
    cout << "7. Generate nucleotide sequences from an amino acid sequence" << endl;
    cout << "8. Find repeat fragments in sequence" << endl << endl;
    choice = getChoice(maxChoice);

    switch(choice)
    {
      case 1: returnLongestSequence(protein); break;
      case 2: nucleotideSeq(protein); break;
      case 3: mRNASeq(protein); break;
      case 4: aminoAcidSeq(protein); break;
      case 5: structuralForm(protein); break;
      case 6: calcMutationProb(AAList, longestChain); break;
      case 7: generateSeqs(); break;
      case 8: findRepeatFrags(protein); break;
    } // switch

  } while (choice > 0);
  
  return 0;
} // main()


// displays the amino acid sequence of the user's choosing
// in one-letter representation
void aminoAcidSeq(const Protein &protein)
{
  int usrChoice = 0;
  int maxChoice = 4;

  do
  {
    cout << "Choose an option:" << endl;
    cout << "0. Go back to Main Menu" << endl;
    cout << "1. View first reading frame amino acid sequence" << endl;
    cout << "2. View second reading frame amino acid sequence" << endl;
    cout << "3. View third reading frame amino acid sequence" << endl;
    cout << "4. View all amino acid sequences" << endl << endl;
    usrChoice = getChoice(maxChoice);

    switch(usrChoice)
    {
      case 1: 
      {
        cout << "Original amino acid strand:" << endl << endl;
        protein.chain1.printAA(); 
        cout << "Complementary amino acid strand:" << endl << endl;
        protein.chain1c.printAA();
        break;
      } // reading frame 1
      case 2: 
      {
        cout << "Original amino acid strand:" << endl << endl;
        protein.chain2.printAA(); 
        cout << "Complementary amino acid strand:" << endl << endl;
        protein.chain2c.printAA();
        break;
      } // reading frame 2
      case 3: 
      {
        cout << "Original amino acid strand:" << endl << endl;
        protein.chain3.printAA(); 
        cout << "Complementary amino acid strand:" << endl << endl;
        protein.chain3c.printAA();
        break;
      } // reading frame 3
      case 4: protein.printAAs(); // all reading frames
              break; 
    } // switch

  } while(usrChoice != 0);

} // aminoAcidSeq()


// calculates the probability that the amino acid sequence will be changed
// based on the percentage of base mutation in the cDNA sequence and
// the number of times the sequence is mutated
void calcMutationProb(const LinearHashTable &l, AminoAcidChain* chain )
{
  int usrPercentage = 0;   // percentage of how much to mutate the sequence
  int usrNumMutations = 0; // how many mutated sequences to generate
  int numMutations = 0;    // number of random sequences that have a mutation
  int mutationProb = 0;

  do
  {
    cout << "Enter mutation percentage: ";
    cin >> usrPercentage;
    cout << endl;

    if(!cin || usrPercentage >= 100 || usrPercentage <= 0)
    {
      cin.clear();
      getline(cin, buffer); // removes rest of user's input
      usrPercentage = 0;
      cout << "Error. Please enter an integer between 0 and 100." 
           << endl << endl;
    } // if usr input isn't valid

  } while(usrPercentage == 0);

  do
  {
    cout << "Enter number of mutations to generate (must be at least 100): ";
    cin >> usrNumMutations;
    cout << endl;

    // chose 100 to be the minimum of random sequences to be generated to get
    // a better probability, but can change if necessary
    if(!cin || usrNumMutations < 100)
    {
      cin.clear();
      getline(cin, buffer); // removes rest of user's input
      usrNumMutations = 0;
      cout << "Error. Please enter an integer greater or equal to 100." 
           << endl << endl;
    } // if usr input isn't valid

  } while(usrNumMutations == 0);

  // generates random sequences
  for(int count = 0; count < usrNumMutations; count++)
  {
    Sequence temp;
    string mutatedAASeq;

    temp.genMutatedSeq(l, chain->sequence, usrPercentage, chain->metPos, chain->stopPos);
    mutatedAASeq = temp.genMutatedAminoAcidSeq(l);
    numMutations += temp.calcNumMutations(chain->aaSeq, mutatedAASeq);
  } // for each mutation count

  mutationProb = ( numMutations * 100 / (usrNumMutations * chain->aaLength) );
  cout << "Probability of amino acid sequence mutations at " 
       << usrPercentage << "% mutated bases out of " << usrNumMutations 
       << " mutated sequences is " << mutationProb << "%" << endl << endl;

} // calcMutationProb()


// finds fragments in the main sequence that are repeated
void findRepeatFrags(const Protein &protein)
{
  protein.chain1.sequence.findRepeats();
} // findRepeatFrags()


// user enters in an amino acid sequence and a DNA sequence will be generated 
// from it
void generateSeqs()
{
  AminoAcid temp; // again, poor practice bc only using it to access class's fcn
  string aa, seq, triplet;
  int num = rand() % 3; // generate a random number s.t. 0 <= rand < 3

  cout << "Enter a valid amino acid sequence with no spaces or any other extra characters. ";
  cout << "Please use the one-letter symbol representation (in uppercase) for each amino acid." << endl;
  cout << endl << "Example: \"MRTGE\" will produce the following amino acid sequence:" << endl;
  cout << "Methionine - Arginine - Threonine - Glycine - Glutamic Acid" << endl << endl;
  cout << "Your sequence >> ";
  cin >> aa;
  cout << endl;

  for(int i = 0; i < aa.length(); i++)
  {
    if(temp.getTriplet(aa[i], triplet))
      seq.append(triplet);

    else // user entered in invalid sequence
    {
      cout << "Error: At least one invalid character was detected. Please try again." << endl;
      return;
    } // else sequence is invalid

  } // for each amino acid in user's sequence

  switch(num)
  {
    case 0: seq.append("ATT"); break;
    case 1: seq.append("ATC"); break;
    case 2: seq.append("ACT"); break;
  } // to append 1 out of the 3 possible stop codons at the end

  cout << "Generated DNA sequence is: \n5' " << seq << " 3'" << endl << endl;

} // generateSeqs()


// reads input from user and validates it; max is the max number user can enter
// as a valid choice
int getChoice(int max)
{
  int usrChoice;

  do
  {
    cout << "Your Choice >> ";
    cin >> usrChoice;

    if(!cin)
    {
      cout << "\nPlease enter a valid number from 0 to " << max << ".";
      cin.clear();
      getline(cin, buffer); // removes rest of user's input
      usrChoice = -1;       // to prompt menu again
    } // if cin fails

    cout << endl;
  } while(usrChoice < 0 || usrChoice > max);

  return usrChoice;
} // getChoice()


// prints the mRNA sequence of user's choice
void mRNASeq(const Protein &protein)
{
  int usrChoice = 0;
  int maxChoice = 4;

  do
  {
    cout << "Choose an option:" << endl;
    cout << "0. Go back to Main Menu" << endl;
    cout << "1. View first reading frame mRNA sequence" << endl;
    cout << "2. View second reading frame mRNA sequence" << endl;
    cout << "3. View third reading frame mRNA sequence" << endl;
    cout << "4. View all mRNA sequences" << endl << endl;
    usrChoice = getChoice(maxChoice);

    switch(usrChoice)
    {
      case 1: 
      {
        cout << "Original mRNA strand:" << endl;
        protein.chain1.printmRNA(); 
        cout << "Complementary mRNA strand:" << endl;
        protein.chain1c.printmRNA();
        break;
      } // case 1
      case 2: 
      {
        cout << "Original mRNA strand:" << endl;
        protein.chain2.printmRNA(); 
        cout << "Complementary mRNA strand:" << endl;
        protein.chain2c.printmRNA();
        break;
      } // case 2
      case 3: 
      {
        cout << "Original mRNA strand:" << endl;
        protein.chain3.printmRNA(); 
        cout << "Complementary mRNA strand:" << endl;
        protein.chain3c.printmRNA();
        break;
      } // case 3
      case 4: protein.printmRNAs(); break;
    } // switch

  } while(usrChoice > 0);

} // mRNASeq()


// prints nucleotide sequence of the user's choice
void nucleotideSeq(const Protein &protein)
{
  int usrChoice = 0;
  int maxChoice = 8;

  do
  {
    cout << "Choose an option:" << endl;
    cout << "0. Go back to Main Menu" << endl;
    cout << "1. View first reading frame FULL nucleotide sequence" << endl;
    cout << "2. View second reading frame FULL nucleotide sequence" << endl;
    cout << "3. View third reading frame FULL nucleotide sequence" << endl;
    cout << "4. View all FULL nucleotide sequences" << endl;
    cout << "5. View first reading frame cDNA nucleotide sequence" << endl;
    cout << "6. View second reading frame cDNA nucleotide sequence" << endl;
    cout << "7. View third reading frame cDNA nucleotide sequence" << endl;
    cout << "8. View all cDNA nucleotide sequences" << endl << endl;
    usrChoice = getChoice(maxChoice);

    switch(usrChoice)
    {
      case 1: 
      {
        cout << "Original strand:" << endl;
        protein.chain1.printSeq(); 
        cout << "Complementary strand:" << endl;
        protein.chain1c.printSeq();
        break;
      } // case 1
      case 2: 
      {
        cout << "Original strand:" << endl;
        protein.chain2.printSeq(); 
        cout << "Complementary strand:" << endl;
        protein.chain2c.printSeq();
        break;
      } // case 2
      case 3: 
      {
        cout << "Original strand:" << endl;
        protein.chain3.printSeq(); 
        cout << "Complementary strand:" << endl;
        protein.chain3c.printSeq();
        break;
      } // case 3

      case 4: protein.printSeqs(); break;
      case 5:
      {
        cout << "Original cDNA strand:" << endl;
        protein.chain1.printCDNA(); 
        cout << "Complementary cDNA strand:" << endl;
        protein.chain1c.printCDNA();
        break;
      } // case 5
      case 6:
      {
        cout << "Original cDNA strand:" << endl;
        protein.chain2.printCDNA(); 
        cout << "Complementary cDNA strand:" << endl;
        protein.chain2c.printCDNA();
        break;
      } // case 6
      case 7:
      {
        cout << "Original cDNA strand:" << endl;
        protein.chain3.printCDNA(); 
        cout << "Complementary cDNA strand:" << endl;
        protein.chain3c.printCDNA();
        break;
      } // case 7
      case 8: protein.printCDNAs(); break;
    } // switch

  } while(usrChoice > 0);

} // nucleotideSeq()


// prints the menu for when the user wants to print the structural form
// of the amino acid sequence of the user's choice
void printStructuralFormMenu()
{
  cout << "Choose an option:" << endl;
  cout << "0. Go back to pH selection" << endl;
  cout << "1. View original first reading frame amino acid structure" << endl;
  cout << "2. View complementary first reading frame amino acid structure" << endl;
  cout << "3. View original second reading frame amino acid structure" << endl;
  cout << "4. View complementary second reading frame amino acid structure" << endl;
  cout << "5. View original third reading frame amino acid structure" << endl;
  cout << "6. View complementary third reading frame amino acid structure" << endl;
  cout << "7. View all amino acid structures" << endl << endl;
} // printStructuralFormMenu


// prints the structural form of the amino acid chain @ pH 7
// all AminoAcidChains' structures are constructed @ pH 7
void printStructuralForms1(const Protein &protein)
{
  int usrChoice = -1;
  int maxChoice = 8;

  do
  {
    printStructuralFormMenu();
    usrChoice = getChoice(maxChoice);

    switch(usrChoice)
    {
      case 1: 
      {
        cout << "Original amino acid strand:" << endl; 
        protein.chain1.printStructure(); 
        break;
      } // case 1
      case 2: 
      {
        cout << "Complementary amino acid strand:" << endl; 
        protein.chain1c.printStructure(); 
        break;
      } // case 2
      case 3: 
      {
        cout << "Original amino acid strand:" << endl; 
        protein.chain2.printStructure(); 
        break;
      } // case 3
      case 4: 
      {
        cout << "Complementary amino acid strand:" << endl; 
        protein.chain2c.printStructure(); 
        break;
      } // case 4
      case 5:
      {
        cout << "Original amino acid strand:" << endl; 
        protein.chain3.printStructure(); 
        break;
      } // case 5
      case 6:
      {
        cout << "Complementary amino acid strand:" << endl; 
        protein.chain3c.printStructure(); 
        break;
      } // case 6
      case 7: protein.printStructures(); break;
    } // switch

  } while(usrChoice > 0);

} // printStructuralForms1()


// creates a copy of the AminoAcidChain of user's choice and alters it based on
// the pH the user enters
void printStructuralForms2(const Protein &protein, float pH)
{
  int usrChoice = -1;
  int maxChoice = 7; // from printStructuralFormMenu

  Structure structure1(238, 35); Structure structure1c(238, 35); // chain1 copy
  Structure structure2(238, 35); Structure structure2c(238, 35); // chain2 copy
  Structure structure3(238, 35); Structure structure3c(238, 35); // chain3 copy
  structure1.createAAStructure(protein.chain1.AASeq, pH, protein.chain1.getaaLength());
  structure1c.createAAStructure(protein.chain1c.AASeq, pH, protein.chain1c.getaaLength());
  structure2.createAAStructure(protein.chain2.AASeq, pH, protein.chain2.getaaLength());
  structure2c.createAAStructure(protein.chain2c.AASeq, pH, protein.chain2c.getaaLength());
  structure3.createAAStructure(protein.chain3.AASeq, pH, protein.chain3.getaaLength()); 
  structure3c.createAAStructure(protein.chain3c.AASeq, pH, protein.chain3c.getaaLength());

  do 
  {
    printStructuralFormMenu();
    usrChoice = getChoice(maxChoice);

    switch(usrChoice)
    {
      case 1: structure1.print(); break;
      case 2: structure1c.print(); break;
      case 3: structure2.print(); break;
      case 4: structure2c.print(); break;
      case 5: structure3.print(); break;
      case 6: structure3c.print(); break;
      case 7: 
      {
        cout << "Reading frame 1" << endl;
        structure1.print(); cout << endl; structure1c.print();
        cout << endl << endl << "Reading frame 2" << endl;
        structure2.print(); cout << endl; structure2c.print();
        cout << endl << endl << "Reading frame 3" << endl;
        structure3.print(); cout << endl; structure3c.print(); 
        break;
      } // case 7

    } // switch()

  } while(usrChoice > 0);
  
} // printStructuralForms2()


// outputs information of the longest possible gene of sequence 
void returnLongestSequence(const Protein &protein)
{
  int num = protein.returnLongestChainNum();

  cout << "Longest amino acid sequence obtained is: ";

  if(num % 2 == 0)
  {
    cout << "original strand of ";
    if(num == 0)
    {
      cout << "reading frame 1" << endl;
      cout << "Nucleotide sequence of longest gene:" << endl;
      protein.chain1.printSeq();
      cout << endl;
      cout << "Complementary sequence of longest gene:" << endl;
      protein.chain1c.printSeq();
      protein.chain1.printAll();
    }

    else if(num == 2)
    {
      cout << "reading frame 2" << endl;
      cout << "Nucleotide sequence of longest gene:" << endl;
      protein.chain2.printSeq();
      cout << endl;
      cout << "Complementary sequence of longest gene:" << endl;
      protein.chain2c.printSeq();
      protein.chain2.printAll();
    }

    else // num == 4
    {
      cout << "reading frame 3" << endl;
      cout << "Nucleotide sequence of longest gene:" << endl;
      protein.chain3.printSeq();
      cout << endl;
      cout << "Complementary sequence of longest gene:" << endl;
      protein.chain3c.printSeq();
      protein.chain3.printAll();
    }

  } // if even

  else // num % 2 == 1
  {
    cout << "complementary strand of ";
    if(num == 1)
    {
      cout << "reading frame 1" << endl;
      cout << "Nucleotide sequence of longest gene:" << endl;
      protein.chain1c.printSeq();
      cout << endl;
      cout << "Complementary sequence of longest gene:" << endl;
      protein.chain1.printSeq();
      protein.chain1c.printAll();
    }

    else if(num == 3)
    {
      cout << "reading frame 2" << endl;
      cout << "Nucleotide sequence:" << endl;
      protein.chain2c.printSeq();
      cout << endl;
      cout << "Complementary sequence of longest gene:" << endl;
      protein.chain2.printSeq();
      protein.chain2c.printAll();
    }

    else // num == 5
    {
      cout << "reading frame 3" << endl;
      cout << "Nucleotide sequence of longest gene:" << endl;
      protein.chain3c.printSeq();
      cout << endl;
      cout << "Complementary sequence of longest gene:" << endl;
      protein.chain3.printSeq();
      protein.chain3c.printAll();
    }

  } // else odd

} // returnLongestSequence


// asks user for which pH to print structure at
void structuralForm(Protein &protein)
{
  int usrChoice = - 1; float pH = - 1;
  int maxChoice = 2;

  do
  {
    cout << "Choose an option:" << endl;
    cout << "0. Go back to Main Menu" << endl;
    cout << "1. View structure(s) at pH 7.0" << endl;
    cout << "2. View structure(s) at a specific pH" << endl;
    usrChoice = getChoice(maxChoice);

    switch(usrChoice)
    {
      case 1: printStructuralForms1(protein); break;
      case 2: 
      {
        do
        {
          cout << "Enter a pH value: "; 
          cin >> pH;
          cout << endl;

          if(cin.fail())
          {
            cin.clear();
            cin.ignore();
            cout << "Error. Please enter a valid decimal number." << endl;
            cout << "Enter a pH value: "; 
            cin >> pH;
            cout << endl;
          } // if user inputs something other than a number
        } while(cin.fail());

        printStructuralForms2(protein, pH);
        break;
      } // case 2

    } // switch

  } while(usrChoice > 0);

} // structuralForm()

