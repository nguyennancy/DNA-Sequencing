#include "AminoAcid.h"

AminoAcid::AminoAcid(): symbol('!'), cpKa(-1), npKa(-1), rpKa(-1)
{}


void AminoAcid::findAminoAcid(const string &cdn, char &s)
{
  switch(cdn[0])
  {
    case 'U':
      caseU(cdn, s);
      break;

    case 'C':
      caseC(cdn, s);
      break;

    case 'A':
      caseA(cdn, s);
      break;

    case 'G':
      caseG(cdn, s);
      break;
  } // switch

} // findAminoAcid()


// first letter in codon is 'A'
void AminoAcid::caseA(const string &cdn, char &s)
{
  switch(cdn[1])
  {
    case 'U':
    {
      if(cdn[2] == 'G')
        s = 'M'; // Met / M

      else // cdn[2] == 'U' || 'C' || 'A'
        s = 'I'; // Ile / I

      break;
    } // case U

    case 'C':
      s = 'T'; // Thr / T
      break;

    case 'A':
    {
      if(cdn[2] == 'U' || cdn[2] == 'C')
        s = 'N'; // Asn / N

      else // cdn[2] == 'A' || 'G'
        s = 'K'; // Lys / K

      break;
    } // case A

    case 'G':
    {
      if(cdn[2] == 'U' || cdn[2] == 'C')
        s = 'S'; // Ser / S

      else // cdn[2] == 'A' || cdn[2] == 'G'
        s = 'R'; // Arg / R

      break;
    } // case G

  } // switch
} // caseA()


// first letter in codon is 'C'
void AminoAcid::caseC(const string &cdn, char &s)
{
  switch(cdn[1])
  {
    case 'U':
      s = 'L'; // Leu / L
      break;

    case 'C':
      s = 'P'; // Pro / P
      break;

    case 'A':
    {
      if(cdn[2] == 'U' || cdn[2] == 'C')
        s = 'H'; // His / H

      else // cdn[2] == 'A' || cdn[2] == 'G'
        s = 'Q'; // Gln / Q

      break;
    } // case A

    case 'G':
      s = 'R'; // Arg / R
      break;

  } // switch
} // caseC()


// first letter in codon is 'G'
void AminoAcid::caseG(const string &cdn, char &s)
{
  switch(cdn[1])
  {
    case 'U':
      s = 'V'; // Val / V
      break;

    case 'C':
      s = 'A'; // Ala / A
      break;

    case 'A':
    {
      if(cdn[2] == 'U' || cdn[2] == 'C')
        s = 'D'; // Asp / D

      else // cdn[2] == 'A' || cdn[2] == 'G'
        s = 'E'; // Glu / E

      break;
    } // case A

    case 'G':
      s = 'G'; // Gly / G
      break;
      
  } // switch
} // caseG()


// first letter in codon is 'U'
void AminoAcid::caseU(const string &cdn, char &s)
{
  switch(cdn[1])
  {
    case 'U':
    {
      if(cdn[2] == 'U' || cdn[2] == 'C')
        s = 'F'; // Phe / F

      else // cdn[2] == 'A' || cdn[2] == 'G'
        s = 'L'; // Leu / L

      break;
    } // case U

    case 'C':
      s = 'S'; // Ser / S
      break;

    case 'A':
    {
      if(cdn[2] == 'U' || cdn[2] == 'C')
        s = 'Y'; // Tyr / Y

      else // cdn[2] == 'A' || 'G'
        s = 'X'; // STOP codon

      break;
    } // case A

    case 'G':
    {
      if(cdn[2] == 'A')
        s = 'X'; // STOP codon

      else if(cdn[2] == 'G')
        s = 'W'; // Trp / W

      else // cdn[2] == 'U' || 'C'
        s = 'C'; // Cys /C

      break;
    } // case G

  } // switch 
} // caseU()


// for generating a random DNA sequence based on user's amino acid sequence
// returns false if user enters in an invalid char; returns true otherwise
bool AminoAcid::getTriplet(char symbol, string &t) const
{
  int num = 0;
  switch(symbol)
  {
    case 'A':
    {
      num = rand() % 4;
      switch(num)
      {
        case 0: t = "GCT"; return true;
        case 1: t = "GCC"; return true;
        case 2: t = "GCA"; return true;
        case 3: t = "GCG"; return true;
      }
    } // case A
    case 'R':
    {
      num = rand() % 6;
      switch(num)
      {
        case 0: t = "AGA"; return true;
        case 1: t = "AGG"; return true;
        case 2: t = "CGT"; return true;
        case 3: t = "CGC"; return true;
        case 4: t = "CGA"; return true;
        case 5: t = "CGG"; return true;
      }
    }  // case R
    case 'N':
    {
      num = rand() % 2;
      switch(num)
      {
        case 0: t = "AAT"; return true;
        case 1: t = "AAC"; return true;
      }
    } // case N
    case 'D': 
    {
      num = rand() % 2;
      switch(num)
      {
        case 0: t = "GAT"; return true;
        case 1: t = "GAC"; return true;
      }
    } // case D
    case 'C':
    {
      num = rand() % 2;
      switch(num)
      {
        case 0: t = "TGT"; return true;
        case 1: t = "TGC"; return true;
      }
    } // case C
    case 'E':
    {
      num = rand() % 2;
      switch(num)
      {
        case 0: t = "GAA"; return true;
        case 1: t = "GAG"; return true;
      }
    }  // case E
    case 'Q':
    {
      num = rand() % 2;
      switch(num)
      {
        case 0: t = "CAA"; return true;
        case 1: t = "CAG"; return true;
      }
    } // case Q
    case 'G': 
    {
      num = rand() % 4;
      switch(num)
      {
        case 0: t = "GGT"; return true;
        case 1: t = "GGC"; return true;
        case 2: t = "GGA"; return true;
        case 3: t = "GGG"; return true;
      }
    } // case G
    case 'H': 
    {
      num = rand() % 2;
      switch(num)
      {
        case 0: t = "CAT"; return true;
        case 1: t = "CAC"; return true;
      }
    } // case H
    case 'I': 
    {
      num = rand() % 3;
      switch(num)
      {
        case 0: t = "ATT"; return true;
        case 1: t = "ATC"; return true;
        case 2: t = "ATA"; return true;
      }
    } // case I
    case 'L':
    {
      num = rand() % 6;
      switch(num)
      {
        case 0: t = "TTA"; return true;
        case 1: t = "TTG"; return true;
        case 2: t = "CTT"; return true;
        case 3: t = "CTC"; return true;
        case 4: t = "CTA"; return true;
        case 5: t = "CTG"; return true;
      }
    } // case L
    case 'K': 
    {
      num = rand() % 2;
      switch(num)
      {
        case 0: t = "AAA"; return true;
        case 1: t = "AAG"; return true;
      }
    } // case K
    case 'M': t = "ATG"; return true;
    case 'F':
    {
      num = rand() % 2;
      switch(num)
      {
        case 0: t = "TTT"; return true;
        case 1: t = "TTC"; return true;
      }
    } // case F 
    case 'P': 
    {
      num = rand() % 4;
      switch(num)
      {
        case 0: t = "CCT"; return true;
        case 1: t = "CCC"; return true;
        case 2: t = "CCA"; return true;
        case 3: t = "CCG"; return true;
      }
    } // case P
    case 'S': 
    {
      num = rand() % 2;
      switch(num)
      {
        case 0: t = "AGT"; return true;
        case 1: t = "AGC"; return true;
      }
    } // case S
    case 'T':
    {
      num = rand() % 4;
      switch(num)
      {
        case 0: t = "ACT"; return true;
        case 1: t = "ACC"; return true;
        case 2: t = "ACA"; return true;
        case 3: t = "ACG"; return true;
      }
    }  // case T
    case 'W': t = "TGG"; return true;
    case 'Y':
    {
      num = rand() % 2;
      switch(num)
      {
        case 0: t = "TAT"; return true;
        case 1: t = "TAC"; return true;
      }
    }  // case Y
    case 'V': 
    {
      num = rand() % 4;
      switch(num)
      {
        case 0: t = "GTT"; return true;
        case 1: t = "GTC"; return true;
        case 2: t = "GTA"; return true;
        case 3: t = "GTG"; return true;
      }
    } // case V
  } // switch

  return false;
} // getTriplet()


void AminoAcid::print()
{
  cout << symbol << " " << abbrev << " " << full << endl;
} // print()


// copies all the necessary information stored AAInfo into Amino Acid equivalent
void AminoAcid::updateAA(LinearHashTable AALst, char symb)
{
  symbol = symb;
  const AAInfo AAinf = AALst.find(symb);
  abbrev.append(AAinf.threeLetter);
  full.append(AAinf.fullName);
  cpKa = AAinf.cpKa;
  npKa = AAinf.npKa;
  rpKa = AAinf.rpKa;
} // AminoAcid()
