#ifndef AMINOACID_H_
#define AMINOACID_H_

#include <string>
#include "LinearProbing.h"

using namespace std;

class AminoAcid
{
  friend class LinearHashTable;

  char symbol; // one letter symbol
  string abbrev; // 3 letter name
  string full; // full name
  float cpKa; // C-terminal pKa
  float npKa; // N-terminal pKa
  float rpKa; // R-group pKa (if R-group has no pKa, == 0)

  void caseU(const string &cdn, char &s);
  void caseC(const string &cdn, char &s);
  void caseA(const string &cdn, char &s);
  void caseG(const string &cdn, char &s);

public:
  AminoAcid();
  void findAminoAcid(const string &cdn, char &s);
  float getCpKa() const { return cpKa; }
  float getNpKa() const { return npKa; }
  float getRpKa() const { return rpKa; }
  char getSymbol() const { return symbol; }
  bool getTriplet(char symbol, string &t) const;
  void print();
  void updateAA(LinearHashTable AALst, char symb); // update info for AA

}; // class AminoAcid

#endif // AMINOACID_H_