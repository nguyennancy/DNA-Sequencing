#ifndef LINEAR_PROBING_H
#define LINEAR_PROBING_H

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

class AAInfo
{
public:
  char letter; // one letter symbol
  int ID; // for hashing / inserting
  string threeLetter;
  string fullName;
  float cpKa; // C-terminal pKa
  float npKa; // T-terminal pKa
  float rpKa; // R-group pKa (if R-group has no pKa, == 0)

  AAInfo(char let, string threeLet, string fullNme, float cP, float tP, 
         float rP);
  AAInfo& operator=(const AAInfo &rhs);
}; // class AAInfo


class LinearHashTable
{
  int count;
  AAInfo **array;
  int size;

  void insert( const AAInfo &x );

public:
  LinearHashTable();
  void create();
  const AAInfo& find(char symbol);
  void print() const;
  
}; // class LinearHashTable

#endif // LINEAR_PROBING_H