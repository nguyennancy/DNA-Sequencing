#ifndef STRUCTURE_H_
#define STRUCTURE_H_

#include <iostream>
#include "AminoAcid.h"

using namespace std;

class Structure
{
  char *box; // implementation of a 2-D array to draw AA chain
  int width;
  int height;
  int size; // total size of box
  int midRow; // row number of where N-C-C backbone is located
  int curCol; // keeps track of where to insert next AA

  int index(int x, int y) const { return x + width * y; }
  void createBackbone1(); // "R-group above" backbone
  void createBackbone2(); // "R-group below" backbone
  void createA(int position);                       // Alanine
  void createR(int position, float pH, float rpKa); // Arginine
  void createN(int position);                       // Asparagine
  void createD(int position, float pH, float rpKa); // Aspartic Acid
  void createC(int position, float pH, float rpKa); // Cysteine
  void createE(int position, float pH, float rpKa); // Glutamic Acid
  void createQ(int position);                       // Glutamine
  void createG(int position);                       // Glycine
  void createH(int position, float pH, float rpKa); // Histidine
  void createI(int position);                       // Isoleucine
  void createL(int position);                       // Leucine
  void createK(int position, float pH, float rpKa); // Lysine
  void createM(int position);                       // Methionine
  void createF(int position);                       // Phenylalanine
  void createP(int position);                       // Proline
  void createS(int position);                       // Serine
  void createT(int position);                       // Tyrosine
  void createW(int position);                       // Tryptophan
  void createY(int position, float pH, float rpKa); // Tyrosine
  void createV(int position);                       // Valine
  void createCH2Above(int row); // creates a CH2 group above backbone
  void createCH2Below(int row); // creates a CH2 group below backbone
  void createNH2Above(int row); // creates a NH2 group above backbone
  void createNH2Below(int row); // creates a NH2 group below backbone

public:
  Structure(int w, int h);
  void createAAStructure(AminoAcid **chain, float pH, int length);
  void print() const;
  void testPrint();

}; // class Structure

#endif // STRUCTURE_H_