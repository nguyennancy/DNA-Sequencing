#include "Structure.h"

// NOTES: 
// 1. midRow is approx. the middle row of the grid;
// 2. order of atoms and bonds are by top-most row to bottom-most row, then
//    left-most column to right-most column
// 3. odd R group linear formula grouping is to help visualize it in grid


Structure::Structure(int w, int h): width(w), height(h), midRow(16), curCol(0)
{
  size = width * height;
  box = new char[size];
  for(int i = 0; i < size; i++)
    box[i] = ' ';
} // Structure()


void Structure::createBackbone1() // above N-C-C backbone
{
  box[ index(curCol, midRow) ] = 'N';         box[ index(curCol + 2, midRow) ] = '-'; 
  box[ index(curCol + 4, midRow) ] = 'C';     box[ index(curCol + 6, midRow) ] = '-'; 
  box[ index(curCol + 8, midRow) ] = 'C';     box[ index(curCol + 10, midRow) ] = '-';
  box[ index(curCol, midRow + 1) ] = '|';     box[ index(curCol + 4, midRow + 1) ] = '|';
  box[ index(curCol + 8, midRow + 1) ] = '|'; box[ index(curCol + 9, midRow + 1) ] = '|';
  box[ index(curCol, midRow + 2) ] = 'H';     box[ index(curCol + 4, midRow + 2) ] = 'H';
  box[ index(curCol + 8, midRow + 2) ] = 'O'; box[ index(curCol + 4, midRow - 1) ] = '|';
  box[ index(curCol + 4, midRow - 2) ] = '|'; box[ index(curCol + 4, midRow - 3) ] = '|';
} // createBackbone1()


void Structure:: createBackbone2() // below N-C-C backbone
{
  box[ index(curCol, midRow - 2) ] = 'H';     box[ index(curCol + 4, midRow - 2) ] = 'H';
  box[ index(curCol + 8, midRow - 2) ] = 'O'; box[ index(curCol, midRow - 1) ] = '|';
  box[ index(curCol + 4, midRow - 1) ] = '|'; box[ index(curCol + 8, midRow - 1) ] = '|';
  box[ index(curCol + 9, midRow - 1) ] = '|';
  box[ index(curCol, midRow) ] = 'N';         box[ index(curCol + 2, midRow) ] = '-';
  box[ index(curCol + 4, midRow) ] = 'C';     box[ index(curCol + 6, midRow) ] = '-';
  box[ index(curCol + 8, midRow) ] = 'C';     box[ index(curCol + 10, midRow) ] = '-';
  box[ index(curCol + 4, midRow + 1) ] = '|'; box[ index(curCol + 4, midRow + 2) ] = '|';
  box[ index(curCol + 4, midRow + 3) ] = '|'; 
} // createBackbone2()


void Structure::createA(int position)
{
  if(position % 2 == 0)
  {
    // R group: (CH3)
    box[ index(curCol + 4, midRow - 6) ] = 'H'; box[ index(curCol + 4, midRow - 5) ] = '|';
    box[ index(curCol, midRow - 4) ] = 'H';     box[ index(curCol + 2, midRow - 4) ] = '-';
    box[ index(curCol + 4, midRow - 4) ] = 'C'; box[ index(curCol + 6, midRow - 4) ] = '-'; 
    box[ index(curCol + 8, midRow - 4) ] = 'H';
    
    // backbone
    createBackbone1();
  } // if position is even --> R group is above backbone

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group (CH3)
    box[ index(curCol, midRow + 4) ] = 'H';     box[ index(curCol + 2, midRow + 4) ] = '-';
    box[ index(curCol + 4, midRow + 4) ] = 'C'; box[ index(curCol + 6, midRow + 4) ] = '-';
    box[ index(curCol + 8, midRow + 4) ] = 'H'; box[ index(curCol + 4, midRow + 5) ] = '|'; 
    box[ index(curCol + 4, midRow + 6) ] = 'H';
  } // else position is odd --> R group is below backbone

  curCol += 12; // shift over current col 
} // createA()


void Structure::createR(int position, float pH, float rpKa)
{
  if(position % 2 == 0)
  {
    // R group: (CH2)(CH2)(CH2)(NH)(C(NH2)(=NH))
    createNH2Above(midRow - 16);
    if((pH - rpKa) < 0) 
    {
      box[ index(curCol + 9, midRow - 13) ] = '+'; 
      box[ index(curCol + 8, midRow - 11) ] = '|';
      box[ index(curCol + 8, midRow - 10) ] = 'H';
    } // if pH < pKa, protonate R group
     
    box[ index(curCol + 4, midRow - 13) ] = '|';  box[ index(curCol + 4, midRow - 12) ] = 'C';
    box[ index(curCol + 6, midRow - 12) ] = '=';  box[ index(curCol + 8, midRow - 12) ] = 'N';
    box[ index(curCol + 10, midRow - 12) ] = '-'; box[ index(curCol + 12, midRow - 12) ] = 'H';
    box[ index(curCol + 4, midRow - 11) ] = '|';  box[ index(curCol, midRow - 10) ] = 'H';  
    box[ index(curCol + 2, midRow - 10) ] = '-';  box[ index(curCol + 4, midRow - 10) ] = 'N'; 
    
    // creates 3 CH2 groups
    for(int i = 0; i < 3; i++)
      createCH2Above(midRow - 8 + (i*2));

    // backbone
    createBackbone1();
  } // if even, R group above backbone

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group: (CH2)(CH2)(CH2)(NH)(C(NH2)(=NH))
    // creates 3 CH2 groups
    for(int i = 0; i < 3; i++)
      createCH2Below(midRow + 4 + (i*2));

    box[ index(curCol + 4, midRow + 10) ] = 'N'; box[ index(curCol + 6, midRow + 10) ] = '-';
    box[ index(curCol + 8, midRow + 10) ] = 'H'; box[ index(curCol + 4, midRow + 11) ] = '|'; 
    box[ index(curCol - 4, midRow + 12) ] = 'H'; box[ index(curCol - 2, midRow + 12) ] = '-';
    box[ index(curCol, midRow + 12) ] = 'N';     box[ index(curCol + 2, midRow + 12) ] = '=';
    box[ index(curCol + 4, midRow + 12) ] = 'C'; box[ index(curCol + 4, midRow +13) ] = '|';

    if((pH - rpKa) < 0) 
    {
      box[ index(curCol - 1, midRow + 11) ] = '+';
      box[ index(curCol, midRow + 13) ] = '|'; 
      box[ index(curCol, midRow + 14) ] = 'H';
    } // if pH < pKa, protonate R group
      
    createNH2Below(midRow + 14);
  } // else odd, R group below backbone

  curCol += 12; // shift over current col
} // createR()


void Structure::createN(int position)
{
  if(position % 2 == 0)
  {
    // R group: (CH2)(C(=O)(NH2))
    createNH2Above(midRow - 10);
    box[ index(curCol + 4, midRow - 7) ] = '|'; box[ index(curCol, midRow - 6) ] = 'O'; 
    box[ index(curCol + 2, midRow - 6) ] = '='; box[ index(curCol + 4, midRow - 6) ] = 'C';
    box[ index(curCol + 4, midRow - 5) ] = '|';
    createCH2Above(midRow - 4);

    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group
    createCH2Below(midRow + 4);
    box[ index(curCol + 4, midRow + 5) ] = '|'; box[ index(curCol + 4, midRow + 6) ] = 'C';
    box[ index(curCol + 6, midRow + 6) ] = '='; box[ index(curCol + 8, midRow + 6) ] = 'O';
    box[ index(curCol + 4, midRow + 7) ] = '|';
    createNH2Below(midRow + 8);
  } // else odd, R group below

  curCol += 12; // shift over current col
} // createN()


void Structure::createD(int position, float pH, float rpKa)
{
  if(position % 2 == 0)
  {
    // R group: (CH2)(C(=O)(OH))
    if((pH - rpKa) < 0)
    {
      box[ index(curCol + 4, midRow - 10) ] = 'H'; 
      box[ index(curCol + 4, midRow - 9) ] = '|';
    } // if pH < pKa, protonate R group
    
    else // pH > pKa, don't protonate and add (-)
      box[ index(curCol + 3, midRow - 9) ] = '-';

    box[ index(curCol + 4, midRow - 8) ] = 'O'; box[ index(curCol + 4, midRow - 7) ] = '|';
    box[ index(curCol + 4, midRow - 5) ] = '|'; box[ index(curCol, midRow - 6) ] = 'O';
    box[ index(curCol + 2, midRow - 6) ] = '='; box[ index(curCol + 4, midRow - 6) ] = 'C';
    createCH2Above(midRow - 4);

    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group (CH2)(C(=O)(OH))
    createCH2Below(midRow + 4);
    box[ index(curCol + 4, midRow + 5) ] = '|'; box[ index(curCol + 4, midRow + 6) ] = 'C';
    box[ index(curCol + 6, midRow + 6) ] = '='; box[ index(curCol + 8, midRow + 6) ] = 'O';
    box[ index(curCol + 4, midRow + 7) ] = '|'; box[ index(curCol + 4, midRow + 8) ] = 'O';

    if((pH - rpKa) < 0)
    {
      box[ index(curCol + 4, midRow + 9) ] = '|'; 
      box[ index(curCol + 4, midRow + 10) ] = 'H';
    } // if pH < pKa, protonate R group

    else // pH > pKa, don't protonate and add (-)
      box[ index(curCol + 5, midRow + 9) ] = '-';

  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createD()


void Structure::createC(int position, float pH, float rpKa)
{
  if(position % 2 == 0)
  {
    // R group: (CH2)(SH)
    if((pH - rpKa) < 0)
    {
      box[ index(curCol + 4, midRow - 8) ] = 'H'; 
      box[ index(curCol + 4, midRow - 7) ] = '|';
    } // if pH < pKa, protonate R group

    else // pH > pKa, don't protonate and add (-)
      box[ index(curCol + 3, midRow - 5) ] = '-';
    
    box[ index(curCol + 4, midRow - 6) ] = 'S';
    createCH2Above(midRow - 4);

    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group: (CH2)(SH)
    createCH2Below(midRow + 4);
    box[ index(curCol + 4, midRow + 6) ] = 'S'; 

    if((pH - rpKa) < 0)
    {
      box[ index(curCol + 4, midRow + 7) ] = '|'; 
      box[ index(curCol + 4, midRow + 8) ] = 'H';
    } // if pH < pKa, protonate R group
    
    else // pH > pKa, don't protonate and add (-)
      box[ index(curCol + 5, midRow + 7) ] = '-';

  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createC()


void Structure::createE(int position, float pH, float rpKa)
{
  if(position % 2 == 0)
  {
    // R group: (CH2)(CH2)((C=0)(OH))
    if((pH - rpKa) < 0)
    {
      box[ index(curCol + 4, midRow - 12) ] = 'H'; 
      box[ index(curCol + 4, midRow - 11) ] = '|';
    } // if pH < pKa, protonate R group
    
    else // pH > pKa, don't protonate and add (-)
      box[ index(curCol + 3, midRow - 11) ] = '-';

    box[ index(curCol + 4, midRow - 10) ] = 'O'; box[ index(curCol + 4, midRow - 9) ] = '|';
    box[ index(curCol, midRow - 8) ] = 'O';      box[ index(curCol + 2, midRow - 8) ] = '='; 
    box[ index(curCol + 4, midRow - 8) ] = 'C';  box[ index(curCol + 4, midRow - 7) ] = '|';
    createCH2Above(midRow - 6); createCH2Above(midRow - 4);
    
    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group:(CH2)(CH2)((C=0)(OH))
    createCH2Below(midRow + 4); createCH2Below(midRow + 6);
    box[ index(curCol + 4, midRow + 7) ] = '|'; box[ index(curCol + 4, midRow + 8) ] = 'C';
    box[ index(curCol + 6, midRow + 8) ] = '='; box[ index(curCol + 8, midRow + 8) ] = 'O';
    box[ index(curCol + 4, midRow + 9) ] = '|'; box[ index(curCol + 4, midRow + 10) ] = 'O';

    if((pH - rpKa) < 0)
    {
      box[ index(curCol + 4, midRow + 11) ] = '|'; 
      box[ index(curCol + 4, midRow + 12) ] = 'H';
    } // if pH < pKa, protonate R group
    
    else // pH > pKa, don't protonate and add (-)
      box[ index(curCol + 5, midRow + 11) ] = '-';
  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createE()


void Structure::createQ(int position)
{
  if(position % 2 == 0)
  {
    // R group: (CH2)(CH2)(C(=O)(NH2)) 
    createNH2Above(midRow - 12);
    box[ index(curCol + 4, midRow - 9) ] = '|'; box[ index(curCol, midRow - 8) ] = 'O'; 
    box[ index(curCol + 2, midRow - 8) ] = '='; box[ index(curCol + 4, midRow - 8) ] = 'C';
    box[ index(curCol + 4, midRow - 7) ] = '|';
    createCH2Above(midRow - 6); createCH2Above(midRow - 4);

    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group
    createCH2Below(midRow + 4); createCH2Below(midRow + 6);
    box[ index(curCol + 4, midRow + 7) ] = '|'; box[ index(curCol + 4, midRow + 8) ] = 'C';
    box[ index(curCol + 6, midRow + 8) ] = '='; box[ index(curCol + 8, midRow + 8) ] = 'O';
    box[ index(curCol + 4, midRow + 9) ] = '|';
    createNH2Below(midRow + 10);
  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createQ()


void Structure::createG(int position)
{
  if(position % 2 == 0)
  {
    // R group: H
    box[ index(curCol + 4, midRow - 5) ] = 'H'; 
    box[ index(curCol + 4, midRow - 4) ] = '|'; 

    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group: H
    box[ index(curCol + 4, midRow + 4) ] = '|'; 
    box[ index(curCol + 4, midRow + 5) ] = 'H';
  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createG()


void Structure::createH(int position, float pH, float rpKa)
{
  if(position % 2 == 0)
  {
    // R group: (CH2)C(CH)N(CH)(NH)
    if((pH - rpKa) < 0)
    {
      box[ index(curCol - 2, midRow - 12) ] = 'H'; 
      box[ index(curCol, midRow - 11) ] = '\\';
      box[ index(curCol + 2, midRow - 11) ] = '+';
    } // if pH < pKa, protonate R group
    
    box[ index(curCol + 10, midRow - 12) ] = 'H'; box[ index(curCol + 8, midRow - 11) ] = '/';
    box[ index(curCol + 1, midRow - 10) ] = 'N';  box[ index(curCol + 4, midRow - 10) ] = '=';
    box[ index(curCol + 7, midRow - 10) ] = 'C';  box[ index(curCol, midRow - 9) ] = '/';
    box[ index(curCol + 8, midRow - 9) ] = '\\';
    box[ index(curCol - 5, midRow - 8) ] = 'H';   box[ index(curCol - 3, midRow - 8) ] = '-'; 
    box[ index(curCol - 1, midRow - 8) ] = 'C';   box[ index(curCol + 9, midRow - 8) ] = 'N';
    box[ index(curCol + 11, midRow - 8) ] = '-';  box[ index(curCol + 13, midRow - 8) ] = 'H';
    box[ index(curCol, midRow - 7) ] = '\\';      box[ index(curCol + 1, midRow - 7) ] = '\\';
    box[ index(curCol + 7, midRow - 7) ] = '/';   box[ index(curCol + 4, midRow - 6) ] = 'C';
    createCH2Above(midRow - 4);

    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group
    createCH2Below(midRow + 4);
    box[ index(curCol + 4, midRow + 6) ] = 'C';   box[ index(curCol + 1, midRow + 7) ] = '/';
    box[ index(curCol + 7, midRow + 7) ] = '\\';  box[ index(curCol + 8, midRow + 7) ] = '\\';
    box[ index(curCol - 5, midRow + 8) ] = 'H';   box[ index(curCol - 3, midRow + 8) ] = '-'; 
    box[ index(curCol - 1, midRow + 8) ] = 'N';   box[ index(curCol + 9, midRow + 8) ] = 'C';
    box[ index(curCol + 11, midRow + 8) ] = '-';  box[ index(curCol + 13, midRow + 8) ] = 'H'; 
    box[ index(curCol, midRow + 9) ] = '\\';      box[ index(curCol + 8, midRow + 9) ] = '/';
    box[ index(curCol + 1, midRow + 10) ] = 'N'; 
    box[ index(curCol + 4, midRow + 10) ] = '=';  box[ index(curCol + 7, midRow + 10) ] = 'C';
    box[ index(curCol + 8, midRow + 11) ] = '\\'; box[ index(curCol + 10, midRow + 12) ] = 'H';

    if((pH - rpKa) < 0)
    {
      box[ index(curCol, midRow + 11) ] = '/'; 
      box[ index(curCol - 2, midRow + 12) ] = 'H';
      box[ index(curCol + 2, midRow + 11) ] = '+';
    } // if pH < pKa, protonate R group

  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createH()


void Structure::createI(int position)
{
  if(position % 2 == 0)
  {
    // R group: (CH(CH3))(CH2)(CH3)
    box[ index(curCol + 4, midRow - 10) ] = 'H';  createCH2Above(midRow - 8);
    box[ index(curCol - 3, midRow - 6) ] = 'H';   createCH2Above(midRow - 6);
    box[ index(curCol - 1, midRow - 5) ] = '\\';  box[ index(curCol + 4, midRow - 5) ] = '|';
    box[ index(curCol - 4, midRow - 4) ] = 'H';   box[ index(curCol - 2 , midRow - 4) ] = '-';
    box[ index(curCol, midRow - 4) ] = 'C';       box[ index(curCol + 2, midRow - 4) ] = '-';
    box[ index(curCol + 4, midRow - 4) ] = 'C';   box[ index(curCol + 6, midRow - 4) ] = '-';
    box[ index(curCol + 8, midRow - 4) ] = 'H';   box[ index(curCol, midRow - 3) ] = '|';
    box[ index(curCol, midRow - 2) ] = 'H';

    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group
    box[ index(curCol + 8, midRow + 2) ] = 'H';  box[ index(curCol + 8, midRow + 3) ] = '|';
    box[ index(curCol, midRow + 4) ] = 'H';      box[ index(curCol + 2, midRow + 4) ] = '-';
    box[ index(curCol + 4, midRow + 4) ] = 'C';  box[ index(curCol + 6, midRow + 4) ] = '-';
    box[ index(curCol + 8, midRow + 4) ] = 'C';  box[ index(curCol + 10, midRow + 4) ] = '-';
    box[ index(curCol + 12, midRow + 4) ] = 'H'; box[ index(curCol + 4, midRow + 5) ] = '|';
    box[ index(curCol + 9, midRow + 5) ] = '\\';
    createCH2Below(midRow + 6); box[ index(curCol + 11, midRow + 6) ] = 'H';
    createCH2Below(midRow + 8); box[ index(curCol + 4, midRow + 10) ] = 'H';
  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createI()


void Structure::createL(int position)
{
  if(position % 2 == 0)
  {
    // R group: (CH2)(CH(CH3))(CH3)
    box[ index(curCol + 4, midRow - 10) ] = 'H'; box[ index(curCol - 3, midRow - 8) ] = 'H';
    createCH2Above(midRow - 8);                  box[ index(curCol - 1, midRow - 7) ] = '\\';
    box[ index(curCol + 4, midRow - 7) ] = '|';  box[ index(curCol - 4, midRow - 6) ] = 'H'; 
    box[ index(curCol - 2 , midRow - 6) ] = '-'; box[ index(curCol, midRow - 6) ] = 'C'; 
    box[ index(curCol + 2, midRow - 6) ] = '-';  box[ index(curCol + 4, midRow - 6) ] = 'C'; 
    box[ index(curCol + 6, midRow - 6) ] = '-';  box[ index(curCol + 8, midRow - 6) ] = 'H';
    box[ index(curCol - 1, midRow - 5) ] = '/';  box[ index(curCol - 3, midRow - 4) ] = 'H';
    createCH2Above(midRow - 4);

    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group
    createCH2Below(midRow + 4);                  box[ index(curCol + 11, midRow + 4) ] = 'H';
    box[ index(curCol + 9, midRow + 5) ] = '/';  box[ index(curCol, midRow + 6) ] = 'H';
    box[ index(curCol + 2, midRow + 6) ] = '-';  box[ index(curCol + 4, midRow + 6) ] = 'C'; 
    box[ index(curCol + 6, midRow + 6) ] = '-';  box[ index(curCol + 8, midRow + 6) ] = 'C'; 
    box[ index(curCol + 10, midRow + 6) ] = '-'; box[ index(curCol + 12, midRow + 6) ] = 'H';
    box[ index(curCol + 9, midRow + 7) ] = '\\'; box[ index(curCol + 4, midRow + 7) ] = '|';
    box[ index(curCol + 11, midRow + 8) ] = 'H'; createCH2Below(midRow + 8);
    box[ index(curCol + 4, midRow + 10) ] = 'H';
  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createL()


void Structure::createK(int position, float pH, float rpKa)
{
  if(position % 2 == 0)
  {
    // R group: (CH2)(CH2)(CH2)(CH2)(NH2)
    if((pH - rpKa) < 0)
    {
      box[ index(curCol + 6, midRow - 13) ] = '+';
      box[ index(curCol + 6, midRow - 12) ] = '-'; 
      box[ index(curCol + 8, midRow - 12) ] = 'H';
    } // if pH < pKa, protonate R group

    createNH2Above(midRow - 14); 
    createCH2Above(midRow - 10); createCH2Above(midRow - 8);
    createCH2Above(midRow - 6); createCH2Above(midRow - 4);

    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group
    createCH2Below(midRow + 4); createCH2Below(midRow + 6);
    createCH2Below(midRow + 8); createCH2Below(midRow + 10);
    createNH2Below(midRow + 12);

    if((pH - rpKa) < 0)
    {
      box[ index(curCol, midRow + 12) ] = 'H'; 
      box[ index(curCol + 2, midRow + 12) ] = '-';
      box[ index(curCol + 2, midRow + 13) ] = '+';
    } // if pH < pKa, protonate R group
    
  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createK()


void Structure::createM(int position)
{
  if(position % 2 == 0)
  {
    // R group: (CH2)(CH2)(S)(CH3)
    box[ index(curCol + 4, midRow - 12) ] = 'H'; createCH2Above(midRow - 10);
    box[ index(curCol + 4, midRow - 9) ] = '|';  box[ index(curCol + 4, midRow - 8) ] = 'S';
    createCH2Above(midRow - 6); createCH2Above(midRow - 4); 

    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group
    createCH2Below(midRow + 4); createCH2Below(midRow + 6);
    box[ index(curCol + 4, midRow + 8) ] = 'S'; box[ index(curCol + 4, midRow + 9) ] = '|';
    createCH2Below(midRow + 10);                box[ index(curCol + 4, midRow + 12) ] = 'H';
  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createM()


void Structure::createF(int position)
{
  if(position % 2 == 0)
  {
    // R group: (CH2)(C6H6)
    box[ index(curCol - 5, midRow - 14) ] = 'H';  box[ index(curCol + 7, midRow - 14) ] = 'H';
    box[ index(curCol - 3, midRow - 13) ] = '\\'; box[ index(curCol + 5, midRow - 13) ] = '/';
    box[ index(curCol - 2, midRow - 12) ] = 'C';  box[ index(curCol + 1, midRow - 12) ] = '=';
    box[ index(curCol + 4, midRow - 12) ] = 'C';  box[ index(curCol - 3, midRow - 11) ] = '/';
    box[ index(curCol + 5, midRow - 11) ] = '\\'; box[ index(curCol - 8, midRow - 10) ] = 'H';
    box[ index(curCol - 6, midRow - 10) ] = '-';  box[ index(curCol - 4, midRow - 10) ] = 'C';  
    box[ index(curCol + 5, midRow - 10) ] = 'C';  box[ index(curCol + 7, midRow - 10) ] = '-';  
    box[ index(curCol + 9, midRow - 10) ] = 'H';  box[ index(curCol - 3, midRow - 9) ] = '\\';  
    box[ index(curCol - 2, midRow - 9) ] = '\\';  box[ index(curCol + 4, midRow - 9) ] = '/';   
    box[ index(curCol + 5, midRow - 9) ] = '/';   box[ index(curCol - 2, midRow - 8) ] = 'C';   
    box[ index(curCol + 1, midRow - 8) ] = '-';   box[ index(curCol + 4, midRow - 8) ] = 'C';
    box[ index(curCol - 3, midRow - 7) ] = '/';   box[ index(curCol - 5, midRow - 6) ] = 'H'; 
    createCH2Above(midRow - 6);
    box[ index(curCol + 4, midRow - 5) ] = '|';   box[ index(curCol + 4, midRow - 4) ] = '|';

    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group
    box[ index(curCol + 4, midRow + 4) ] = '|';   box[ index(curCol + 4, midRow + 5) ] = '|';
    createCH2Below(midRow + 6);
    box[ index(curCol + 13, midRow + 6) ] = 'H';  box[ index(curCol + 11, midRow + 7) ] = '/'; 
    box[ index(curCol + 4, midRow + 8) ] = 'C';   box[ index(curCol + 7, midRow + 8) ] = '-';
    box[ index(curCol + 10, midRow + 8) ] = 'C';  box[ index(curCol + 3, midRow + 9) ] = '/'; 
    box[ index(curCol + 4, midRow + 9) ] = '/';   box[ index(curCol + 10, midRow + 9) ] = '\\'; 
    box[ index(curCol + 11, midRow + 9) ] = '\\'; box[ index(curCol - 2, midRow + 10) ] = 'H'; 
    box[ index(curCol, midRow + 10) ] = '-';      box[ index(curCol + 2, midRow + 10) ] = 'C'; 
    box[ index(curCol + 11, midRow + 10) ] = 'C'; box[ index(curCol + 13, midRow + 10) ] = '-'; 
    box[ index(curCol + 15, midRow + 10) ] = 'H'; box[ index(curCol + 3, midRow + 11) ] = '\\'; 
    box[ index(curCol + 11, midRow + 11) ] = '/'; box[ index(curCol + 4, midRow + 12) ] = 'C'; 
    box[ index(curCol + 7, midRow + 12) ] = '=';  box[ index(curCol + 10, midRow + 12) ] = 'C';
    box[ index(curCol + 3, midRow + 13) ] = '/';  box[ index(curCol + 11, midRow + 13) ] = '\\';
    box[ index(curCol + 1, midRow + 14) ] = 'H';  box[ index(curCol + 13, midRow + 14) ] = 'H';
  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createF()


void Structure::createP(int position)
{
  if(position % 2 == 0)
  {
    // backbone; doesn't use generic backbone function
    box[ index(curCol + 10, midRow - 2) ] = 'H'; box[ index(curCol + 16, midRow - 2) ] = 'O';
    box[ index(curCol + 10, midRow - 1) ] = '|'; box[ index(curCol + 16, midRow - 1) ] = '|';
    box[ index(curCol + 17, midRow - 1) ] = '|'; 
    box[ index(curCol, midRow) ] = '-';          box[ index(curCol - 1, midRow) ] = '-';
    box[ index(curCol + 1, midRow) ] = '-';      box[ index(curCol + 2, midRow) ] = '-';
    box[ index(curCol + 4, midRow) ] = 'N';      box[ index(curCol + 6, midRow) ] = '-';
    box[ index(curCol + 7, midRow) ] = '-';      box[ index(curCol + 8, midRow) ] = '-';
    box[ index(curCol + 10, midRow) ] = 'C';     box[ index(curCol + 12, midRow) ] = '-';
    box[ index(curCol + 13, midRow) ] = '-';     box[ index(curCol + 14, midRow) ] = '-';
    box[ index(curCol + 16, midRow) ] = 'C';     box[ index(curCol + 18, midRow) ] = '-';
    box[ index(curCol + 19, midRow) ] = '-';     box[ index(curCol + 20, midRow) ] = '-';

    // R group
    box[ index(curCol + 4, midRow + 1) ] = '|';  box[ index(curCol + 10, midRow + 1) ] = '|';
    box[ index(curCol, midRow + 2) ] = 'H';      box[ index(curCol + 2, midRow + 2) ] = '-';
    box[ index(curCol + 4, midRow + 2) ] = 'C';  box[ index(curCol + 10, midRow + 2) ] = 'C';
    box[ index(curCol + 12, midRow + 2) ] = '-'; box[ index(curCol + 14, midRow + 2) ] = 'H';
    box[ index(curCol + 3, midRow + 3) ] = '/';  box[ index(curCol + 5, midRow + 3) ] = '\\';
    box[ index(curCol + 9, midRow + 3) ] = '/';  box[ index(curCol + 11, midRow + 3) ] = '\\';
    box[ index(curCol + 1, midRow + 4) ] = 'H';  box[ index(curCol + 7, midRow + 4) ] = 'C';
    box[ index(curCol + 13, midRow + 4) ] = 'H'; box[ index(curCol + 6, midRow + 5) ] = '/';
    box[ index(curCol + 8, midRow + 5) ] = '\\'; box[ index(curCol + 4, midRow + 6) ] = 'H';
    box[ index(curCol + 10, midRow + 6) ] = 'H';
  } // if even, R group below (special case)

  else // position % 2 != 0
  {
    // R group
    box[ index(curCol + 4, midRow - 6) ] = 'H';  box[ index(curCol + 10, midRow - 6) ] = 'H';
    box[ index(curCol + 6, midRow - 5) ] = '\\'; box[ index(curCol + 8, midRow - 5) ] = '/';
    box[ index(curCol + 1, midRow - 4) ] = 'H';  box[ index(curCol + 7, midRow - 4) ] = 'C';
    box[ index(curCol + 13, midRow - 4) ] = 'H'; box[ index(curCol + 3, midRow - 3) ] = '\\';
    box[ index(curCol + 5, midRow - 3) ] = '/';  box[ index(curCol + 9, midRow - 3) ] = '\\';
    box[ index(curCol + 11, midRow - 3) ] = '/'; box[ index(curCol, midRow - 2) ] = 'H';
    box[ index(curCol + 2, midRow - 2) ] = '-';  box[ index(curCol + 4, midRow - 2) ] = 'C';
    box[ index(curCol + 10, midRow - 2) ] = 'C'; box[ index(curCol + 12, midRow - 2) ] = '-';
    box[ index(curCol + 14, midRow - 2) ] = 'H'; box[ index(curCol + 4, midRow - 1) ] = '|';
    box[ index(curCol + 10, midRow - 1) ] = '|';
  
    // backbone
    box[ index(curCol - 1, midRow) ] = '-';      box[ index(curCol, midRow) ] = '-'; 
    box[ index(curCol + 1, midRow) ] = '-';      box[ index(curCol + 2, midRow) ] = '-'; 
    box[ index(curCol + 4, midRow) ] = 'N';      box[ index(curCol + 6, midRow) ] = '-'; 
    box[ index(curCol + 7, midRow) ] = '-';      box[ index(curCol + 8, midRow) ] = '-'; 
    box[ index(curCol + 10, midRow) ] = 'C';     box[ index(curCol + 12, midRow) ] = '-'; 
    box[ index(curCol + 13, midRow) ] = '-';     box[ index(curCol + 14, midRow) ] = '-'; 
    box[ index(curCol + 16, midRow) ] = 'C';     box[ index(curCol + 18, midRow) ] = '-';
    box[ index(curCol + 19, midRow) ] = '-';     box[ index(curCol + 20, midRow) ] = '-';
    box[ index(curCol + 10, midRow + 1) ] = '|'; box[ index(curCol + 16, midRow + 1) ] = '|';
    box[ index(curCol + 17, midRow + 1) ] = '|'; box[ index(curCol + 10, midRow + 2) ] = 'H'; 
    box[ index(curCol + 16, midRow + 2) ] = 'O';
  } // else odd, R group above (special case)

  curCol += 22; // shift over current col  
} // createP()


void Structure::createS(int position)
{
  if(position % 2 == 0)
  {
    // R group: (CH2)(OH)
    box[ index(curCol + 4, midRow - 8) ] = 'H'; box[ index(curCol + 4, midRow - 7) ] = '|';
    box[ index(curCol + 4, midRow - 6) ] = 'O'; createCH2Above(midRow - 4);

    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group
    createCH2Below(midRow + 4);                 box[ index(curCol + 4, midRow + 6) ] = 'O';
    box[ index(curCol + 4, midRow + 7) ] = '|'; box[ index(curCol + 4, midRow + 8) ] = 'H';
  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createS()


void Structure::createT(int position)
{
  if(position % 2 == 0)
  {
    // R group: (CH(OH))(CH3)
    box[ index(curCol + 4, midRow - 8) ] = 'H'; createCH2Above(midRow - 6);
    box[ index(curCol + 4, midRow - 5) ] = '|'; box[ index(curCol - 4, midRow - 4) ] = 'H';
    box[ index(curCol - 2, midRow - 4) ] = '-'; box[ index(curCol, midRow - 4) ] = 'O';
    box[ index(curCol + 2, midRow - 4) ] = '-'; box[ index(curCol + 4, midRow - 4) ] = 'C';
    box[ index(curCol + 6, midRow - 4) ] = '-'; box[ index(curCol + 8, midRow - 4) ] = 'H';
    
    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group
    box[ index(curCol, midRow + 4) ] = 'H';      box[ index(curCol + 2, midRow + 4) ] = '-';
    box[ index(curCol + 4, midRow + 4) ] = 'C';  box[ index(curCol + 6, midRow + 4) ] = '-';
    box[ index(curCol + 8, midRow + 4) ] = 'O';  box[ index(curCol + 10, midRow + 4) ] = '-';
    box[ index(curCol + 12, midRow + 4) ] = 'H'; box[ index(curCol + 4, midRow + 5) ] = '|';
    createCH2Below(midRow + 6);                  box[ index(curCol + 4, midRow + 8) ] = 'H';
  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createT()


void Structure::createW(int position)
{
  if(position % 2 == 0)
  {
    // R group: (phenyl)(NH)(CH=C)(CH2)
    box[ index(curCol - 9, midRow - 16) ] = 'H';  box[ index(curCol + 3, midRow - 16) ] = 'H';
    box[ index(curCol - 7, midRow - 15) ] = '\\'; box[ index(curCol + 1, midRow - 15) ] = '/';
    box[ index(curCol - 6, midRow - 14) ] = 'C';  box[ index(curCol - 3, midRow - 14) ] = '=';
    box[ index(curCol, midRow - 14) ] = 'C';      box[ index(curCol + 10, midRow - 14) ] = 'H';
    box[ index(curCol - 7, midRow - 13) ] = '/';  box[ index(curCol + 1, midRow - 13) ] = '\\';
    box[ index(curCol + 8, midRow - 13) ] = '/';  box[ index(curCol - 12, midRow - 12) ] = 'H';
    box[ index(curCol - 10, midRow - 12) ] = '-'; box[ index(curCol - 8, midRow - 12) ] = 'C';  
    box[ index(curCol + 2, midRow - 12) ] = 'C';  box[ index(curCol + 4, midRow - 12) ] = '-';  
    box[ index(curCol + 5, midRow - 12) ] = '-';  box[ index(curCol + 7, midRow - 12) ] = 'N'; 
    box[ index(curCol - 7, midRow - 11) ] = '\\'; box[ index(curCol - 6, midRow - 11) ] = '\\';
    box[ index(curCol, midRow - 11) ] = '/';      box[ index(curCol + 1, midRow - 11) ] = '/';
    box[ index(curCol + 8, midRow - 11) ] = '\\'; box[ index(curCol - 6, midRow - 10) ] = 'C'; 
    box[ index(curCol - 3, midRow - 10) ] = '-';  box[ index(curCol, midRow - 10) ] = 'C'; 
    box[ index(curCol + 9, midRow - 10) ] = 'C';  box[ index(curCol - 7, midRow - 9) ] = '/'; 
    box[ index(curCol + 1, midRow - 9) ] = '\\';  box[ index(curCol + 7, midRow - 9) ] = '/'; 
    box[ index(curCol + 8, midRow - 9) ] = '/';   box[ index(curCol - 8, midRow - 8) ] = 'H'; 
    box[ index(curCol + 4, midRow - 8) ] = 'C';

    // extends R group so that it doesn't overlap with other R groups
    for(int i = 0; i < 4; i++)
      box[ index(curCol + 4, midRow - 7 + i) ] = '|';

    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group
    // extends R group so that it doesn't overlap with other R groups
    for(int i = 0; i < 4; i++)
      box[ index(curCol + 4, midRow + 4 + i) ] = '|';

    box[ index(curCol + 4, midRow + 8) ] = 'C';    box[ index(curCol + 17, midRow + 8) ] = 'H';
    box[ index(curCol, midRow + 9) ] = '/';        box[ index(curCol + 1, midRow + 9) ] = '/';
    box[ index(curCol + 7, midRow + 9) ] = '\\';   box[ index(curCol + 15, midRow + 9) ] = '/';
    box[ index(curCol - 5, midRow + 10) ] = 'H';   box[ index(curCol - 3, midRow + 10) ] = '-';
    box[ index(curCol - 1, midRow + 10) ] = 'C';   box[ index(curCol + 8, midRow + 10) ] = 'C'; 
    box[ index(curCol + 11, midRow + 10) ] = '-';  box[ index(curCol + 14, midRow + 10) ] = 'C'; 
    box[ index(curCol, midRow + 11) ] = '\\';      box[ index(curCol + 7, midRow + 11) ] = '/'; 
    box[ index(curCol + 8, midRow + 11) ] = '/';   box[ index(curCol + 14, midRow + 11) ] = '\\';
    box[ index(curCol + 15, midRow + 11) ] = '\\'; box[ index(curCol + 1 , midRow + 12) ] = 'N';
    box[ index(curCol + 4, midRow + 12) ] = '-';   box[ index(curCol + 7, midRow + 12) ] = 'C'; 
    box[ index(curCol + 16, midRow + 12) ] = 'C';  box[ index(curCol + 18, midRow + 12) ] = '-'; 
    box[ index(curCol + 20, midRow + 12) ] = 'H';  box[ index(curCol, midRow + 13) ] = '/';
    box[ index(curCol + 7, midRow + 13) ] = '\\';  box[ index(curCol + 16, midRow + 13) ] = '/';
    box[ index(curCol - 2, midRow + 14) ] = 'H';   box[ index(curCol + 8, midRow + 14) ] = 'C'; 
    box[ index(curCol + 11, midRow + 14) ] = '=';  box[ index(curCol + 14, midRow + 14) ] = 'C';
    box[ index(curCol + 7, midRow + 15) ] = '/';   box[ index(curCol + 15, midRow + 15) ] = '\\';
    box[ index(curCol + 5, midRow + 16) ] = 'H';   box[ index(curCol + 17, midRow + 16) ] = 'H';
  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createW()


void Structure::createY(int position, float pH, float rpKa)
{
  if(position % 2 == 0)
  {
    // R group: (CH2)(phenol)
    if((pH - rpKa) < 0) // if pH < pKa, protonate R group
      box[ index(curCol - 6, midRow - 15) ] = 'H';
    
    else // pH > pKa, don't protonate and add (-)
      box[ index(curCol - 6, midRow - 16) ] = '-';

    box[ index(curCol - 5, midRow - 15) ] = 'O';  box[ index(curCol + 7, midRow - 15) ] = 'H';
    box[ index(curCol - 3, midRow - 14) ] = '\\'; box[ index(curCol + 5, midRow - 14) ] = '/';
    box[ index(curCol - 2, midRow - 13) ] = 'C';  box[ index(curCol + 1, midRow - 13) ] = '-';
    box[ index(curCol + 4, midRow - 13) ] = 'C';  box[ index(curCol - 3, midRow - 12) ] = '/';  
    box[ index(curCol + 5, midRow - 12) ] = '\\'; box[ index(curCol - 8, midRow - 11) ] = 'H';  
    box[ index(curCol - 6, midRow - 11) ] = '-';  box[ index(curCol - 4, midRow - 11) ] = 'C';  
    box[ index(curCol + 5, midRow - 11) ] = 'C';  box[ index(curCol + 7, midRow - 11) ] = '-';  
    box[ index(curCol + 9, midRow - 11) ] = 'H';  box[ index(curCol - 3, midRow - 10) ] = '\\'; 
    box[ index(curCol - 2, midRow - 10) ] = '\\'; box[ index(curCol + 4, midRow - 10) ] = '/';  
    box[ index(curCol + 5, midRow - 10) ] = '/';  box[ index(curCol - 2, midRow - 9) ] = 'C';   
    box[ index(curCol + 1, midRow - 9) ] = '-';   box[ index(curCol + 4, midRow - 9) ] = 'C';
    box[ index(curCol - 3, midRow - 8) ] = '/';   box[ index(curCol - 5, midRow - 7) ] = 'H'; 
    createCH2Above(midRow - 7);

    // extends R group so that it doesn't overlap with other R groups
    for(int i = 0; i < 3; i++)
      box[ index(curCol + 4, midRow - 6 + i) ] = '|';

    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group
    // extends R group so that it doesn't overlap with other R groups
    for(int i = 0; i < 3; i++)
      box[ index(curCol + 4, midRow + 4 + i) ] = '|';

    createCH2Below(midRow + 7);
    box[ index(curCol + 13, midRow + 7) ] = 'H';   box[ index(curCol + 11, midRow + 8) ] = '/'; 
    box[ index(curCol + 4, midRow + 9) ] = 'C';    box[ index(curCol + 7, midRow + 9) ] = '-';
    box[ index(curCol + 10, midRow + 9) ] = 'C';   box[ index(curCol + 3, midRow + 10) ] = '/'; 
    box[ index(curCol + 4, midRow + 10) ] = '/';   box[ index(curCol + 10, midRow + 10) ] = '\\'; 
    box[ index(curCol + 11, midRow + 10) ] = '\\'; box[ index(curCol - 2, midRow + 11) ] = 'H';   
    box[ index(curCol, midRow + 11) ] = '-';       box[ index(curCol + 2, midRow + 11) ] = 'C';   
    box[ index(curCol + 11, midRow + 11) ] = 'C';  box[ index(curCol + 13, midRow + 11) ] = '-';  
    box[ index(curCol + 15, midRow + 11) ] = 'H';  box[ index(curCol + 3, midRow + 12) ] = '\\';  
    box[ index(curCol + 11, midRow + 12) ] = '/';  box[ index(curCol + 4, midRow + 13) ] = 'C';   
    box[ index(curCol + 7, midRow + 13) ] = '=';   box[ index(curCol + 10, midRow + 13) ] = 'C';
    box[ index(curCol + 3, midRow + 14) ] = '/';   box[ index(curCol + 11, midRow + 14) ] = '\\';
    box[ index(curCol + 1, midRow + 15) ] = 'H';   box[ index(curCol + 13, midRow + 15) ] = 'O';

    if((pH - rpKa) < 0) // if pH < pKa, protonate R group
      box[ index(curCol + 14, midRow + 15) ] = 'H';
    
    else // pH > pKa, don't protonate and add (-)
      box[ index(curCol + 14, midRow + 16) ] = '-';
    
  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createY()


void Structure::createV(int position)
{
  if(position % 2 == 0)
  {
    // R group: (CH(CH3))(CH3)
    box[ index(curCol + 4, midRow - 8) ] = 'H';   createCH2Above(midRow - 6);
    box[ index(curCol - 3, midRow - 6) ] = 'H';   box[ index(curCol - 1, midRow - 5) ] = '\\'; 
    box[ index(curCol + 4, midRow - 5) ] = '|';   box[ index(curCol - 4, midRow - 4) ] = 'H'; 
    box[ index(curCol - 2 , midRow - 4) ] = '-';  box[ index(curCol, midRow - 4) ] = 'C';  
    box[ index(curCol + 2, midRow - 4) ] = '-';   box[ index(curCol + 4, midRow - 4) ] = 'C'; 
    box[ index(curCol + 6, midRow - 4) ] = '-';   box[ index(curCol + 8, midRow - 4) ] = 'H'; 
    box[ index(curCol, midRow - 3) ] = '|';       box[ index(curCol, midRow - 2) ] = 'H';
    
    // backbone
    createBackbone1();
  } // if even, R group above

  else // position % 2 != 0
  {
    // backbone
    createBackbone2();

    // R group
    box[ index(curCol + 8, midRow + 2) ] = 'H';  box[ index(curCol + 8, midRow + 3) ] = '|';
    box[ index(curCol, midRow + 4) ] = 'H';      box[ index(curCol + 2, midRow + 4) ] = '-';
    box[ index(curCol + 4, midRow + 4) ] = 'C';  box[ index(curCol + 6, midRow + 4) ] = '-';
    box[ index(curCol + 8, midRow + 4) ] = 'C';  box[ index(curCol + 10, midRow + 4) ] = '-';
    box[ index(curCol + 12, midRow + 4) ] = 'H'; box[ index(curCol + 4, midRow + 5) ] = '|';
    box[ index(curCol + 9, midRow + 5) ] = '\\'; box[ index(curCol + 11, midRow + 6) ] = 'H';
    createCH2Below(midRow + 6);                  box[ index(curCol + 4, midRow + 8) ] = 'H';
  } // else odd, R group below

  curCol += 12; // shift over current col  
} // createV()


// create a CH2 group above the backbone
void Structure::createCH2Above(int row)
{
  box[ index(curCol + 4, row - 1) ] = '|'; box[ index(curCol, row) ] = 'H'; 
  box[ index(curCol + 2, row) ] = '-';     box[ index(curCol + 4, row) ] = 'C';
  box[ index(curCol + 6, row) ] = '-';     box[ index(curCol + 8, row) ] = 'H';
} // createCH2Above()


// creates a CH2 group below the backbone
void Structure::createCH2Below(int row)
{
  box[ index(curCol, row) ] = 'H';     box[ index(curCol + 2, row) ] = '-';
  box[ index(curCol + 4, row) ] = 'C'; box[ index(curCol + 6, row) ] = '-';
  box[ index(curCol + 8, row) ] = 'H'; box[ index(curCol + 4, row + 1) ] = '|';
} // createCH2Below()


// creates an NH2 group above the backbone
void Structure::createNH2Above(int row) // row starts at top H
{
  box[ index(curCol + 4, row) ] = 'H'; box[ index(curCol + 4, row + 1) ] = '|';
  box[ index(curCol, row + 2) ] = 'H'; box[ index(curCol + 2, row + 2) ] = '-';
  box[ index(curCol + 4, row + 2) ] = 'N';
} // create NH2Above()


// creates an NH2 group below the backbone
void Structure::createNH2Below(int row) // row starts at N
{
  box[ index(curCol + 4, row) ] = 'N'; box[ index(curCol + 6, row) ] = '-';
  box[ index(curCol + 8, row) ] = 'H'; box[ index(curCol + 4, row + 1) ] = '|';
  box[ index(curCol + 4, row + 2) ] = 'H';
} // createNH2Below()


// creates the structure based on the amino acid in chain
void Structure::createAAStructure(AminoAcid **chain, float pH, int length)
{
  if(length == 0) // if chain is empty
    return;
  
  // inserting the amino terminal (NH3+)
  if( (pH - chain[0]->getNpKa()) < 0 )
  {
    box[ index(1, midRow - 2) ] = 'H';
    box[ index(3, midRow - 1) ] = '\\'; 
    box[ index(5, midRow - 1) ] = '+'; 
  } // if pH < N-terminal pKa

  box[ index(0, midRow) ] = 'H'; box[ index(2, midRow) ] = '-';

  curCol = 4;

  for(int pos = 0; pos < length; pos++)
  {
    switch(chain[pos]->getSymbol())
    {
      case 'A': createA(pos); break;
      case 'R': createR( pos, pH, chain[pos]->getRpKa() ); break;
      case 'N': createN(pos); break;
      case 'D': createD( pos, pH, chain[pos]->getRpKa() ); break;
      case 'C': createC( pos, pH, chain[pos]->getRpKa() ); break;
      case 'E': createE( pos, pH, chain[pos]->getRpKa() ); break;
      case 'Q': createQ(pos); break;
      case 'G': createG(pos); break;
      case 'H': createH( pos, pH, chain[pos]->getRpKa() ); break;
      case 'I': createI(pos); break;
      case 'L': createL(pos); break;
      case 'K': createK( pos, pH, chain[pos]->getRpKa() ); break;
      case 'M': createM(pos); break;
      case 'F': createF(pos); break;
      case 'P': createP(pos); break;
      case 'S': createS(pos); break;
      case 'T': createT(pos); break;
      case 'W': createW(pos); break;
      case 'Y': createY( pos, pH, chain[pos]->getRpKa() ); break;
      case 'V': createV(pos); break;
      case 'X': break; // STOP codon
    } // switch

  } // for each amino acid in chain

  // insert the carboxyl terminal
  box[ index(curCol, midRow) ] = 'O';

  if( (pH - chain[length - 1]->getCpKa()) < 0 )
  {
    box[ index(curCol + 2, midRow) ] = '-'; 
    box[ index(curCol + 4, midRow) ] = 'H';
  } // if pH < C-terminal pKa, protonate it
  
  else // pH > C-terminal pKa, don't protonate it
    box[ index(curCol + 1, midRow - 1) ] = '-';
} // createAAStructure()


void Structure::print() const
{
  for(int i = 0; i < size; i++)
  {
    if(i % width == 0)
      cout << endl;

    cout << box[i];
  } // for each square in grid
} // print()


// function doesn't get called but this function can be used to print a structure
// of choice. Was used to test the printing of structures
void Structure::testPrint()
{
  box[ index(1, midRow - 2) ] = 'H';
  box[ index(3, midRow - 1) ] = '\\'; box[ index(5, midRow - 1) ] = '+';
  box[ index(0, midRow) ] = 'H';      box[ index(2, midRow) ] = '-';
  curCol = 4;

  createR(0, 7.0, 12.48);
  createR(1, 13.0, 12.48);
  createD(0, 7.0, 3.65); 
  createD(1, 2.49, 3.65);
  createC(0, 7.0, 8.18); 
  createC(1, 9.73, 8.18);
  createE(0, 7.0, 4.25); 
  createE(1, 4.1, 4.25);
  createH(0, 7.0, 6.00); 
  createH(1, 5.96, 6.00);
  createK(0, 7.0, 10.53); 
  createK(1, 11.09, 10.53);
  createY(0, 7.0, 10.07); 
  createY(1, 10.08, 10.07);

  createE(0, 7.0, 4.25); 
  //createE(1, 4.1, 4.25);
  //createH(0, 7.0, 6.00); 
  createH(1, 5.96, 6.00);
  //createK(0, 7.0, 10.53); 
  createK(1, 11.09, 10.53);
  createY(0, 7.0, 10.07); 
  //createY(1, 10.08, 10.07);

  box[ index(curCol, midRow) ] = 'O'; box[ index(curCol + 1, midRow - 1) ] = '-'; 
  
  for(int i = 0; i < size; i++)
  {
    if(i % width == 0)
      cout << endl;

    cout << box[i];
  } // for

} // testPrint()
