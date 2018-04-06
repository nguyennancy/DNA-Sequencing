#include "LinearProbing.h"


AAInfo::AAInfo(char let, string threeLet, string fullNme, float cP, float tP, 
               float rP): letter(let), cpKa(cP), npKa(tP), rpKa(rP)
{
  ID = int(letter) - 16; // convert letters to number equivalents
  threeLetter.append(threeLet);
  fullName.append(fullNme);
} // AAInfo()


AAInfo& AAInfo::operator= (const AAInfo &rhs)
{
  if(&rhs == this)
    return *this;

  letter = rhs.letter;
  ID = rhs.ID;
  threeLetter.append(rhs.threeLetter);
  fullName.append(rhs.fullName);
  cpKa = rhs.cpKa;
  npKa = rhs.npKa;
  rpKa = rhs.rpKa;

  return *this;
} // operator= ()


// inserts entry based on the hash fcn: ID % size
void LinearHashTable::insert( const AAInfo & x )
{
  int pos = x.ID % size;

  while(array[pos]->ID > 0 && array[pos]->ID != x.ID)
    if(++pos == size)
      pos = 0;

  *array[pos] = x;

} // insert()


LinearHashTable::LinearHashTable(): size(23) // next prime number after 20
{
  array = new AAInfo*[size];

  for(int i = 0; i < size; i++)
  {
    array[i] = new AAInfo(' ', "", "", 0.0, 0.0, 0.0);
    array[i]->ID = 0;
  } // for
} // LinearHashTable()


void LinearHashTable::create()
{
  char *lettr;
  string s, num, buffer, abbrev, full;
  float cPka, nPka, rPka;
  int hasR;

  ifstream inf("AminoAcids.txt");
  getline(inf, buffer); // gets header

  while(inf)
  {
    getline(inf, full, ',');
    getline(inf, abbrev, ',');
    getline(inf, s, ',');
    lettr = (char*) s.c_str();
    getline(inf, num, ',');
    cPka = stof(num);
    getline(inf, num, ',');
    nPka = stof(num);
    getline(inf, num, ',');
    hasR = stoi(num);

    if(hasR == 0)
    {
      rPka = 0.00;
      getline(inf, buffer, '\n'); // read in '\n'
    } // if amino acid's r-group doesn't have a pKa

    else // hasR == 1
    {
      getline(inf, buffer, '\n');
      rPka = stof(buffer);
    } // if amino acid's r-group has a pKa


    AAInfo temp((char)lettr[0], abbrev, full, cPka, nPka, rPka);
    insert(temp);
    count++;
  } //while there are lines to read in the file

  inf.close();
} // create()


// returns a stop codon if isn't found; but should always be able to find the correct
// AAInfo
const AAInfo& LinearHashTable::find(char symbol)
{
  int tempID = (int(symbol) - 16);
  int pos = tempID % size;

  while(array[pos]->ID != 0 && array[pos]->ID != tempID)
    if(++pos == size)
      pos = 0;

  if(array[pos]->ID == tempID)
    return *array[pos];

  else 
  {
    int stopID = (int('X') - 16);
    return *array[stopID % size];
  } // else isn't found

} // find()


// function isn't called; can be used to print out the LinearHashTable if user
// wishes to
void LinearHashTable::print() const
{
  for(int i = 0; i < 23; i++)
  {
    cout << setw(16) << setfill(' ') << left << array[i]->fullName
    << setw(6) << array[i]->letter << setw(6) << array[i]->threeLetter 
    << array[i]->cpKa << " " << array[i]->npKa << " " 
    << array[i]->rpKa << endl;
  } // for
} // print()

