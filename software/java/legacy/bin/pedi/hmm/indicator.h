

#ifndef INDICATOR
#define INDICATOR



class Indicator
{
 public:
  Indicator(int i, int p, double gf, double gm)
    {
      indiv_index = i;
      parent = p;
      grandfather = gf;
      grandmother = gm;
    }

  int indiv_index;
  int parent;
  double grandfather;
  double grandmother;
};


#endif
