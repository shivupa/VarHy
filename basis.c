#include <math.h>
#include "basis.h"
float alpha[4];
float get_alpha(char i){
  return alpha[i];
}
void set_alpha(char i,float val){
    alpha[i] = val;
}
