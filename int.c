#include <math.h>
#include "int.h"
#include "basis.h"
float overlap_elem(char i,char j){
  return powf(M_PI/(get_alpha(i) + get_alpha(j)), 1.5);
}
float coulomb_elem(char i,char j){
  return -((2.0 * M_PI) / (get_alpha(i) + get_alpha(j)));
}
float kinetic_elem(char i,char j){
  return 3.0*((get_alpha(i) * get_alpha(j) * powf(M_PI,1.5)) / (powf((get_alpha(i) + get_alpha(j)), 2.5)));
}
float hamiltonian_elem(char i,char j){
  return coulomb_elem(i,j)+kinetic_elem(i,j);
}
