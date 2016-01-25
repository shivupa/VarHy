/* varH.c Calculates the variational ground state of the Hydrogen atom */
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_eigen.h>
char num_eigenvalues = 4;
float alpha[4] = {13.00773, 1.962079, 0.444529, 0.1219492};
double hamiltonian[4*4];
/*float hamiltonian[num_eigenvalues][num_eigenvalues];
float overlap[num_eigenvalues][num_eigenvalues];*/
double overlap[4*4];
float overlap_elem(char i,char j){
  return powf(M_PI/(alpha[i] + alpha[j]), 1.5);
}
float coulomb_elem(char i,char j){
  return -((2.0 * M_PI) / (alpha[i] + alpha[j]));
}
float kinetic_elem(char i,char j){
  return 3.0*((alpha[i] * alpha[j] * powf(M_PI,1.5)) / (powf((alpha[i] + alpha[j]), 2.5)));
}
float hamiltonian_elem(char i,char j){
  return coulomb_elem(i,j)+kinetic_elem(i,j);
}
void variational_H(){
  /* Fill matricies */
  gsl_matrix *coefficents = gsl_matrix_alloc (num_eigenvalues, num_eigenvalues);
  gsl_vector *energy = gsl_vector_alloc (num_eigenvalues);
  for(char i=num_eigenvalues;i>=0;i--){
    for(char j=num_eigenvalues;j>=i;j--){
      overlap[num_eigenvalues*i + j] = overlap_elem(i,j);
      overlap[num_eigenvalues*j + i] = overlap[num_eigenvalues*i + j];
      hamiltonian[num_eigenvalues*i + j] = hamiltonian_elem(i,j);
      hamiltonian[num_eigenvalues*j + i] = hamiltonian[num_eigenvalues*i + j];
    }
  }
  gsl_matrix_view s = gsl_matrix_view_array (overlap, 4, 4);
  gsl_matrix_view h = gsl_matrix_view_array (hamiltonian, 4, 4);
  gsl_eigen_gensymmv_workspace * eigen_sys = gsl_eigen_gensymmv_alloc(num_eigenvalues);
  gsl_eigen_gensymmv(&h.matrix, &s.matrix, energy, coefficents, eigen_sys);
  gsl_eigen_gensymmv_free(eigen_sys);

	gsl_eigen_gensymmv_sort(energy, coefficents, GSL_EIGEN_SORT_VAL_DESC);
  for(char i=num_eigenvalues-1;i>=0;i--){
    double energy_i = gsl_vector_get (energy, i);
      printf("eigenvalue = %g\n", energy_i);
  }
}
int main(){
  printf("***********************\n");
  printf("Variation Hydrogen Atom\n   By: Shiv Upadhyay\n");
  printf("***********************\n");
  variational_H();
  return 0;
}
