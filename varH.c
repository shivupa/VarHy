/* varH.c Calculates the variational ground state of the Hydrogen atom */
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_eigen.h>
#include "int.h"
#include "basis.h"
char num_eigenvalues = 4;
double hamiltonian[4*4];
double overlap[4*4];
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
      printf("Energy level %d = %g\n", abs(i-num_eigenvalues) , energy_i);
  }
  printf("Coefficent Matrix for plotting\n");
  for(char i=num_eigenvalues-1;i>=0;i--){
    for(char j=num_eigenvalues-1;j>=0;j--){
        printf("%g\t", gsl_matrix_get(coefficents,i,j));
    }
    printf("\n");
  }
  printf("Basis set\n");
  for(char i=num_eigenvalues;i>0;i--){
        printf("%g\t", get_alpha(abs(i-num_eigenvalues)));
  }
  printf("\n");
}
int main(){
  printf("***********************\n");
  printf("Variation Hydrogen Atom\n   By: Shiv Upadhyay\n");
  printf("***********************\n");
  variational_H();
  return 0;
}
