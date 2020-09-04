#include "msolver.hpp"
#include <iostream>
#include <fstream>
#include "time.h"

//f(x) test function
double matrix_solver::f(double x){
  return 100 * exp(-10*x);
}

double matrix_solver::analytic(double x){
   return 1 - (1 - exp(-10))*x - exp(-10*x);
}

//Setting up random a,b,c vectors
void matrix_solver::random_vectors(){
  //Initializing the uniformly randomly generated a,b,c vectors
  a_vec = arma::randu(n);
  b_vec = arma::randu(n);
  c_vec = arma::randu(n);
  }

  void matrix_solver::specified_vectors(int a_,int b_,int c_){
    a = a_; b = b_; c = c_;
    //Initializing the specified vectors
    a_vec = arma::vec(n);
    b_vec = arma::vec(n);
    c_vec = arma::vec(n);
    for(int i=0;i<n;i++){
      a_vec(i) = a;
      b_vec(i) = b;
      c_vec(i) = c;
    }
    }


//General forward solver
void matrix_solver::forward_solver_general(){
  //Setting inital values for the tilde vectors
  b_tilde(0) = b_vec(0); c_tilde(0) = c_vec(0); g_tilde(0) = g_vec(0);

  //Starting clock
  runtime = 0;
  start = clock();

  //Solving the forward algorythm
  for(int i=1;i<n;i++){
    b_tilde(i) = b_vec(i)*b_tilde(i-1)*1./a_vec(i-1) - c_tilde(i-1);
    c_tilde(i) = c_vec(i)*b_tilde(i-1)*1./a_vec(i-1);
    g_tilde(i) = g_vec(i)*b_tilde(i-1)*1./a_vec(i-1) - g_tilde(i-1);
  }
  //ending clock and adding to runtime
  finish = clock();
  runtime += ( (finish - start)*1./CLOCKS_PER_SEC );
}

void matrix_solver::forward_solver_specialized(){
  //Setting inital values for the tilde vectors
  b_tilde(0) = b_vec(0); c_tilde(0) = c_vec(0); g_tilde(0) = g_vec(0);

  //Starting clock
  runtime = 0;
  start = clock();
  //Solving the forward algorythm
  for(int i=1;i<n;i++){
    c_tilde(i) = b_tilde(i-1);
    b_tilde(i) = -2*b_tilde(i-1)-c_tilde(i-1);
    g_tilde(i) = -g_vec(i)*b_tilde(i-1)-g_tilde(i-1);

  }
  //ending clock and adding to runtime
  finish = clock();
  runtime += ( (finish - start)*1./CLOCKS_PER_SEC );
}

//General backward_solver
void matrix_solver::backward_solver(){
  //Setting up the endpoint
  v_vec(n) = g_tilde(n-1)*1./b_tilde(n-1);

  //Starting clock
  start = clock();
  //Solving the backward soloution
  for(int i=2;i<n+1;i++){
    int j = (n-i);
    v_vec(j+1) = (g_tilde(j) - c_tilde(j) *v_vec(j+2))*1./b_tilde(j);
  }
  //ending clock and adding to runtime
  finish = clock();
  runtime += ( (finish - start)*1./CLOCKS_PER_SEC );
}

//Writing x,u to outfile
void matrix_solver::write_file(std::string filename_){
  filename = filename_;
  std::ofstream outfile (filename);
  outfile << "| N=" << n<< " | epsilon_max=" << epsilon_max<<" | epsilon_tot="<<epsilon
  <<" | CPU-time:" << runtime << " | log_10(h)=" << log10(h) << " |\n";

  outfile << "| x | v | u |\n ";
  for(int i=0;i<n+2;i++){
    outfile << x_vec(i) << " " << v_vec(i) << " " << u_vec(i) << "\n";
  }
  outfile.close();
}

void matrix_solver::error(){
  epsilon = 0;
  epsilon_max = log10(abs((1./u_vec(5))*(v_vec(5)-u_vec(5))));
  for(int i=1;i<n+1;i++){
    epsilon += log10(abs((1./u_vec(i))*(v_vec(i)-u_vec(i))));
  }
}

double matrix_solver::ret_eps_max(){
  return epsilon_max;
}

double matrix_solver::ret_h(){
  return h;
}

double matrix_solver::ret_runtime(){
  return runtime;
}

void matrix_solver::LU_decomp(int a_, int b_, int c_, std::string filename_){
  a = a_; b = b_; c = c_;
  A = arma::mat(n,n,arma::fill::zeros);
  P = arma::mat(n,n,arma::fill::zeros);
  L = arma::mat(n,n,arma::fill::zeros);
  U = arma::mat(n,n,arma::fill::zeros);


  for(int i=0;i<n;i++){
    if (i!=0){
        A(i,i-1) = a;
      }
    A(i,i) = b;
    if(i!=n-1){
        A(i,i+1) = c;
      }
    }
    //Starting clock
    runtime = 0;
    start = clock();

    arma::lu(L,U,P,A);
    arma::vec b_LU;
    b_LU = arma::solve(trimatu(U), arma::solve(trimatl(L), P*g_vec));

    //ending clock and adding to runtime
    finish = clock();
    runtime += ( (finish - start)*1./CLOCKS_PER_SEC );

    for(int i=1;i<n+1;i++){
      v_vec(i) = b_LU(i-1);
    }


    filename = filename_;
    std::ofstream outfile (filename);

    outfile << "| x | v | u | CPU-time:" << runtime << std::endl;
    for(int i=0;i<n+2;i++){
      outfile << x_vec(i) << " " << v_vec(i) << " " << u_vec(i) << "\n";
    }
    outfile.close();
    runtime = 0;
}
