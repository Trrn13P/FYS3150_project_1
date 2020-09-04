#include <armadillo>
#include <iostream>
#include <cmath>
#include "msolver.cpp"
#include "time.h"
#include <fstream>

void write_stats(std::string filename_,int n_start, int n_end,std::string type_, int a, int b, int c){
  std::string filename = filename_;
  std::string type = type_;
  std::ofstream outfile (filename);
  int end = n_end;
  int start = n_start;
  for(int i=start;i<=end;i++){
    int n = int(pow(10,i));
    matrix_solver run_(n);
    run_.specified_vectors(a,b,c);
    if(type=="general"){
      run_.forward_solver_general();
    }
    if(type=="specialized"){
      run_.forward_solver_specialized();
    }
    run_.backward_solver();
    //std::cout << "done" << std::endl;
    run_.error();

    double epsilon_max = run_.ret_eps_max();
    double h = run_.ret_h();
    long double runtime = run_.ret_runtime();

    outfile << "N=1E" << i<< " epsilon_max=" << epsilon_max<<" log_10(h)="<<log10(h)<<" CPU-time=" << runtime << std::endl;

    }
    outfile.close();
}



int main() {
  //Setting up variables
  int n[3] = {10,100,1000};

  int a = -1;
  int b = 2;
  int c = -1;

std::string filename;

    for(int i=0; i< 3;i++){
      matrix_solver run_(n[i]);
      run_.specified_vectors(a,b,c);
      run_.forward_solver_general();
      run_.backward_solver();
      run_.error();

      //Making a filename
      filename = "./data/genN" + std::to_string(n[i]) + ".txt";
      run_.write_file(filename);
    }

    for(int i=0; i< 3;i++){
      matrix_solver run_(n[i]);
      run_.specified_vectors(a,b,c);
      run_.forward_solver_specialized();
      run_.backward_solver();
      run_.error();

      //Making a filename
      filename = "./data/speN" + std::to_string(n[i]) + ".txt";
      run_.write_file(filename);
    }

  filename = "./data/spe_stats.txt";
  write_stats(filename,1,7,"specialized", a, b, c);

  filename = "./data/gen_stats.txt";
  write_stats(filename,1,7,"general", a, b, c);

  delete [] a; delete [] b; delete [] c;
}
