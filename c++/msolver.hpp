class matrix_solver {

  private:
    //Setting up integers
    int n;
    double h;
    int a, b, c;
    //Setting up vectors
    arma::vec a_vec, b_vec, c_vec, g_vec;
    arma::vec b_tilde, c_tilde, g_tilde;
    arma::vec v_vec, x_vec, u_vec;


    //Setting up matrix
    arma::mat A, P, L, U;


    //error
    double epsilon, epsilon_max;

    long double runtime;

    std::string filename;

    clock_t start, finish;

    void Initialize(int n_){
      //Initialing stepsize, numbers ect.
      n = n_;
      h = 1./(n+2);

      //Initializing x, u(x), g(x)
      x_vec = arma::linspace(0,1,n+2);
      v_vec = arma::vec(n+2,arma::fill::zeros);
      u_vec = arma::vec(n+2,arma::fill::zeros);
      
      g_vec = arma::vec(n,arma::fill::zeros);

      //filling in values for g(x)
      for(int i=0;i<n;i++){
        g_vec(i) = pow(h,2) * f(x_vec(i+1));
      }

      //filling in analytic values for u(x)
      for(int i=0;i<n+2;i++){
        u_vec(i) = analytic(x_vec(i));
      }


      //setting up the tilde vectors
      b_tilde = arma::vec(n,arma::fill::zeros);
      c_tilde = arma::vec(n,arma::fill::zeros);
      g_tilde = arma::vec(n,arma::fill::zeros);


    }

    //overload function
    void Initialize(){
      Initialize(100);
    }


public:
  //setting up functions
  void random_vectors();
  void specified_vectors(int a_,int b_,int c_);

  double f(double x);
  double analytic(double x);

  void error();

  double ret_eps_max();
  double ret_h();
  double ret_runtime();

  void forward_solver_general();
  void forward_solver_specialized();
  void backward_solver();

  void write_file(std::string filename_);

  void write_stats(std::string filename_);

  void LU_decomp(int a_, int b_, int c_, std::string filename_);

  //setting up the overload
  matrix_solver(int n_){
    Initialize(n_);
  }
  matrix_solver(){
    Initialize();
  }
};
