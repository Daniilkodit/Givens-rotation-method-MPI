#ifndef HEADER
#define HEADER
#include "mpi.h"
#include <fenv.h>
#include<cstdio>
#include<fstream>
#include<string>
#include<cmath>
#include <string>
#include <cstring>
#include<time.h>
#include <unistd.h>
#include <sched.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
int feenableexcept(int excepts);
int fedisableexcept(int excepts);
int fegetexcept(void);

void Ainit(double* a, int n, int m, int p, int k, int s);

void print_matrix(double* A, int n, int m, int p, int k, double* buf, int max_print, MPI_Comm com);
int print_array(double* A, int n, int m, int printed_rows, int max_print);

int read_matrix(double* A, int n, int m, int p, int k, const std::string& name, double* buf, MPI_Comm com);
int read_array(std::ifstream& in, double* buf, int n);
double f(int i,int j,int n,int s);

int my_main(int argc, char*argv[]);

void Block_mult_add(double* Result, double* Block_A, double* Block_B, int row_A, int col_B, int col_row, const int m);

int Rotation (double *A, double *X, int n, int m, int id_process, int p, double eps,
              double *block, double *blockE, double *blockC,
              double *blockCE, double *rcos, double *rsin,
              double norm_mat, MPI_Comm com, double *buf1, double *buf2);

void get_block(double *, double *, int, int, int, int);
void put_block(double *, double *, int, int, int, int);

void Xmul_rotate (double *blockE, int m,int m1, int j, double *rcos,
                  double *rsin, int l, int k, double w, int stride);
void Xmul_rotate2 (double *blockE, double *blockCE,
                   int m1, int m2, int j, double *rcos, double *rsin, int l, int k, double w, int stride);
void mul_rotate (double *block, double *blockE, int m, int j, double *rcos,
                 double *rsin, int l, int k, double w, int stride);
int mul_fill_rotate2 (double *block, double *blockE, double *blockC,
                      double *blockCE, int m1, int m2, double *rcos, double *rsin,
                      double w);
void mul_rotate2 (double *block, double *blockE,
                  double *blockC, double *blockCE,
                  int m1, int m2, int j, double *rcos, double *rsin, int l, int k, double w, int stride);
int mul_fill_rotate (double *block, double *blockE, int m,int m2, double *rcos,
                     double *rsin, double w);
int
Inverse(double* X, double* A, int n, int l, double Anorm, double eps);

void Disr_block(double*Res,double*rmax,int m,int m1,int m2,int);
double Discrepancy(double*A, double*X, int n,int m,int k,int p, double* block,
                    double*blockE,double* Res,double* rmax, MPI_Comm com, double* buf);
double  Norm_matrix(double* A, int n, int m, int p, int id_process,  double* rmax, MPI_Comm com, double* rmax_dop );

void Block_multiplication(double* Result, double* Block_A, double* Block_B,int row_A,int col_B,int col_row,const int m);
void Block_mult_diff(double* Result, double* Block_A, double* Block_B, int row_A, int col_B,int col_row,const int m);

int l_to_g(int /* n */, int m, int p, int k, int i_loc);
int l_to_g_block(int /* n */, int /* m */, int p, int k, int i_loc_m);
int g_to_l(int /* n */, int m, int p, int k, int i_glob);
int g_to_l_block(int /* n */, int /* m */, int p, int k, int i_glob_m);

int get_max_block_rows(int n, int m, int p);
int get_k(int /*n*/, int m, int p, int i_glob);
int get_k_block(int /*n*/, int /*m*/, int p, int i_glob_m);
int get_block_rows(int n, int m, int p, int k);
int get_rows(int n, int m, int p, int k);

double get_full_time();
double get_cpu_time();
#endif