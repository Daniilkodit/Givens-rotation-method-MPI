#include "header.h"
double
  Discrepancy(double*A, double*X, int n,int m,int k,int p, double* block,
        double*blockE,double* res,double* rmax, MPI_Comm com, double* buf)
{
  if (n > 11000)
    return 0.0;

  int b_rows = get_block_rows(n, m, p, k);
  int max_block_rows = get_max_block_rows(n, m, p);
  int max_b = (n + m - 1) / m;
  int dst = (k - 1 + p) % p;
  int src = (k + 1 + p) % p;
  memset (rmax, 0, sizeof(double) * n);
  for (int s = 0; s < p; s++)
    {
      int sk = (k + s) % p;
      int sk_b_rows = get_block_rows (n, m, p, sk);
      // Сit = Ai,sk_i *Bsk_i,t  первый цикл по t так как хочу положить в res
      // сумму первого столбца произведения двух матриц Сit рамера h * v,
      // Ai,sk_i h*w, Bsk_i,t w*v
      for (int t = 0; t < max_b; t++)
        {
          memset (res, 0, m * m * sizeof(double));
          int t_glob_b = t; // у столбцов локальная и глобальная нумерациии совпадают
          int v = (t_glob_b * m + m < n ? m : n - m * t_glob_b); // нууууу должно работать, матрица же квадратная
          int max_h = 0;
          for (int i = 0; i < b_rows; i++)
            {
              int i_glob_b = l_to_g_block (n, m, p, k, i);
              int h = (i_glob_b * m + m < n ? m : n - m * i_glob_b);
              max_h = std::max (max_h, h);
              for (int sk_i = 0; sk_i < sk_b_rows; sk_i++) // по чужим
                {
                  int i_sk_glob = l_to_g_block (n, m, p, sk, sk_i);
                  int w = i_sk_glob * m + m < n ? m : n - i_sk_glob * m;
                  get_block (A, block, n, m, i, i_sk_glob);
                  get_block (X, blockE, n, m, sk_i, t);
                  Block_mult_add (res, block, blockE, h, v, w, m); // int row_A, int col_B, int col_row, const int m
                }
            }
          Disr_block (res, rmax, m, v, max_h, t * m); // h*v размер res именно такой порядок, потому что неадекватно написана функция
        }
      MPI_Status st;
      MPI_Sendrecv_replace (X, max_block_rows * m * n, MPI_DOUBLE, dst, 0, src, 0, com, &st);
    }
  // обмен всех процессов rmax, хочу чтобы у 0 го процесса было rmax
  // окончательное

  for (int ii = 0; ii < n; ii++)
    buf[ii] = rmax[ii];
  for (int s = 1; s < p; s++)
    {
      MPI_Status st;
      MPI_Sendrecv_replace (buf, n, MPI_DOUBLE, dst, 0, src, 0, com, &st);
      for (int ii = 0; ii < n; ii++)
        rmax[ii] += buf[ii];
    }
  double x = -1;
  for (int i = 0; i < n; i++)
    {
      if (x < rmax[i] - 1) // -1 т.к. надо единичную матрицу вычесть
        x = fabs(rmax[i] - 1);
    }
  return x;
}

double
  Norm_matrix (double* A, int n, int m, int p, int id_process, double* rmax, MPI_Comm com, double* rmax_dop )
{
  int col, row, rows = get_max_block_rows (n, m, p) * m;
  int dst = (id_process - 1 + p) % p;
  int src = (id_process + 1 + p) % p;

  for (col = 0; col < n; col++)
    {
      for (row = 0; row < rows; row++)
        {
          rmax[col] += fabs (A[row * n + col]);
        }
    }

  for (int ii = 0; ii < n; ii++)
      rmax_dop[ii] = rmax[ii];

   for (int s = 1; s < p; s++)
     {
       MPI_Status st;
       MPI_Sendrecv_replace(rmax_dop, n, MPI_DOUBLE, dst, 0, src, 0, com, &st);
       for (int ii = 0; ii < n; ii++)
         rmax[ii] += rmax_dop[ii];
     }
   double x = -1;
   for (int i = 0; i < n; i++)
     {
       if (x < rmax[i])
         x = rmax[i];
     }
   return x;
}

void
  Disr_block (double* Res, double* rmax, int m, int m1, int m2, int t)
{
  int row, col;
  double r1 = 0;
  for (col = 0; col < m1; col++)
    {
      for (row = 0; row < m2; row++)
        {
          r1 += Res[row * m + col];
        }
      rmax[t + col] += r1;
      r1 = 0;
    }
}
