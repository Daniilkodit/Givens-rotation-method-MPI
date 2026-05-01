#include "header.h"

int
Rotation (double *A, double *X, int n, int m, int id_process, int p, double eps,
          double *block, double *blockE, double *blockC,
          double *blockCE, double *rcos, double *rsin,
          double norm_mat,  MPI_Comm com, double *buf1, double *buf2)
{
  int k = n / m, l = n - k * m, k_l = (l != 0 ? k + 1 : k);
  int r, j_loc, i_loc, c_loc, t_loc;
  int i_glob;
  double Anorm = norm_mat;
  double w = eps * Anorm;
  if (Anorm <= 1.e-64)
    {
      return -1;
    }

  for (int r = 0; r < k; r++)
    {
      int start = 0;

      if (id_process < (r % p))
        {
          start = g_to_l_block (n,m,p,id_process, r - (r % p) + p + id_process);
        }
      else
        {
          start = g_to_l_block (n,m,p,id_process,r);
        }

      int start_glob = l_to_g_block (n,m,p,id_process,start);

      if (start_glob < k)
        {
          j_loc = r;
          get_block (A, block, n, m, start, j_loc);
          get_block (X, blockE, n, m, start, j_loc);
          mul_fill_rotate (block, blockE, m, m, rcos, rsin, w);
          put_block (A, block, n, m, start, j_loc);
          put_block (X, blockE, n, m, start, j_loc);

          for (c_loc = 0; c_loc < r; c_loc++)
            {
              Xmul_rotate (X + start * m * n + c_loc * m, m, m, c_loc, rcos, rsin, l, k, eps, n);
            }

          for (j_loc = r + 1; j_loc < k_l; j_loc++)
            {
              mul_rotate (A + start * m * n + j_loc * m, X + start * m * n + j_loc * m, m, j_loc, rcos, rsin, l, k, eps, n);
            }
        }

      int K = get_block_rows (n, m, p, id_process);

      for (i_loc = start + 1; i_loc < K; i_loc++)
        {
          j_loc = r;
          int i_glob = l_to_g_block(n, m, p, id_process, i_loc);
          int v = i_glob * m + m < n ? m : n - m * i_glob;

          get_block (A, block, n, m, start, j_loc);
          get_block (A, blockC, n, m, i_loc, j_loc);
          get_block (X, blockE, n, m, start, j_loc);
          get_block (X, blockCE, n, m, i_loc, j_loc);
          mul_fill_rotate2 (block, blockE, blockC, blockCE, m, v, rcos, rsin, w);
          put_block (A, block, n, m, start, j_loc);
          put_block (A, blockC, n, m, i_loc, j_loc);
          put_block (X, blockE, n, m, start, j_loc);
          put_block (X, blockCE, n, m, i_loc, j_loc);

          for (c_loc = 0; c_loc < r; c_loc++)
            {
              Xmul_rotate2 (X + start * m * n + c_loc * m, X + i_loc * m * n + c_loc * m, m, v, c_loc, rcos, rsin, l, k, eps, n);
            }

          for (j_loc = r + 1; j_loc < k_l; j_loc++)
            {
              mul_rotate2 (A + start * m * n + j_loc * m, X + start * m * n + j_loc * m,
                           A + i_loc * m * n + j_loc * m, X + i_loc * m * n + j_loc * m,
                           m, v, j_loc, rcos, rsin, l, k, eps, n);
            }
        }

      for(int ii = 1; ii < p && ii < k - r + 1; ii *= 2)
        {
          int start_glob = l_to_g_block (n, m, p, id_process, start);
          int id_new = start_glob - r;
          bool main = false, assist = false;
          int src = -1, dst = -1;

          if ((id_new % (ii * 2) == 0 || id_new % (ii * 2) == ii) && start_glob < k_l)
            {
              if (id_new % (ii * 2) == 0)
                {
                  main = true;
                  src = (start_glob + ii) % p;
                  dst = src;
                }
              else
                {
                  assist = true;
                  dst = (start_glob - ii) % p;
                  src = dst;
                }
            }

          if (main)
            {
              j_loc = r;

              if (start_glob + ii < k_l && start_glob + ii <= r + p - 1)
                {
                  int v_dst = (start_glob + ii) * m + m < n ? m : n - m * (start_glob + ii);
                  int v_src = m;

                  MPI_Status st;

                  MPI_Sendrecv (A + start * m * n, n * v_src, MPI_DOUBLE, dst, 0,
                                buf1,              n * v_dst, MPI_DOUBLE, src, 0,
                                com, &st);

                  MPI_Sendrecv (X + start * m * n, n * v_src, MPI_DOUBLE, dst, 0,
                                buf2,              n * v_dst, MPI_DOUBLE, src, 0,
                                com, &st);

                  get_block (A, block, n, m, start, r);
                  get_block (buf1, blockC, n, m, 0, r);
                  get_block (X, blockE, n, m, start, r);
                  get_block (buf2, blockCE, n, m, 0, r);
                  mul_fill_rotate2 (block, blockE, blockC, blockCE, m, v_dst, rcos, rsin, w);
                  put_block (A, block, n, m, start, r);
                  put_block (buf1, blockC, n, m, 0, r);
                  put_block (X, blockE, n, m, start, r);
                  put_block (buf2, blockCE, n, m, 0, r);

                  int a_mid = (r + k_l) / 2;
                  int x_mid = k_l / 2;

                  for (c_loc = x_mid; c_loc < r; c_loc++)
                    {
                      Xmul_rotate2 (X + start * m * n + c_loc * m, buf2 + c_loc * m, m, v_dst, c_loc, rcos, rsin, l, k, eps, n);
                    }

                  for (c_loc = std::max (x_mid, r + 1); c_loc < a_mid; c_loc++)
                    {
                      Xmul_rotate2 (X + start * m * n + c_loc * m, buf2 + c_loc * m, m, v_dst, c_loc, rcos, rsin, l, k, eps, n);
                    }

                  for (j_loc = std::max(a_mid, r + 1); j_loc < k_l; j_loc++)
                    {
                      mul_rotate2 (A + start * m * n + j_loc * m, X + start * m * n + j_loc * m,
                                   buf1 + j_loc * m, buf2 + j_loc * m,
                                   m, v_dst, j_loc, rcos, rsin, l, k, eps, n);
                    }

                  int a_right_len = n - a_mid * m;
                  int a_left_len  = a_mid * m - r * m;
                  int x_right_len = n - x_mid * m;
                  int x_left_len  = x_mid * m;

                  MPI_Datatype a_send_type, a_recv_type, x_send_type, x_recv_type;
                  MPI_Type_vector(v_dst, a_right_len, n, MPI_DOUBLE, &a_send_type);
                  MPI_Type_vector(v_src, a_left_len,  n, MPI_DOUBLE, &a_recv_type);
                  MPI_Type_vector(v_dst, x_right_len, n, MPI_DOUBLE, &x_send_type);
                  MPI_Type_vector(v_src, x_left_len,  n, MPI_DOUBLE, &x_recv_type);
                  MPI_Type_commit(&a_send_type);
                  MPI_Type_commit(&a_recv_type);
                  MPI_Type_commit(&x_send_type);
                  MPI_Type_commit(&x_recv_type);

                  MPI_Sendrecv(buf1 + a_mid * m,           1, a_send_type, dst, 2,
                                A + start * m * n + r * m,  1, a_recv_type, src, 2,
                                com, &st);

                  MPI_Sendrecv(buf2 + x_mid * m,           1, x_send_type, dst, 3,
                                X + start * m * n,          1, x_recv_type, src, 3,
                                com, &st);

                  MPI_Type_free(&a_send_type);
                  MPI_Type_free(&a_recv_type);
                  MPI_Type_free(&x_send_type);
                  MPI_Type_free(&x_recv_type);
                }
            }

          if (assist)
            {
              j_loc = r;
              if (start_glob < k_l)
                {
                  int v_src = start_glob * m + m < n ? m : n - m * start_glob;
                  int v_dst = m;

                  MPI_Status st;

                  MPI_Sendrecv (A + start * m * n, n * v_src, MPI_DOUBLE, dst, 0,
                                buf1,              n * v_dst, MPI_DOUBLE, src, 0,
                                com, &st);

                  MPI_Sendrecv (X + start * m * n, n * v_src, MPI_DOUBLE, dst, 0,
                                buf2,              n * v_dst, MPI_DOUBLE, src, 0,
                                com, &st);

                  get_block (buf1, block, n, m, 0, r);
                  get_block (A, blockC, n, m, start, r);
                  get_block (buf2, blockE, n, m, 0, r);
                  get_block (X, blockCE, n, m, start, r);
                  mul_fill_rotate2 (block, blockE, blockC, blockCE, m, v_src, rcos, rsin, w);
                  put_block (buf1, block, n, m, 0, r);
                  put_block (A, blockC, n, m, start, r);
                  put_block (buf2, blockE, n, m, 0, r);
                  put_block (X, blockCE, n, m, start, r);

                  int a_mid = (r + k_l) / 2;
                  int x_mid = k_l / 2;

                  for (c_loc = 0; c_loc < std::min(r, x_mid); c_loc++)
                    {
                      Xmul_rotate2 (buf2 + c_loc * m, X + start * m * n + c_loc * m, m, v_src, c_loc, rcos, rsin, l, k, eps, n);
                    }

                  for (j_loc = r + 1; j_loc < std::min(a_mid, k_l); j_loc++)
                    {
                      mul_rotate2 (buf1 + j_loc * m, buf2 + j_loc * m,
                                   A + start * m * n + j_loc * m, X + start * m * n + j_loc * m,
                                   m, v_src, j_loc, rcos, rsin, l, k, eps, n);
                    }

                  int a_right_len = n - a_mid * m;
                  int a_left_len  = a_mid * m - r * m;
                  int x_right_len = n - x_mid * m;
                  int x_left_len  = x_mid * m;

                  MPI_Datatype a_send_type, a_recv_type, x_send_type, x_recv_type;
                  MPI_Type_vector(v_dst, a_left_len,  n, MPI_DOUBLE, &a_send_type);
                  MPI_Type_vector(v_src, a_right_len, n, MPI_DOUBLE, &a_recv_type);
                  MPI_Type_vector(v_dst, x_left_len,  n, MPI_DOUBLE, &x_send_type);
                  MPI_Type_vector(v_src, x_right_len, n, MPI_DOUBLE, &x_recv_type);
                  MPI_Type_commit(&a_send_type);
                  MPI_Type_commit(&a_recv_type);
                  MPI_Type_commit(&x_send_type);
                  MPI_Type_commit(&x_recv_type);

                  MPI_Sendrecv(buf1 + r * m,                  1, a_send_type, dst, 2,
                                A + start * m * n + a_mid * m, 1, a_recv_type, src, 2,
                                com, &st);

                  MPI_Sendrecv(buf2,                          1, x_send_type, dst, 3,
                                X + start * m * n + x_mid * m, 1, x_recv_type, src, 3,
                                com, &st);

                  MPI_Type_free(&a_send_type);
                  MPI_Type_free(&a_recv_type);
                  MPI_Type_free(&x_send_type);
                  MPI_Type_free(&x_recv_type);
                }
            }
        }
    }

  if (l != 0
      && (k % p) == id_process)
    {
      r = g_to_l_block (n, m, p, id_process, k);
      j_loc = k;

      get_block (A, block, n, m, r, j_loc);
      get_block (X, blockE, n, m, r, j_loc);
      mul_fill_rotate (block, blockE, m, l, rcos, rsin, w);
      put_block (A, block, n, m, r, j_loc);
      put_block (X, blockE, n, m, r, j_loc);

      for (c_loc = 0; c_loc < k; c_loc++)
        {
          Xmul_rotate (X + r * m * n + c_loc * m, m, l, c_loc, rcos, rsin, l, k, eps, n);
        }
    }

  int err = 0, res = 0;
  int K = get_block_rows (n, m, p, id_process);
  int main_k = 0;
  for (i_loc = 0; i_loc < K; i_loc++)
    {
      i_glob = l_to_g_block(n, m, p, id_process,i_loc);
      int max_b = (n + m - 1) / m;

      get_block (A, block, n, m, i_loc, i_glob );
      int v = i_glob * m + m < n ? m : n - m * i_glob;

      if (Inverse (blockCE, block, m, v, Anorm, eps) != 0)
        {
          err = 1;
          MPI_Allreduce (&err, &res, 1, MPI_INT, MPI_SUM, com);
          return res;
        }

      for (j_loc = i_glob  + 1; j_loc < max_b; j_loc++)
        {
          int j_glob = j_loc;
          int h = j_glob * m + m < n ? m : n - m * j_glob;

          get_block (A, blockE, n, m, i_loc, j_loc);
          Block_multiplication (block, blockCE, blockE, v, h, v, m);
          put_block (A, block, n, m, i_loc, j_loc);
        }

      for (t_loc = 0; t_loc < max_b; t_loc++)
        {
          int t_glob = t_loc;
          int h = t_glob * m + m < n ? m : n - m * t_glob;

          get_block (X, blockE, n, m, i_loc, t_loc);
          Block_multiplication (block, blockCE, blockE, v, h, v, m);
          put_block (X, block, n, m, i_loc, t_loc);
        }
    }

  MPI_Allreduce (&err, &res, 1, MPI_INT, MPI_SUM, com);
  if (res)
    return res;

  for (int step = k_l - 2; step >= 0; step--)
    {
      main_k = get_k_block (n, m, p, step);
      int last_str = get_block_rows(n, m, p, id_process) - 1;
      int last_glob_str = l_to_g_block(n, m, p, id_process, last_str);
      int max_b = (n + m - 1) / m;

      if (main_k == id_process)
        {
          i_loc = g_to_l_block (n, m, p, id_process, step);
          MPI_Bcast (A + i_loc * m * n, n * m, MPI_DOUBLE, main_k, com);

          for (j_loc = 0; j_loc < max_b; j_loc++)
            {
              get_block (X, blockE, n, m, i_loc, j_loc);
              memset (blockC, 0, m * m * sizeof(double));

              int j_glob = j_loc;
              int w_val = (j_glob * m + m < n ? m : n - m * j_glob);

              for(t_loc = i_loc + 1; t_loc < K; t_loc++)
                {
                  int t_glob = l_to_g_block(n, m, p, id_process,t_loc);
                  int h = t_glob * m + m < n ? m : n - t_glob * m;

                  get_block (A, block, n, m, i_loc, t_glob);
                  get_block (X, blockCE, n, m, t_loc, j_loc);
                  Block_mult_diff (blockE, block, blockCE, m, w_val, h, m);
                }

              MPI_Reduce (blockE, blockC, m * m, MPI_DOUBLE, MPI_SUM, main_k, com);
              put_block(X, blockC, n, m, i_loc, j_loc);
            }
        }
      else if (last_glob_str > step)
        {
          MPI_Bcast (buf1, n * m, MPI_DOUBLE, main_k, com);

          while (last_str > 0 && last_glob_str  > step)
            {
              last_str--;
              last_glob_str = l_to_g_block(n, m, p, id_process, last_str);
            }

          if (last_glob_str  <= step)
            {
              last_str++;
              last_glob_str = l_to_g_block(n, m, p, id_process, last_str);
            }

          for (j_loc = 0; j_loc < max_b; j_loc++)
            {
              memset(blockE, 0, sizeof(double) * m * m);

              int j_glob = j_loc;
              int w_val = (j_glob * m + m < n ? m : n - m * j_glob);

              for(t_loc = last_str; t_loc < K; t_loc++)
                {
                  int t_glob = l_to_g_block(n, m, p, id_process,t_loc);
                  int h = t_glob * m + m < n ? m : n - t_glob * m;

                  get_block (buf1, block, n, m, 0, t_glob);
                  get_block (X, blockCE, n, m, t_loc, j_loc);
                  Block_mult_diff (blockE, block, blockCE, m, w_val, h, m);
                }
              MPI_Reduce (blockE, block, m * m, MPI_DOUBLE, MPI_SUM, main_k, com);
            }
        }
      else
        {
          MPI_Bcast (buf1, n * m, MPI_DOUBLE, main_k, com);
          memset(blockE, 0, sizeof(double) * m * m);

          for (j_loc = 0; j_loc < max_b; j_loc++)
            {
              MPI_Reduce (blockE, block, m * m, MPI_DOUBLE, MPI_SUM, main_k, com);
            }
        }
    }

  return 0;
}

void
get_block(double* matr, double* block, int real_n, int block_n, int i, int j)
{
  int k = real_n / block_n;
  int l = real_n - k * block_n;
  int width = (j < k ? block_n : l);
  int height = (i < k ? block_n : l);
  int row, col;

  double* source_block = matr + i * real_n * block_n + j * block_n;

  for (row = 0; row < height; row++)
    {
      for (col = 0; col < width; col++)
        {
          block[row * block_n + col] = source_block[row * real_n + col];
        }
    }
}

int
Inverse(double* X, double* A, int n, int l, double Anorm, double eps)
{
  int i, j, t;
  for (i = 0; i < l; i++)
    {
      for (j = 0; j < l; j++)
        {
          X[i * n + j] = (i == j) ? 1.0 : 0.0;
        }
    }

  for (t = 0; t < l; t++)
    {
      for (i = l - 1; i >= 0; i--)
        {
          double sum = 0.0;

          if (fabs(A[i * n + i]) < eps * Anorm)
            return -1;

          j = i + 1;
          int remainder = (l - j) % 4;

          for (int k = 0; k < remainder; k++)
            {
              sum += A[i * n + j] * X[j * n + t];
              j++;
            }

          for (; j < l; j += 4)
            {
              sum += A[i * n + j] * X[j * n + t];
              sum += A[i * n + j + 1] * X[(j + 1) * n + t];
              sum += A[i * n + j + 2] * X[(j + 2) * n + t];
              sum += A[i * n + j + 3] * X[(j + 3) * n + t];
            }

          X[i * n + t] = (X[i * n + t] - sum) / A[i * n + i];
        }
    }
  return 0;
}

void
Xmul_rotate (double *blockE, int m,int m1, int j, double *rcos,
             double *rsin, int l, int k, double w, int stride)
{
  double cosk, sink, x0, y0, x1, y1;
  int p, q, c;
  int z = (j == k) ? l : m;
  if(j==k && l==0) return;
  for (q = 0; q < m1 - 1; q++)
    {
      for (p = q + 1; p < m1; p++)
        {
          cosk = rcos[m * q + p];
          sink = rsin[m * q + p];
          if (fabs (sink) < w && fabs (cosk - 1) < w)
            continue;
          double *rq = blockE + q * stride;
          double *rp = blockE + p * stride;
          int rem = z % 2;
          for (c = 0; c < rem; c++)
            {
              x0 = cosk * rq[c] - sink * rp[c];
              y0 = sink * rq[c] + cosk * rp[c];
              rq[c] = x0;
              rp[c] = y0;
            }
          for (c = rem; c < z; c += 2)
            {
              x0 = cosk * rq[c]   - sink * rp[c];
              y0 = sink * rq[c]   + cosk * rp[c];
              x1 = cosk * rq[c+1] - sink * rp[c+1];
              y1 = sink * rq[c+1] + cosk * rp[c+1];
              rq[c]   = x0; rp[c]   = y0;
              rq[c+1] = x1; rp[c+1] = y1;
            }
        }
    }
}

void
Xmul_rotate2 (double *blockE, double *blockCE,
              int m1, int m2, int j, double *rcos, double *rsin, int l, int k, double w, int stride)
{
  double cosk, sink, x0, y0, x1, y1;
  int p, q, c;
  int z = (j == k) ? l : m1;
  if(j==k && l==0) return;
  for (q = 0; q < m1; q++)
    {
      for (p = 0; p < m2; p++)
        {
          cosk = rcos[m1 * q + p];
          sink = rsin[m1 * q + p];
          if (fabs (sink) < w && fabs (cosk - 1) < w)
            continue;
          double *rq = blockE  + q * stride;
          double *rp = blockCE + p * stride;
          int rem = z % 2;
          for (c = 0; c < rem; c++)
            {
              x0 = cosk * rq[c] - sink * rp[c];
              y0 = sink * rq[c] + cosk * rp[c];
              rq[c] = x0;
              rp[c] = y0;
            }
          for (c = rem; c < z; c += 2)
            {
              x0 = cosk * rq[c]   - sink * rp[c];
              y0 = sink * rq[c]   + cosk * rp[c];
              x1 = cosk * rq[c+1] - sink * rp[c+1];
              y1 = sink * rq[c+1] + cosk * rp[c+1];
              rq[c]   = x0; rp[c]   = y0;
              rq[c+1] = x1; rp[c+1] = y1;
            }
        }
    }
}

void mul_rotate (double *block, double *blockE, int m, int j, double *rcos,
            double *rsin, int l, int k, double w, int stride)
{
  double cosk, sink, x0, y0, x1, y1;
  double xe0, ye0, xe1, ye1;
  int p, q, c;
  int z = (j == k) ? l : m;

  if(j==k &&l==0)return;
  for (q = 0; q < m - 1; q++)
    {
      for (p = q + 1; p < m; p++)
        {
          cosk = rcos[m * q + p];
          sink = rsin[m * q + p];
          if (fabs (sink) < w && fabs (cosk - 1) < w)
            continue;
          double *bq = block  + q * stride;
          double *bp = block  + p * stride;
          double *eq = blockE + q * stride;
          double *ep = blockE + p * stride;
          int rem = z % 2;
          for (c = 0; c < rem; c++)
            {
              x0 = cosk * bq[c] - sink * bp[c];
              y0 = sink * bq[c] + cosk * bp[c];
              bq[c] = x0;
              bp[c] = y0;
              xe0 = cosk * eq[c] - sink * ep[c];
              ye0 = sink * eq[c] + cosk * ep[c];
              eq[c] = xe0;
              ep[c] = ye0;
            }
          for (c = rem; c < z; c += 2)
            {
              x0 = cosk * bq[c]   - sink * bp[c];
              y0 = sink * bq[c]   + cosk * bp[c];
              x1 = cosk * bq[c+1] - sink * bp[c+1];
              y1 = sink * bq[c+1] + cosk * bp[c+1];
              bq[c]   = x0; bp[c]   = y0;
              bq[c+1] = x1; bp[c+1] = y1;

              xe0 = cosk * eq[c]   - sink * ep[c];
              ye0 = sink * eq[c]   + cosk * ep[c];
              xe1 = cosk * eq[c+1] - sink * ep[c+1];
              ye1 = sink * eq[c+1] + cosk * ep[c+1];
              eq[c]   = xe0; ep[c]   = ye0;
              eq[c+1] = xe1; ep[c+1] = ye1;
            }

        }
    }
}

int mul_fill_rotate2 (double *block, double *blockE, double *blockC,
                  double *blockCE, int m1, int m2, double *rcos, double *rsin,
                  double w)
{
  double cosk, sq, sink, x0, y0, x1, y1;
  double xe0, ye0, xe1, ye1;
  int p, q, c;
  for (q = 0; q < m1; q++)
    {
      for (p = 0; p < m2; p++)
        {
          x0 = block[q * m1 + q];
          y0 = blockC[p * m1 + q];
          sq = (sqrt (x0 * x0 + y0 * y0));
          if (sq < w)
            {
              rcos[m1 * q + p] = 1;
              rsin[m1 * q + p] = 0;
              block[q * m1 + q] = 0;
              blockC[p * m1 + q] = 0;
              continue;
            }
          rcos[m1 * q + p] = x0 / sq;
          rsin[m1 * q + p] = -y0 / sq;
          cosk = rcos[m1 * q + p];
          sink = rsin[m1 * q + p];

          double *bq  = block   + q * m1;
          double *bp  = blockC  + p * m1;
          double *eq  = blockE  + q * m1;
          double *ep  = blockCE + p * m1;
          int rem = m1 % 2;
          for (c = 0; c < rem; c++)
            {
              x0 = cosk * bq[c] - sink * bp[c];
              y0 = sink * bq[c] + cosk * bp[c];
              bq[c] = x0;
              bp[c] = y0;
              xe0 = cosk * eq[c] - sink * ep[c];
              ye0 = sink * eq[c] + cosk * ep[c];
              eq[c] = xe0;
              ep[c] = ye0;
            }
          for (c = rem; c < m1; c += 2)
            {
              x0 = cosk * bq[c]   - sink * bp[c];
              y0 = sink * bq[c]   + cosk * bp[c];
              x1 = cosk * bq[c+1] - sink * bp[c+1];
              y1 = sink * bq[c+1] + cosk * bp[c+1];
              bq[c]   = x0; bp[c]   = y0;
              bq[c+1] = x1; bp[c+1] = y1;

              xe0 = cosk * eq[c]   - sink * ep[c];
              ye0 = sink * eq[c]   + cosk * ep[c];
              xe1 = cosk * eq[c+1] - sink * ep[c+1];
              ye1 = sink * eq[c+1] + cosk * ep[c+1];
              eq[c]   = xe0; ep[c]   = ye0;
              eq[c+1] = xe1; ep[c+1] = ye1;
            }
          blockC[p * m1 + q] = 0.0;

        }
    }
  return 0;

}

void mul_rotate2 (double *block, double *blockE,
             double *blockC, double *blockCE,
             int m1, int m2, int j, double *rcos, double *rsin, int l, int k, double w, int stride)
{
  double cosk, sink, x0, y0, x1, y1;
  double xe0, ye0, xe1, ye1;
  int p, q, c;
  int z = (j == k) ? l : m1;
  if(j==k && l==0) return;
  for (q = 0; q < m1; q++)
    {
      for (p = 0; p < m2; p++)
        {
          cosk = rcos[m1 * q + p];
          sink = rsin[m1 * q + p];
          if (fabs (sink) < w && fabs (cosk - 1) < w)
            continue;
          double *bq  = block   + q * stride;
          double *bp  = blockC  + p * stride;
          double *eq  = blockE  + q * stride;
          double *ep  = blockCE + p * stride;
          int rem = z % 2;
          for (c = 0; c < rem; c++)
            {
              x0 = cosk * bq[c] - sink * bp[c];
              y0 = sink * bq[c] + cosk * bp[c];
              bq[c] = x0;
              bp[c] = y0;
              xe0 = cosk * eq[c] - sink * ep[c];
              ye0 = sink * eq[c] + cosk * ep[c];
              eq[c] = xe0;
              ep[c] = ye0;
            }
          for (c = rem; c < z; c += 2)
            {
              x0 = cosk * bq[c]   - sink * bp[c];
              y0 = sink * bq[c]   + cosk * bp[c];
              x1 = cosk * bq[c+1] - sink * bp[c+1];
              y1 = sink * bq[c+1] + cosk * bp[c+1];
              bq[c]   = x0; bp[c]   = y0;
              bq[c+1] = x1; bp[c+1] = y1;

              xe0 = cosk * eq[c]   - sink * ep[c];
              ye0 = sink * eq[c]   + cosk * ep[c];
              xe1 = cosk * eq[c+1] - sink * ep[c+1];
              ye1 = sink * eq[c+1] + cosk * ep[c+1];
              eq[c]   = xe0; ep[c]   = ye0;
              eq[c+1] = xe1; ep[c+1] = ye1;
            }

        }
    }
}

int mul_fill_rotate (double *block, double *blockE, int m,int m2, double *rcos,
                 double *rsin, double w)
{
  int q, p, c;
  double cosk, sink, x0, y0, x1, y1, sq;
  double xe0, ye0, xe1, ye1;
  for (q = 0; q < m2 - 1; q++)
    {
      for (p = q + 1; p < m2; p++)
        {
          x0 = block[q * m + q];
          y0 = block[p * m + q];
          sq = (sqrt (x0 * x0 + y0 * y0));
          if (sq < w)
            {
              rcos[m * q + p] = 1;
              rsin[m * q + p] = 0;
              block[q * m + q] = 0;
              block[p * m + q] = 0;
              continue;
            }
          rcos[m * q + p] = x0 / sq;
          rsin[m * q + p] = -y0 / sq;
          cosk = rcos[m * q + p];
          sink = rsin[m * q + p];
          double *bq = block  + q * m;
          double *bp = block  + p * m;
          double *eq = blockE + q * m;
          double *ep = blockE + p * m;
          int rem = m2 % 2;
          for (c = 0; c < rem; c++)
            {
              x0 = cosk * bq[c] - sink * bp[c];
              y0 = sink * bq[c] + cosk * bp[c];
              bq[c] = x0;
              bp[c] = y0;
              xe0 = cosk * eq[c] - sink * ep[c];
              ye0 = sink * eq[c] + cosk * ep[c];
              eq[c] = xe0;
              ep[c] = ye0;
            }
          for (c = rem; c < m2; c += 2)
            {
              x0 = cosk * bq[c]   - sink * bp[c];
              y0 = sink * bq[c]   + cosk * bp[c];
              x1 = cosk * bq[c+1] - sink * bp[c+1];
              y1 = sink * bq[c+1] + cosk * bp[c+1];
              bq[c]   = x0; bp[c]   = y0;
              bq[c+1] = x1; bp[c+1] = y1;

              xe0 = cosk * eq[c]   - sink * ep[c];
              ye0 = sink * eq[c]   + cosk * ep[c];
              xe1 = cosk * eq[c+1] - sink * ep[c+1];
              ye1 = sink * eq[c+1] + cosk * ep[c+1];
              eq[c]   = xe0; ep[c]   = ye0;
              eq[c+1] = xe1; ep[c+1] = ye1;

            }
          block[p * m + q] = 0.0;

        }
    }
  return 0;
}


void
put_block(double* matr, double* block, int real_n, int block_n, int i, int j)
{
  int k = real_n / block_n;
  int l = real_n - k * block_n;
  int width = (j < k ? block_n : l);
  int height = (i < k ? block_n : l);
  int row, col;

  double* target_block = matr + i * real_n * block_n + j * block_n;

  for (row = 0; row < height; row++)
    {
      for (col = 0; col < width; col++)
        {
          target_block[row * real_n + col] = block[row * block_n + col];
        }
    }
}

void
Block_mult_diff(double* Result,
                 double* Block_A,
                 double* Block_B,
                 int row_A,
                 int col_B,
                 int col_row,
                 const int m)
{
  int row_l = row_A % 3;
  int col_l = col_B % 3;
  int row_k = (row_A - row_l) / 3;
  int col_k = (col_B - col_l) / 3;
  double res_00 = 0, res_01 = 0, res_02 = 0;
  double res_10 = 0, res_11 = 0, res_12 = 0;
  double res_20 = 0, res_21 = 0, res_22 = 0;
  int col_row1 = col_row;
  col_B = m;
  row_A = m;
  col_row = m;
  for (int b_i = 0; b_i < row_k; b_i++)
    {
      for (int b_j = 0; b_j < col_k; b_j++)
        {
          res_00 = 0, res_01 = 0, res_02 = 0;
          res_10 = 0, res_11 = 0, res_12 = 0;
          res_20 = 0, res_21 = 0, res_22 = 0;
          for (int s = 0; s < col_row1; s++)
            {
              res_00 += Block_A[b_i * 3 * col_row + s]
                        * Block_B[s * col_B + b_j * 3];
              res_01 += Block_A[b_i * 3 * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 1];
              res_02 += Block_A[b_i * 3 * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 2];

              res_10 += Block_A[(b_i * 3 + 1) * col_row + s]
                        * Block_B[s * col_B + b_j * 3];
              res_11 += Block_A[(b_i * 3 + 1) * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 1];
              res_12 += Block_A[(b_i * 3 + 1) * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 2];

              res_20 += Block_A[(b_i * 3 + 2) * col_row + s]
                        * Block_B[s * col_B + b_j * 3];
              res_21 += Block_A[(b_i * 3 + 2) * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 1];
              res_22 += Block_A[(b_i * 3 + 2) * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 2];
            }

          Result[b_i * 3 * col_B + b_j * 3] -= res_00;
          Result[b_i * 3 * col_B + b_j * 3 + 1] -= res_01;
          Result[b_i * 3 * col_B + b_j * 3 + 2] -= res_02;
          Result[(b_i * 3 + 1) * col_B + b_j * 3] -= res_10;
          Result[(b_i * 3 + 1) * col_B + b_j * 3 + 1] -= res_11;
          Result[(b_i * 3 + 1) * col_B + b_j * 3 + 2] -= res_12;
          Result[(b_i * 3 + 2) * col_B + b_j * 3] -= res_20;
          Result[(b_i * 3 + 2) * col_B + b_j * 3 + 1] -= res_21;
          Result[(b_i * 3 + 2) * col_B + b_j * 3 + 2] -= res_22;
        }

      if (col_l != 0)
        {
          res_00 = 0, res_01 = 0, res_10 = 0;
          res_11 = 0, res_20 = 0, res_21 = 0;

          for (int s = 0; s < col_row1; s++)
            {
              if (col_l > 1)
                {
                  res_01 += Block_A[b_i * 3 * col_row + s]
                            * Block_B[s * col_B + col_k * 3 + 1];
                  res_11 += Block_A[(b_i * 3 + 1) * col_row + s]
                            * Block_B[s * col_B + col_k * 3 + 1];
                  res_21 += Block_A[(b_i * 3 + 2) * col_row + s]
                            * Block_B[s * col_B + col_k * 3 + 1];
                }

              res_00 += Block_A[b_i * 3 * col_row + s]
                        * Block_B[s * col_B + col_k * 3];
              res_10 += Block_A[(b_i * 3 + 1) * col_row + s]
                        * Block_B[s * col_B + col_k * 3];
              res_20 += Block_A[(b_i * 3 + 2) * col_row + s]
                        * Block_B[s * col_B + col_k * 3];
            }

          Result[b_i * 3 * col_B + col_k * 3] -= res_00;
          Result[(b_i * 3 + 1) * col_B + col_k * 3] -= res_10;
          Result[(b_i * 3 + 2) * col_B + col_k * 3] -= res_20;

          if (col_l > 1)
            {
              Result[b_i * 3 * col_B + col_k * 3 + 1] -= res_01;
              Result[(b_i * 3 + 1) * col_B + col_k * 3 + 1] -= res_11;
              Result[(b_i * 3 + 2) * col_B + col_k * 3 + 1] -= res_21;
            }
        }
    }

  if (row_l != 0)
    {
      for (int b_j = 0; b_j < col_k; b_j++)
        {
          res_00 = 0, res_01 = 0, res_02 = 0;
          res_10 = 0, res_11 = 0, res_12 = 0;
          for (int s = 0; s < col_row1; s++)
            {
              res_00 += Block_A[(row_k * 3) * col_row + s]
                        * Block_B[s * col_B + b_j * 3];
              res_01 += Block_A[(row_k * 3) * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 1];
              res_02 += Block_A[(row_k * 3) * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 2];

              if (row_l > 1)
                {
                  res_10 += Block_A[(row_k * 3 + 1) * col_row + s]
                            * Block_B[s * col_B + b_j * 3];
                  res_11 += Block_A[(row_k * 3 + 1) * col_row + s]
                            * Block_B[s * col_B + b_j * 3 + 1];
                  res_12 += Block_A[(row_k * 3 + 1) * col_row + s]
                            * Block_B[s * col_B + b_j * 3 + 2];
                }
            }

          Result[row_k * 3 * col_B + b_j * 3] -= res_00;
          Result[row_k * 3 * col_B + b_j * 3 + 1] -= res_01;
          Result[row_k * 3 * col_B + b_j * 3 + 2] -= res_02;

          if (row_l > 1)
            {
              Result[(row_k * 3 + 1) * col_B + b_j * 3] -= res_10;
              Result[(row_k * 3 + 1) * col_B + b_j * 3 + 1] -= res_11;
              Result[(row_k * 3 + 1) * col_B + b_j * 3 + 2] -= res_12;
            }
        }

      if (col_l != 0)
        {
          res_00 = 0, res_01 = 0;
          res_10 = 0, res_11 = 0;
          for (int s = 0; s < col_row1; s++)
            {
              res_00 += Block_A[row_k * 3 * col_row + s]
                        * Block_B[s * col_B + col_k * 3];
              if (col_l > 1)
                {
                  res_01 += Block_A[row_k * 3 * col_row + s]
                            * Block_B[s * col_B + col_k * 3 + 1];
                }
              if (row_l > 1)
                {
                  res_10 += Block_A[(row_k * 3 + 1) * col_row + s]
                            * Block_B[s * col_B + col_k * 3];
                }
              if (col_l > 1 && row_l > 1)
                {
                  res_11 += Block_A[(row_k * 3 + 1) * col_row + s]
                            * Block_B[s * col_B + col_k * 3 + 1];
                }
            }

          Result[row_k * 3 * col_B + col_k * 3] -= res_00;

          if (col_l > 1)
            {
              Result[row_k * 3 * col_B + col_k * 3 + 1] -= res_01;
            }

          if (row_l > 1)
            {
              Result[(row_k * 3 + 1) * col_B + col_k * 3] -= res_10;
            }

          if (row_l > 1 && col_l > 1)
            {
              Result[(row_k * 3 + 1) * col_B + col_k * 3 + 1] -= res_11;
            }
        }
    }
}
void
Block_mult_add(double* Result,
                double* Block_A,
                double* Block_B,
                int row_A,
                int col_B,
                int col_row,
                const int m)
{
  if (row_A == 0 || col_B == 0 || col_row == 0)
    return;
  int row_l = row_A % 3;
  int col_l = col_B % 3;
  int row_k = (row_A - row_l) / 3;
  int col_k = (col_B - col_l) / 3;
  double res_00 = 0, res_01 = 0, res_02 = 0;
  double res_10 = 0, res_11 = 0, res_12 = 0;
  double res_20 = 0, res_21 = 0, res_22 = 0;
  int col_row1 = col_row;
  col_B = m;
  row_A = m;
  col_row = m;
  for (int b_i = 0; b_i < row_k; b_i++)
    {
      for (int b_j = 0; b_j < col_k; b_j++)
        {
          res_00 = 0, res_01 = 0, res_02 = 0;
          res_10 = 0, res_11 = 0, res_12 = 0;
          res_20 = 0, res_21 = 0, res_22 = 0;
          for (int s = 0; s < col_row1; s++)
            {
              res_00 += Block_A[b_i * 3 * col_row + s]
                        * Block_B[s * col_B + b_j * 3];
              res_01 += Block_A[b_i * 3 * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 1];
              res_02 += Block_A[b_i * 3 * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 2];

              res_10 += Block_A[(b_i * 3 + 1) * col_row + s]
                        * Block_B[s * col_B + b_j * 3];
              res_11 += Block_A[(b_i * 3 + 1) * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 1];
              res_12 += Block_A[(b_i * 3 + 1) * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 2];

              res_20 += Block_A[(b_i * 3 + 2) * col_row + s]
                        * Block_B[s * col_B + b_j * 3];
              res_21 += Block_A[(b_i * 3 + 2) * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 1];
              res_22 += Block_A[(b_i * 3 + 2) * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 2];
            }

          Result[b_i * 3 * col_B + b_j * 3] += res_00;
          Result[b_i * 3 * col_B + b_j * 3 + 1] += res_01;
          Result[b_i * 3 * col_B + b_j * 3 + 2] += res_02;
          Result[(b_i * 3 + 1) * col_B + b_j * 3] += res_10;
          Result[(b_i * 3 + 1) * col_B + b_j * 3 + 1] += res_11;
          Result[(b_i * 3 + 1) * col_B + b_j * 3 + 2] += res_12;
          Result[(b_i * 3 + 2) * col_B + b_j * 3] += res_20;
          Result[(b_i * 3 + 2) * col_B + b_j * 3 + 1] += res_21;
          Result[(b_i * 3 + 2) * col_B + b_j * 3 + 2] += res_22;
        }

      if (col_l != 0)
        {
          res_00 = 0, res_01 = 0, res_10 = 0;
          res_11 = 0, res_20 = 0, res_21 = 0;

          for (int s = 0; s < col_row1; s++)
            {
              if (col_l > 1)
                {
                  res_01 += Block_A[b_i * 3 * col_row + s]
                            * Block_B[s * col_B + col_k * 3 + 1];
                  res_11 += Block_A[(b_i * 3 + 1) * col_row + s]
                            * Block_B[s * col_B + col_k * 3 + 1];
                  res_21 += Block_A[(b_i * 3 + 2) * col_row + s]
                            * Block_B[s * col_B + col_k * 3 + 1];
                }

              res_00 += Block_A[b_i * 3 * col_row + s]
                        * Block_B[s * col_B + col_k * 3];
              res_10 += Block_A[(b_i * 3 + 1) * col_row + s]
                        * Block_B[s * col_B + col_k * 3];
              res_20 += Block_A[(b_i * 3 + 2) * col_row + s]
                        * Block_B[s * col_B + col_k * 3];
            }

          Result[b_i * 3 * col_B + col_k * 3] += res_00;
          Result[(b_i * 3 + 1) * col_B + col_k * 3] += res_10;
          Result[(b_i * 3 + 2) * col_B + col_k * 3] += res_20;

          if (col_l > 1)
            {
              Result[b_i * 3 * col_B + col_k * 3 + 1] += res_01;
              Result[(b_i * 3 + 1) * col_B + col_k * 3 + 1] += res_11;
              Result[(b_i * 3 + 2) * col_B + col_k * 3 + 1] += res_21;
            }
        }
    }

  if (row_l != 0)
    {
      for (int b_j = 0; b_j < col_k; b_j++)
        {
          res_00 = 0, res_01 = 0, res_02 = 0;
          res_10 = 0, res_11 = 0, res_12 = 0;
          for (int s = 0; s < col_row1; s++)
            {
              res_00 += Block_A[(row_k * 3) * col_row + s]
                        * Block_B[s * col_B + b_j * 3];
              res_01 += Block_A[(row_k * 3) * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 1];
              res_02 += Block_A[(row_k * 3) * col_row + s]
                        * Block_B[s * col_B + b_j * 3 + 2];

              if (row_l > 1)
                {
                  res_10 += Block_A[(row_k * 3 + 1) * col_row + s]
                            * Block_B[s * col_B + b_j * 3];
                  res_11 += Block_A[(row_k * 3 + 1) * col_row + s]
                            * Block_B[s * col_B + b_j * 3 + 1];
                  res_12 += Block_A[(row_k * 3 + 1) * col_row + s]
                            * Block_B[s * col_B + b_j * 3 + 2];
                }
            }

          Result[row_k * 3 * col_B + b_j * 3] += res_00;
          Result[row_k * 3 * col_B + b_j * 3 + 1] += res_01;
          Result[row_k * 3 * col_B + b_j * 3 + 2] += res_02;

          if (row_l > 1)
            {
              Result[(row_k * 3 + 1) * col_B + b_j * 3] += res_10;
              Result[(row_k * 3 + 1) * col_B + b_j * 3 + 1] += res_11;
              Result[(row_k * 3 + 1) * col_B + b_j * 3 + 2] += res_12;
            }
        }

      if (col_l != 0)
        {
          res_00 = 0, res_01 = 0;
          res_10 = 0, res_11 = 0;
          for (int s = 0; s < col_row1; s++)
            {
              res_00 += Block_A[row_k * 3 * col_row + s]
                        * Block_B[s * col_B + col_k * 3];
              if (col_l > 1)
                {
                  res_01 += Block_A[row_k * 3 * col_row + s]
                            * Block_B[s * col_B + col_k * 3 + 1];
                }
              if (row_l > 1)
                {
                  res_10 += Block_A[(row_k * 3 + 1) * col_row + s]
                            * Block_B[s * col_B + col_k * 3];
                }
              if (col_l > 1 && row_l > 1)
                {
                  res_11 += Block_A[(row_k * 3 + 1) * col_row + s]
                            * Block_B[s * col_B + col_k * 3 + 1];
                }
            }

          Result[row_k * 3 * col_B + col_k * 3] += res_00;

          if (col_l > 1)
            {
              Result[row_k * 3 * col_B + col_k * 3 + 1] += res_01;
            }

          if (row_l > 1)
            {
              Result[(row_k * 3 + 1) * col_B + col_k * 3] += res_10;
            }

          if (row_l > 1 && col_l > 1)
            {
              Result[(row_k * 3 + 1) * col_B + col_k * 3 + 1] += res_11;
            }
        }
    }
}