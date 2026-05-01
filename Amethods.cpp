#include "header.h"
void
  Ainit(double* a, int n, int m, int p, int k, int s)
{
  int i_loc, j_loc, i_glob, j_glob;
  int rows = get_rows(n, m, p, k);
  for (i_loc = 0; i_loc < rows; ++i_loc)
    {
      i_glob = l_to_g(n, m, p, k, i_loc);
      for (j_loc = 0; j_loc < n; ++j_loc)
        {
          j_glob = j_loc;
          a[i_loc * n + j_loc] = f(i_glob, j_glob, n, s);
        }
    }
}
double
  f(int i, int j, int n, int s)
{
  switch (s)
    {
    case 1:
      return (n - std::max(i, j));
    case 2:
      return std::max(i, j) + 1;
    case 3:
      return fabs(i - j);
    case 4:
      return (1.0 / (i + j + 1.0));
    default:
      return 0.0;
    }
}

void
  print_matrix(double* A, int n, int m, int p, int k, double* buf, int max_print, MPI_Comm com)
{
  int main_k = 0;
  int b, b_max = (n + m - 1) / m, printed_rows = 0;
  for (b = 0; b < b_max; b++)
    {
      int owner = b % p, b_loc = b / p;
      int rows = std::min(m, n - b * m);
      if (k == main_k)
        {
          if (owner == main_k)
            {
              printed_rows += print_array(
                  A + b_loc * n * m, n, rows, printed_rows, max_print);
            }
          else
            {
              MPI_Status st;
              MPI_Recv(buf, n * rows, MPI_DOUBLE, owner, 0, com, &st);
              printed_rows
                  += print_array(buf, n, rows, printed_rows, max_print);
            }
        }
      else if (owner == k)
        {
          MPI_Send(A + b_loc * n * m, n * rows, MPI_DOUBLE, main_k, 0, com);
        }
    }
}
int
  print_array(double* A, int n, int m, int printed_rows, int max_print)
{
  if (printed_rows >= max_print)
    {
      return 0;
    }
  int p_n = (n > max_print ? max_print : n);
  int p_m = (printed_rows + m <= max_print ? m : max_print - printed_rows);
  for (int i = 0; i < p_m; ++i)
    {
      for (int j = 0; j < p_n; ++j)
        {
          printf(" %10.3e", A[i * n + j]);
        }
      printf("\n");
    }
  return p_m;
}

int
  read_matrix(double* A, int n, int m, int p, int k, const std::string& name, double* buf, MPI_Comm com)
{
  int main_k = 0, err = 0;
  std::ifstream in;

  if (k == main_k)
    {
      in.open (name);
      if (!in.is_open ())
        {
          err = 1;
        }
    }
  MPI_Bcast(&err, 1, MPI_INT, main_k, com);
  if (err)
    return err;

  int b, max_b = (n + m - 1) / m;
  memset(buf, 0, n * m * sizeof(double));
  for (b = 0; b < max_b; b++)
    {
      int owner = b % p;
      int rows = (b * m + m < n ? m : n - b * m); // b*m +m??
      int b_loc = b / p;

      if (k == main_k)
        {
          err += read_array (in, buf, n * rows);
          if (owner == main_k)
            {
              memcpy (A + b_loc * n * m, buf, n * rows * sizeof(double));
            }
          else
            {
              MPI_Send (buf, n * rows, MPI_DOUBLE, owner, 0, com);
            }
        }
      else if (owner == k)
        {
          MPI_Status st;
          MPI_Recv (A + b_loc * n * m, n * rows, MPI_DOUBLE, main_k, 0, com, &st);
        }
    }
  if (k == 0)
    in.close();

  MPI_Bcast(&err, 1, MPI_INT, main_k, com);

  if (err)
    return err;

  return 0;
}

int
  read_array (std::ifstream& in, double* buf, int n)
{
  double x;
  for (int i = 0; i < n; i++)
    {
      if (!(in >> x))
        return 10;
      buf[i] = x;
    }
  return 0;
}
