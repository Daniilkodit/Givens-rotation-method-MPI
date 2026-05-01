#include "header.h"

int
  l_to_g(int /* n */, int m, int p, int k, int i_loc)
{
  int i_loc_m = i_loc / m;
  int i_glob_m = i_loc_m * p + k;
  return i_glob_m * m + (i_loc % m);
}

int
  l_to_g_block(int /* n */, int /* m */, int p, int k, int i_loc_m)
{
  return i_loc_m * p + k;
}

int
  g_to_l(int /* n */, int m, int p, int /* k */, int i_glob)
{
  int i_glob_m = i_glob / m;
  int i_loc_m = i_glob_m / p;
  return i_loc_m * m + (i_glob % m);
}

int
  g_to_l_block(int /* n */, int /* m */, int p, int /* k */, int i_glob_m)
{
  return i_glob_m / p;
}

int
get_max_block_rows(int n, int m, int p)
{
  // count block rows
  int b = (m + n - 1) / m;
  // ( b % p == 0 ? b / p : b / p + 1 )
  return (b + p - 1) / p;
}

int
  get_k(int /*n*/, int m, int p, int i_glob)
{
  int i_glob_m = i_glob / m;
  return i_glob_m % p;
}

int
  get_k_block(int /*n*/, int /*m*/, int p, int i_glob_m)
{
  return i_glob_m % p;
}

int
  get_block_rows(int n, int m, int p, int k)
{
  int b = (n + m - 1) / m;
  return ((b % p) <= k ? b / p : b / p + 1);
}

int
  get_rows(int n, int m, int p, int k)
{
  int b = (n + m - 1) / m;
  int b_last = (b - 1) % p;
  int b_loc = ((b % p) <= k ? b / p : b / p + 1);

  if (b_last != k)
    {
      return b_loc * m;
    }

  int l = n % m;
  if (l == 0)
    {
      return b_loc * m;
    }

  if (b_loc == 0)
    {
      return 0;
    }

  return (b_loc - 1) * m + l;
}
