#include "header.h"
int
  my_main (int argc, char* argv[])
{
  int n, m, p, r, s, task = 22, res = 0, id_process = 0;
  std::string name;
  double eps = 1.e-15;
  double* A = nullptr;
  double* X = nullptr;
  double* buf1 = nullptr;
  double* buf2 = nullptr;
  double t1 = 0, t2 = 0, r1 = -1, r2 = -1, norm_mat = 1;

  MPI_Comm com = MPI_COMM_WORLD;
  MPI_Comm_rank (com, &id_process);
  MPI_Comm_size (com, &p);

  feenableexcept (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

  if (1
        && !((argc == 5 || argc == 6)
        && sscanf(argv[1], "%d", &n) == 1
        && sscanf(argv[2], "%d", &m) == 1
        && sscanf(argv[3], "%d", &r) == 1
        && sscanf(argv[4], "%d", &s) == 1))
    {
      if(id_process == 0)
        {
          printf ("Usage: %s n m r s\n", argv[0]);
        }

      return 0;
    }

  if (n <= 0 || m <= 0 || p <= 0 || r <= 0 || m > n)
    {
      if(id_process == 0)
        {
          printf ("Usage: n>0 m>0 p>0 r>0 s and m<=n\n");
        }

      return 0;
    }

  if (s == 0)
    {
      if (argc != 6)
        {
          if(id_process == 0)
            {
              printf ("Usage: %s n m r 0 namefile.txt\n", argv[0]);
            }

          return 0;
        }
      else
        name = argv[5];
    }

  int rows = get_max_block_rows (n, m, p) * m;
  int rows_id_proc = get_rows (n, m, p, id_process);
  A = new double[rows * n]{};

  if (!A)
    {
      printf ("failed to allocate memory A\n");
      return -1;
    }
  X = new double[rows * n]{};

  if (!X)
    {
      printf ("failed to allocate memory X\n");
      delete[] A;
      return -1;
    }

  buf1 = new double[n * m]; // one block row
  buf2 = new double[n * m];

  if (s == 4)
    eps = 1.e-20;

  if (s == 0)
    {
      int err;
      err = read_matrix (A, n, m, p, id_process, name, buf1, com);

      if (err == 1)
        {
          if (id_process == 0)
            {
              printf ("Not found file %s\n", argv[5]);
            }

          delete[] A;
          delete[] X;
          delete[] buf1;
          delete[] buf2;

          return 0;
        }
      else if (err > 1)
        {
          if (id_process == 0)
            {
              printf ("File %s has trash\n", argv[5]);
            }

          delete[] A;
          delete[] X;
          delete[] buf1;
          delete[] buf2;

          return 0;
        }
    }
  else
    {
      Ainit (A, n, m, p, id_process, s);
    }

  print_matrix (A, n, m, p, id_process,  buf1, r, com);

  double* block = nullptr;
  double* rcos = nullptr;
  double* rsin = nullptr;
  double* blockE = nullptr;
  double* blockC = nullptr;
  double* blockCE = nullptr;
  double* rmax = nullptr;
  double* rmax_dop = nullptr;

  for (int i_loc = 0; i_loc < rows_id_proc; i_loc++)
    {
      int j_glob = l_to_g (n, m, p, id_process, i_loc);
      X[i_loc * n + j_glob] = 1;
    }

  block = new double[m * m];
  blockE = new double[m * m];
  blockCE = new double[m * m];
  blockC = new double[m * m];
  rcos = new double[m * m];
  rsin = new double[m * m];
  rmax = new double[n]{};
  rmax_dop = new double[n]{};
  norm_mat = Norm_matrix (A, n, m, p, id_process, rmax, com, rmax_dop);

  t1 = get_full_time ();

  res = Rotation(A, X, n, m, id_process, p, eps, block, blockE, blockC, blockCE, rcos,
                 rsin, norm_mat, com, buf1, buf2);

  t1 = get_full_time () - t1;


  if (res > 0)
    {
      r1 = -1;
      r2 = -1;
    }
  else
    {
      if (s == 0)
        {
          read_matrix (A, n, m, p, id_process, name, buf1, com);
        }
      else
        {
          Ainit (A, n, m, p, id_process, s);
        }

      if (id_process == 0)
        printf ("\n");

      print_matrix (X, n, m, p, id_process, buf1, r, com);

      t2 = get_full_time ();
      r1 = Discrepancy (A, X, n, m, id_process, p, blockCE, blockE, block, rmax, com, rmax_dop);
      t2 = get_full_time ()- t2;

      r2 = Discrepancy (X, A, n, m, id_process, p, blockCE, blockE, block, rmax, com, rmax_dop);
    }



  if (id_process == 0)
    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
            argv[0], task, r1, r2, t1, t2, s, n, m, p);

  delete[] rcos;
  delete[] rsin;
  delete[] blockCE;
  delete[] blockE;
  delete[] block;
  delete[] blockC;
  delete[] A;
  delete[] X;
  delete[] rmax;
  delete[] rmax_dop;
  delete[] buf1;
  delete[] buf2;

  return 0;
}
