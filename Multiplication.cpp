#include "header.h"
void
Block_multiplication(double* Result,
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

          Result[b_i * 3 * col_B + b_j * 3] = res_00;
          Result[b_i * 3 * col_B + b_j * 3 + 1] = res_01;
          Result[b_i * 3 * col_B + b_j * 3 + 2] = res_02;
          Result[(b_i * 3 + 1) * col_B + b_j * 3] = res_10;
          Result[(b_i * 3 + 1) * col_B + b_j * 3 + 1] = res_11;
          Result[(b_i * 3 + 1) * col_B + b_j * 3 + 2] = res_12;
          Result[(b_i * 3 + 2) * col_B + b_j * 3] = res_20;
          Result[(b_i * 3 + 2) * col_B + b_j * 3 + 1] = res_21;
          Result[(b_i * 3 + 2) * col_B + b_j * 3 + 2] = res_22;
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
          Result[b_i * 3 * col_B + col_k * 3] = res_00;
          Result[(b_i * 3 + 1) * col_B + col_k * 3] = res_10;
          Result[(b_i * 3 + 2) * col_B + col_k * 3] = res_20;

          if (col_l > 1)
            {
              Result[b_i * 3 * col_B + col_k * 3 + 1] = res_01;
              Result[(b_i * 3 + 1) * col_B + col_k * 3 + 1] = res_11;
              Result[(b_i * 3 + 2) * col_B + col_k * 3 + 1] = res_21;
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

          Result[row_k * 3 * col_B + b_j * 3] = res_00;
          Result[row_k * 3 * col_B + b_j * 3 + 1] = res_01;
          Result[row_k * 3 * col_B + b_j * 3 + 2] = res_02;

          if (row_l > 1)
            {
              Result[(row_k * 3 + 1) * col_B + b_j * 3] = res_10;
              Result[(row_k * 3 + 1) * col_B + b_j * 3 + 1] = res_11;
              Result[(row_k * 3 + 1) * col_B + b_j * 3 + 2] = res_12;
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

          Result[row_k * 3 * col_B + col_k * 3] = res_00;

          if (col_l > 1)
            {
              Result[row_k * 3 * col_B + col_k * 3 + 1] = res_01;
            }

          if (row_l > 1)
            {
              Result[(row_k * 3 + 1) * col_B + col_k * 3] = res_10;
            }

          if (row_l > 1 && col_l > 1)
            {
              Result[(row_k * 3 + 1) * col_B + col_k * 3 + 1] = res_11;
            }
        }
    }
}
