#include "header.h"
int main (int argc, char* argv[])
{

  MPI_Init (&argc, &argv);
  my_main (argc, argv);
  MPI_Finalize ();

  return 0;
}
