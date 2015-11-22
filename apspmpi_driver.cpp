#include"apspmpi.h"

#define MASTER 0

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
	MPI_Apsp(MPI_COMM_WORLD, 8, "8x1vec.bin", "8x1out.bin");
	return 0;
}
