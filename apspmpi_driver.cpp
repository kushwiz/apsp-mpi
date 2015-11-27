#include"apspmpi.h"

#define MASTER 0

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	MPI_Apsp(MPI_COMM_WORLD, "9x1vec.bin", "9x1out.bin");
	MPI_Finalize();
	return 0;
}
