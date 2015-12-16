#include"apspmpi.h"

#define MASTER 0

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	MPI_Apsp(MPI_COMM_WORLD, "100x100vec.bin", "100x100out.bin", 10000);
	MPI_Finalize();
	return 0;
}
