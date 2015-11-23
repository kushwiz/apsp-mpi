#include"apspmpi.h"

#define MASTER 0

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	MPI_Apsp(MPI_COMM_WORLD, "80x10vec.bin", "80x10out.bin");
	MPI_Finalize();
	return 0;
}
