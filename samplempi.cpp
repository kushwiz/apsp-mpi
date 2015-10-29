#include<mpi.h>
#include<iostream>
#include<stdlib.h>
#include<math.h>

#define MASTER 0

int main(int argc, char **argv)
{

	int n, rank, size;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int p = size;
	int rootp = sqrt(size);
	n = 800;
	int nsquare = pow(n, 2);

	MPI_Status status;

	std::cout<<"My Rank:"<<rank<<std::endl;

	MPI_File thefile;

	MPI_File_open(MPI_COMM_WORLD, "80x10vec.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &thefile);

	MPI_Offset filesize;

	MPI_File_get_size(thefile, &filesize);

	MPI_Datatype MATRIX_INT_VECTOR;

	MPI_Type_vector(1, n/rootp, n/rootp, MPI_INT, &MATRIX_INT_VECTOR);
	MPI_Type_commit(&MATRIX_INT_VECTOR);

	int buffr[n/rootp];

	int first_offset = (nsquare / rootp) * floor((rank) / rootp);
	int second_offset = (rank%rootp) * n / rootp;

		MPI_File_read_at(thefile, (first_offset+second_offset)*sizeof(int), buffr, 1, MATRIX_INT_VECTOR, &status);
		std::cout<<"--"<<rank<<"--";
		for(int i=0;i<n/rootp;i++)
		{
			std::cout<<buffr[i]<<" ";
		}
		std::cout<<std::endl;
		MPI_File_close(&thefile);

	MPI_Type_free(&MATRIX_INT_VECTOR);
	MPI_Finalize();
	return 0;
}
