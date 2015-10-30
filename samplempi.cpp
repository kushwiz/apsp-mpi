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
	n = 8;
	int nsquare = pow(n, 2);

	MPI_Status status;

	std::cout<<"My Rank:"<<rank<<std::endl;

	MPI_File thefile;

	MPI_File_open(MPI_COMM_WORLD, "8x1vec.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &thefile);

	MPI_Offset filesize;

	MPI_File_get_size(thefile, &filesize);

	MPI_Datatype MATRIX_INT_VECTOR;

	MPI_Type_vector(n/rootp, n/rootp, n, MPI_INT, &MATRIX_INT_VECTOR);
	MPI_Type_commit(&MATRIX_INT_VECTOR);

	int buffr[(n/rootp)*(n/rootp)];

	int first_offset = (nsquare / rootp) * floor((rank) / rootp);
	int second_offset = (rank%rootp) * n / rootp;

	//if(rank == 0) {
		MPI_File_set_view(thefile, (first_offset + second_offset)*sizeof(int), MPI_INT, MATRIX_INT_VECTOR, "native", MPI_INFO_NULL);
		MPI_File_read_all(thefile, buffr, ((n/rootp)*(n/rootp)), MPI_INT, &status);
		std::cout<<"--"<<rank<<"--";
		for(int i=0;i<((n/rootp) * (n/rootp));i++)
		{
			std::cout<<buffr[i]<<" ";
		}
		std::cout<<std::endl;
	//}

	MPI_File_close(&thefile);
	MPI_Type_free(&MATRIX_INT_VECTOR);
	MPI_Finalize();
	return 0;
}
