#include<mpi.h>
#include<iostream>
#include<stdlib.h>
#include<math.h>

#define MASTER 0

int main(int argc, char **argv)
{

	int n, rank, size;

	MPI::Init(argc, argv);

	size = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();

	int p = size;
	int rootp = sqrt(size);
	n = 800;
	int nsquare = pow(n, 2);

	MPI::Status status;

	std::cout<<"My Rank:"<<rank<<std::endl;

	MPI::File thefile = MPI::File::Open(MPI::COMM_WORLD, "80x10vec.bin", MPI::MODE_RDONLY, MPI::INFO_NULL);

	MPI::Offset filesize = thefile.Get_size();

	MPI::Datatype matrixIntVector = MPI::INT.Create_vector(1, n/rootp, n/rootp);
	matrixIntVector.Commit();

	int buffr[n/rootp];

	int first_offset = (nsquare / rootp) * floor((rank) / rootp);
	int second_offset = (rank%rootp) * n / rootp;

		thefile.Read_at((first_offset + second_offset)*sizeof(int), buffr, 1, matrixIntVector, status);
		for(int i=0;i<n/rootp;i++)
		{
			std::cout<<buffr[i]<<" ";
		}
		std::cout<<std::endl;
	thefile.Close();

	MPI::Finalize();
	return 0;
}
