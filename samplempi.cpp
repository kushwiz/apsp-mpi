#include<mpi.h>
#include<iostream>
#include<stdlib.h>

#define MASTER 0

int main(int argc, char **argv) 
{
	int n, rank, size;

	MPI::Init(argc, argv);

	size = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();

	MPI::Status status;

	std::cout<<"My Rank:"<<rank<<std::endl;

	MPI::File thefile = MPI::File::Open(MPI::COMM_WORLD, "80x10vec.bin", MPI::MODE_RDONLY, MPI::INFO_NULL);

	MPI::Offset filesize = thefile.Get_size();

	int buffr[400];

	if(rank == MASTER)
	{
		thefile.Read(buffr, 400, MPI::INT, status);
		for(int i=0;i<400;i++)
		{
			std::cout<<buffr[i]<<" ";
		}
		std::cout<<std::endl;
	}
	thefile.Close();

	MPI::Finalize();
	return 0;
}
