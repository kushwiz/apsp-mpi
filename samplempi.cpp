#include<mpi.h>
#include<iostream>
#include<vector>
#include<stdlib.h>
#include<math.h>

#define MASTER 0

int main(int argc, char **argv)
{

	int n, rank, size, row_rank, col_rank;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_Comm ROW_COMM, COL_COMM;


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

	//int buffr[(n/rootp)*(n/rootp)];
	std::vector<int> buffr((n/rootp)*(n/rootp));

	int row_offset = (nsquare / rootp) * floor((rank) / rootp);
	int col_offset = (rank%rootp) * n / rootp;

	MPI_Comm_split(MPI_COMM_WORLD, row_offset, 0, &ROW_COMM);
 	MPI_Comm_split(MPI_COMM_WORLD, col_offset, 0, &COL_COMM);	

	MPI_Comm_rank(ROW_COMM, &row_rank);
	MPI_Comm_rank(COL_COMM, &col_rank);

	//if(rank == 0) {
		MPI_File_set_view(thefile, (row_offset + col_offset)*sizeof(int), MPI_INT, MATRIX_INT_VECTOR, "native", MPI_INFO_NULL);
		MPI_File_read_all(thefile, buffr.data(), buffr.size(),  MPI_INT, &status);
		std::cout<<"--"<<rank<<"--"<<row_offset<<"--"<<col_offset<<"--";
		for(int i=0;i<((n/rootp) * (n/rootp));i++)
		{
			std::cout<<buffr[i]<<" ";
		}
		std::cout<<std::endl;

		int pe_layout_row_offset = floor(rank/rootp);
		int pe_layout_col_offset = rank%rootp;


		int k = 0;
		int nbyrootp = n/rootp;

		std::vector<int> rowDataBuffer(n/rootp);
		std::vector<int> colDataBuffer(n/rootp);
		std::vector<int> myRowData(n/rootp);
		std::vector<int> myColData(nbyrootp);

		if(pe_layout_row_offset == floor(k/nbyrootp) && (pe_layout_col_offset == floor(k/nbyrootp))) {
			std::cout<<"world rank:"<<rank<<"---";
			// make copy of row to send
			memcpy(&rowDataBuffer[0], &buffr[(k%4)*nbyrootp], nbyrootp*sizeof(int));
//			for(int i=0;i<rowDataBuffer.size();i++) {
	//			std::cout<<rowDataBuffer[i]<<" ";
	//		}
			std::cout<<std::endl;
			// Broadcast ROW data to ROW_COMM
			MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, pe_layout_row_offset, ROW_COMM);
			// make copy of column to send
			
			for(int i=0, j=(k%nbyrootp); i<n/rootp; i++, j+=n/rootp) {
				colDataBuffer[i] = buffr[j];
			}

			//std::cout<<"Col data:";
			//for(int i=0; i<colDataBuffer.size(); i++) {
			//	std::cout<<colDataBuffer[i]<<" ";
			//}

			// Broadcast COL data to COL_COMM
			MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, pe_layout_col_offset, COL_COMM);

		} else if(pe_layout_row_offset == floor(k/nbyrootp)) {
			// receive ROW from  row_offset(th) RANK.
			MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, pe_layout_row_offset, ROW_COMM);
			//std::cout<<"row Recd. data"<<rank<<"---";
			//for(int i=0;i<rowDataBuffer.size();i++) {
			//	std::cout<<rowDataBuffer[i]<<" ";
			//}
			// keep a copy for myself
			// Broadcast to my COL COMM.
			memcpy(&myRowData[0], &buffr[(k%4)*nbyrootp], nbyrootp*sizeof(int));
			MPI_Bcast(myRowData.data(), myRowData.size(), MPI_INT, col_rank, COL_COMM);
		} else if(pe_layout_col_offset == floor(k/nbyrootp)) {
			// receive COL from col_offset(th) RANK.
			MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, pe_layout_col_offset, COL_COMM);
//			std::cout<<"world rank:"<<rank<<"recv Col data:";
//			for(int i=0; i<colDataBuffer.size(); i++) {
//				std::cout<<colDataBuffer[i]<<" ";
//			}

			// keep a copy for myself
			// Broadcast my Col data to my ROW COMM.
			for(int i=0, j=(k%nbyrootp); i<nbyrootp; i++, j+=nbyrootp) {
				myColData[i] = buffr[j];
			}


			MPI_Bcast(myColData.data(), myColData.size(), MPI_INT, row_rank, ROW_COMM);
		} else {
			// Receive COL from ROW COMM.
			MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, floor(k/nbyrootp), ROW_COMM);
			std::cout<<"new recv:"<<rank<<"---";
			for(int i=0;i<colDataBuffer.size();i++) {
				std::cout<<colDataBuffer[i]<<" ";
			}

			MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, floor(k/nbyrootp), COL_COMM);
			std::cout<<"\nrow new recv:"<<rank<<"---";
			for(int i=0;i<rowDataBuffer.size();i++) {
				std::cout<<rowDataBuffer[i]<<" ";
			}


			// Receive ROW from COL COMM.
		}		
		
		//MPI_Barrier(MPI_COMM_WORLD);

	//}

	MPI_File_close(&thefile);
	MPI_Type_free(&MATRIX_INT_VECTOR);
	MPI_Finalize();
	return 0;
}
