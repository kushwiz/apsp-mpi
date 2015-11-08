#include<mpi.h>
#include<iostream>
#include<vector>
#include<stdlib.h>
#include<math.h>
#include<limits>

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

	MPI_File thefile;

	MPI_File_open(MPI_COMM_WORLD, "8x1vec.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &thefile);

	MPI_Offset filesize;

	MPI_File_get_size(thefile, &filesize);

	MPI_Datatype MATRIX_INT_VECTOR;

	MPI_Type_vector(n/rootp, n/rootp, n, MPI_INT, &MATRIX_INT_VECTOR);
	MPI_Type_commit(&MATRIX_INT_VECTOR);

	std::vector<int> buffr((n/rootp)*(n/rootp));

	int row_offset = (nsquare / rootp) * floor((rank) / rootp);
	int col_offset = (rank%rootp) * n / rootp;
        
	MPI_Comm_split(MPI_COMM_WORLD, row_offset, 0, &ROW_COMM);
 	MPI_Comm_split(MPI_COMM_WORLD, col_offset, 0, &COL_COMM);	
        double t0 = MPI_Wtime();
	MPI_Comm_rank(ROW_COMM, &row_rank);
	MPI_Comm_rank(COL_COMM, &col_rank);

		MPI_File_set_view(thefile, (row_offset + col_offset)*sizeof(int), MPI_INT, MATRIX_INT_VECTOR, "native", MPI_INFO_NULL);
		MPI_File_read_all(thefile, buffr.data(), buffr.size(),  MPI_INT, &status);
		//std::cout<<"--"<<rank<<"--"<<row_offset<<"--"<<col_offset<<"--";
		/*for(int i=0;i<((n/rootp) * (n/rootp));i++)
		{
			std::cout<<buffr[i]<<" ";
		}*/
		//std::cout<<std::endl;

		int pe_layout_row_offset = floor(rank/rootp);
		int pe_layout_col_offset = rank%rootp;

                
		//int k = 0;
		int nbyrootp = n/rootp;

		std::vector<int> rowDataBuffer(nbyrootp);
		std::vector<int> colDataBuffer(nbyrootp);
		std::vector<int> myRowData(nbyrootp);
		std::vector<int> myColData(nbyrootp);
                for( int k=0;k<n;k++)
                {
		if(pe_layout_row_offset == floor(k/nbyrootp) && (pe_layout_col_offset == floor(k/nbyrootp))) {
			// make copy of row to send
			memcpy(&rowDataBuffer[0], &buffr[(k%nbyrootp)*nbyrootp], nbyrootp*sizeof(int));
			// Broadcast ROW data to ROW_COMM
			MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, pe_layout_row_offset, ROW_COMM);
			// make copy of column to send
			
			for(int i=0, j=(k%nbyrootp); i<n/rootp; i++, j+=nbyrootp) {
				colDataBuffer[i] = buffr[j];
			}

			// Broadcast COL data to COL_COMM
			MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, pe_layout_col_offset, COL_COMM);

		} else if(pe_layout_row_offset == floor(k/nbyrootp)) {
			// receive ROW from  row_offset(th) RANK.
			MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, pe_layout_row_offset, ROW_COMM);
			// Store my column in ColDataBuffer
                        /*
		       	for(int i=0,j=(k%nbyrootp);i<nbyrootp;i++,j+=nbyrootp) {
				colDataBuffer[i] = buffr[j];
			}
                        */
			// Broadcast to my COL COMM.
			memcpy(&myRowData[0], &buffr[(k%nbyrootp)*nbyrootp], nbyrootp*sizeof(int));
			MPI_Bcast(myRowData.data(), myRowData.size(), MPI_INT, col_rank, COL_COMM);
		} else if(pe_layout_col_offset == floor(k/nbyrootp)) {
			// receive COL from col_offset(th) RANK.
			MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, pe_layout_col_offset, COL_COMM);

			// Store my Row in RowDataBuffer	
			memcpy(&rowDataBuffer[0], &buffr[(k%nbyrootp)*nbyrootp], nbyrootp*sizeof(int));
			
			// Broadcast my Col data to my ROW COMM.
			for(int i=0, j=(k%nbyrootp); i<nbyrootp; i++, j+=nbyrootp) {
				myColData[i] = buffr[j];
			}


			MPI_Bcast(myColData.data(), myColData.size(), MPI_INT, row_rank, ROW_COMM);
		} else {
			// Receive COL from ROW COMM.
			MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, floor(k/nbyrootp), ROW_COMM);
			/*std::cout<<"new recv:"<<rank<<"---";
	  		for(int i=0;i<colDataBuffer.size();i++) {
	 			std::cout<<colDataBuffer[i]<<" ";
			}*/

			// Receive ROW from COL COMM.
			MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, floor(k/nbyrootp), COL_COMM);
		/*	std::cout<<"\nrow new recv:"<<rank<<"---";
			for(int i=0;i<rowDataBuffer.size();i++) {
				std::cout<<rowDataBuffer[i]<<" ";
			}*/
		}		
		
		MPI_Barrier(MPI_COMM_WORLD);

			//std::cout<<"Bellman Ford for rank:"<<rank;
			int temp, temp1;
			for(int i=0;i<nbyrootp;i++) 
			{
				for(int j=0;j<nbyrootp;j++) 
				{
					//std::cout<<"--min("<<buffr[i*nbyrootp+j]<<","<<colDataBuffer[i]<<"+"<<rowDataBuffer[j]<<") ";
					if(rowDataBuffer[j]==0 || colDataBuffer[i]==0) {
						temp = std::numeric_limits<int>::max();
					} else {
						temp = rowDataBuffer[j] + colDataBuffer[i];
					}
					if(buffr[i*nbyrootp+j]==0)
					{
						temp1 = std::numeric_limits<int>::max();
					} else {
						temp1 = buffr[i*nbyrootp+j];
					}
					//std::cout<<((std::min(temp1, temp)==std::numeric_limits<int>::max())?0:std::min(temp1, temp))<<">>";
					buffr[i*nbyrootp+j] = ((std::min(temp1, temp)==std::numeric_limits<int>::max())?0:std::min(temp1, temp));
				}
			}
			//std::cout<<std::endl;
			MPI_Barrier(MPI_COMM_WORLD);
		} //end of k for loop
        MPI_Barrier(MPI_COMM_WORLD);
        double t1 = MPI_Wtime();
        if(rank==0)
        {
            std::cout<<"Running time:" <<(t1-t0)<<std::endl;
        }
	//}

	MPI_File_close(&thefile);
	MPI_Type_free(&MATRIX_INT_VECTOR);
	MPI_Finalize();
	return 0;
}
