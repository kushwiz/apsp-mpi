#include"apspmpi.h"

#define MASTER 0

void MPI_Apsp(MPI_Comm COMM, std::string infile_path, std::string outfile_path)
{
	int rank, size, row_rank, col_rank, root_pe_count, n_by_rootp;
	long nodes, nsquare, nsquare_by_rootp;

	std::ios_base::sync_with_stdio(false);
	std::cin.tie(0);

	MPI_Comm_rank(COMM, &rank);
	MPI_Comm_size(COMM, &size);

	MPI_Comm ROW_COMM, COL_COMM;
	MPI_Status status;
	MPI_File infile, outfile;
	MPI_Offset filesize;
	MPI_Datatype MATRIX_INT_VECTOR;

	MPI_File_open(COMM, infile_path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &infile);
	MPI_File_get_size(infile, &filesize);
	nodes = sqrt(filesize/sizeof(int));
	if(rank==MASTER) {
		std::cout<<"No. of nodes:"<<nodes<<std::endl;
	}

	root_pe_count = sqrt(size);
	nsquare = nodes * nodes;
	nsquare_by_rootp = nsquare / root_pe_count;
	n_by_rootp = nodes / root_pe_count;

	MPI_Type_vector(n_by_rootp, n_by_rootp, nodes, MPI_INT, &MATRIX_INT_VECTOR);
	MPI_Type_commit(&MATRIX_INT_VECTOR);

	std::vector<int> matrixBuffer(n_by_rootp * n_by_rootp);

	int pe_layout_row_offset = floor(rank/root_pe_count);
	int pe_layout_col_offset = rank%root_pe_count;
	long row_offset = (nsquare / root_pe_count) * floor(rank / root_pe_count);
	long col_offset = (rank%root_pe_count) * n_by_rootp;

	std::cout<<"row offset:"<<row_offset<<" coloffset:"<<col_offset<<std::endl;

	MPI_Comm_split(COMM, pe_layout_row_offset, rank, &ROW_COMM);
	MPI_Comm_split(COMM, pe_layout_col_offset, rank, &COL_COMM);	

	MPI_Comm_rank(ROW_COMM, &row_rank);
	MPI_Comm_rank(COL_COMM, &col_rank);

	double t0 = MPI_Wtime();

	MPI_File_set_view(infile, (row_offset + col_offset)*sizeof(int), MPI_INT, MATRIX_INT_VECTOR, "native", MPI_INFO_NULL);
	MPI_File_read(infile, matrixBuffer.data(), matrixBuffer.size(),  MPI_INT, &status);

	double t1 = MPI_Wtime();

	if(rank==MASTER) {
		std::cout<<"File reading time:"<<(t1-t0)<<std::endl;
	}

	std::vector<int> rowDataBuffer(n_by_rootp);
	std::vector<int> colDataBuffer(n_by_rootp);
	std::vector<int> myRowData(n_by_rootp);
	std::vector<int> myColData(n_by_rootp);

	for(long k=0;k<nodes;k++)
	{
		if(pe_layout_row_offset == floor(k/n_by_rootp) && (pe_layout_col_offset == floor(k/n_by_rootp))) {
			// make copy of row to send
			memcpy(&rowDataBuffer[0], &matrixBuffer[(k%n_by_rootp)*n_by_rootp], n_by_rootp*sizeof(int));
			// Broadcast ROW data to ROW_COMM
			MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, pe_layout_row_offset, ROW_COMM);
			// make copy of column to send
			for(long i=0, j=(k%n_by_rootp); i<n_by_rootp; i++, j+=n_by_rootp) {
				colDataBuffer[i] = matrixBuffer[j];
			}
			// Broadcast COL data to COL_COMM
			MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, pe_layout_col_offset, COL_COMM);
		} else if(pe_layout_row_offset == floor(k/n_by_rootp)) {
			// receive ROW from  row_offset(th) RANK.
			MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, pe_layout_row_offset, ROW_COMM);
			// Store my column in ColDataBuffer
			memcpy(&myRowData[0], &matrixBuffer[(k%n_by_rootp)*n_by_rootp], n_by_rootp*sizeof(int));
			// Broadcast to my COL COMM.
			MPI_Bcast(myRowData.data(), myRowData.size(), MPI_INT, col_rank, COL_COMM);
		} else if(pe_layout_col_offset == floor(k/n_by_rootp)) {
			// receive COL from col_offset(th) RANK.
			MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, pe_layout_col_offset, COL_COMM);

			// Store my Row in RowDataBuffer	
			memcpy(&rowDataBuffer[0], &matrixBuffer[(k%n_by_rootp)*n_by_rootp], n_by_rootp*sizeof(int));

			// Broadcast my Col data to my ROW COMM.
			for(long i=0, j=(k%n_by_rootp); i<n_by_rootp; i++, j+=n_by_rootp) {
				myColData[i] = matrixBuffer[j];
			}
			MPI_Bcast(myColData.data(), myColData.size(), MPI_INT, row_rank, ROW_COMM);
		} else {
			// Receive COL from ROW COMM.
			MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, floor(k/n_by_rootp), ROW_COMM);

			// Receive ROW from COL COMM.
			MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, floor(k/n_by_rootp), COL_COMM);
		}			

		MPI_Barrier(MPI_COMM_WORLD);

		int ab, axb, min;
		for(long i=0;i<n_by_rootp;i++) 
		{
			for(long j=0;j<n_by_rootp;j++) 
			{
				if(rowDataBuffer[i]==0 || colDataBuffer[j]==0) {
					axb = std::numeric_limits<int>::max();
				} else {
					axb = rowDataBuffer[i] + colDataBuffer[j];
				}
				if(matrixBuffer[i*n_by_rootp+j]==0)
				{
					ab = std::numeric_limits<int>::max();
				} else {
					ab = matrixBuffer[i*n_by_rootp+j];
				}
				min = std::min(ab, axb);
				matrixBuffer[i*n_by_rootp+j] = (min==std::numeric_limits<int>::max())?0:min;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	} //end of k for loop

	double t2 = MPI_Wtime();
	if(rank==MASTER) {
		std::cout<<"Running time:" <<(t2-t0)<<std::endl;
	}

	MPI_File_open(COMM, outfile_path.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &outfile);
	MPI_File_set_view(outfile, (row_offset + col_offset)*sizeof(int), MPI_INT, MATRIX_INT_VECTOR, "native", MPI_INFO_NULL);
	if(rank==0) {
		std::cout<<"matrix size:"<<matrixBuffer.size()<<std::endl;
	}
	MPI_File_write(outfile, matrixBuffer.data(), matrixBuffer.size(),  MPI_INT, &status);

	double t3 = MPI_Wtime();
	if(rank==MASTER) {
		std::cout<<"File Writing time:"<<t3-t2<<std::endl;
	}

	MPI_File_close(&infile);
	MPI_File_close(&outfile);
	MPI_Comm_free(&ROW_COMM);
	MPI_Comm_free(&COL_COMM);
	MPI_Type_free(&MATRIX_INT_VECTOR);
}
