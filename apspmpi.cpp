#include"apspmpi.h"

#define MASTER 0

enum Color {INNER, RIGHT, LOWER, CORNER};

void MPI_Apsp(MPI_Comm COMM, std::string infile_path, std::string outfile_path, long nodes)
{
	int rank, size, row_rank, col_rank, root_pe_count, n_by_rootp, n_mod_rootp, buffer_row_size, buffer_col_size;
	long nsquare, nsquare_by_rootp;

	std::ios_base::sync_with_stdio(false);
	std::cin.tie(0);

	MPI_Comm_rank(COMM, &rank);
	MPI_Comm_size(COMM, &size);

	MPI_Comm ROW_COMM, COL_COMM;
	MPI_Status status;
	MPI_File infile, outfile;
	MPI_Offset filesize;
	MPI_Datatype MATRIX_INT_INNER_VECTOR, MATRIX_INT_RIGHT_EDGE_VECTOR, MATRIX_INT_LOWER_RIGHT_CORNER_VECTOR, MATRIX_INT_LOWER_EDGE_VECTOR;
	Color pe_position;

	MPI_File_open(COMM, infile_path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &infile);
	if(rank==MASTER) {
		std::cout<<"No. of nodes:"<<nodes<<std::endl;
	}
	root_pe_count = sqrt(size);
	nsquare = nodes * nodes;
	nsquare_by_rootp = nsquare / root_pe_count;
	n_by_rootp = nodes / root_pe_count;
	n_mod_rootp = nodes % root_pe_count;

	int pe_layout_row_offset = floor(rank/root_pe_count);
	int pe_layout_col_offset = rank%root_pe_count;
	long row_offset = (nodes / root_pe_count) * nodes * floor(rank / root_pe_count);
	long col_offset = (rank%root_pe_count) * n_by_rootp;

	//std::cout<<"rank: "<<rank<<" | row offset:"<<row_offset<<" coloffset:"<<col_offset<<std::endl;

	if(pe_layout_row_offset == root_pe_count-1 && pe_layout_col_offset == root_pe_count-1) {
		// Lower right corner PE that takes extra data in both col and row.
		buffer_row_size = n_by_rootp + n_mod_rootp;
		buffer_col_size = n_by_rootp + n_mod_rootp;

		pe_position = CORNER;


	} else if(pe_layout_row_offset == root_pe_count-1) {
		// The last row of PEs that take extra data in col.
		buffer_col_size = n_by_rootp + n_mod_rootp;
		buffer_row_size = n_by_rootp;

		pe_position = LOWER;


	} else if(pe_layout_col_offset == root_pe_count-1) {
		// the right most PEs that take the extra data
		buffer_col_size = n_by_rootp;
		buffer_row_size = n_by_rootp + n_mod_rootp;

		pe_position = RIGHT;

	} else {
		// Inner PEs, which have equally divided rows and columns of data
		buffer_row_size = n_by_rootp;
		buffer_col_size = n_by_rootp;

		pe_position = INNER;

	}

	MPI_Barrier(MPI_COMM_WORLD);

	std::cout<<"rank:"<<rank<<" col_size:"<<buffer_col_size<<" row_size:"<<buffer_row_size<<std::endl;
	std::vector<int> matrixBuffer(buffer_row_size * buffer_col_size);

	MPI_Comm_split(COMM, pe_layout_row_offset, rank, &ROW_COMM);
	MPI_Comm_split(COMM, pe_layout_col_offset, rank, &COL_COMM);	

	MPI_Comm_rank(ROW_COMM, &row_rank);
	MPI_Comm_rank(COL_COMM, &col_rank);

	double t0 = MPI_Wtime();

	switch(pe_position) {
		case INNER:
		MPI_Type_vector(buffer_col_size, buffer_row_size, nodes, MPI_INT, &MATRIX_INT_INNER_VECTOR);
		MPI_Type_commit(&MATRIX_INT_INNER_VECTOR);
		MPI_File_set_view(infile, (row_offset + col_offset)*sizeof(int), MPI_INT, MATRIX_INT_INNER_VECTOR, "native", MPI_INFO_NULL);
		MPI_File_read(infile, matrixBuffer.data(), matrixBuffer.size(),  MPI_INT, &status);
		break;

		case RIGHT:
		MPI_Type_vector(buffer_col_size, buffer_row_size, nodes, MPI_INT, &MATRIX_INT_RIGHT_EDGE_VECTOR);
		MPI_Type_commit(&MATRIX_INT_RIGHT_EDGE_VECTOR);
		MPI_File_set_view(infile, (row_offset + col_offset)*sizeof(int), MPI_INT, MATRIX_INT_RIGHT_EDGE_VECTOR, "native", MPI_INFO_NULL);
		MPI_File_read(infile, matrixBuffer.data(), matrixBuffer.size(),  MPI_INT, &status);
		break;

		case LOWER:
		MPI_Type_vector(buffer_col_size, buffer_row_size, nodes, MPI_INT, &MATRIX_INT_LOWER_EDGE_VECTOR);
		MPI_Type_commit(&MATRIX_INT_LOWER_EDGE_VECTOR);
		MPI_File_set_view(infile, (row_offset + col_offset)*sizeof(int), MPI_INT, MATRIX_INT_LOWER_EDGE_VECTOR, "native", MPI_INFO_NULL);
		MPI_File_read(infile, matrixBuffer.data(), matrixBuffer.size(),  MPI_INT, &status);
		break;

		case CORNER:
		MPI_Type_vector(buffer_col_size, buffer_row_size, nodes, MPI_INT, &MATRIX_INT_LOWER_RIGHT_CORNER_VECTOR);
		MPI_Type_commit(&MATRIX_INT_LOWER_RIGHT_CORNER_VECTOR);
		MPI_File_set_view(infile, (row_offset + col_offset)*sizeof(int), MPI_INT, MATRIX_INT_LOWER_RIGHT_CORNER_VECTOR, "native", MPI_INFO_NULL);
		MPI_File_read(infile, matrixBuffer.data(), matrixBuffer.size(),  MPI_INT, &status);
		break;
	}


	double t1 = MPI_Wtime();

	if(rank==MASTER) {
		std::cout<<"File reading time:"<<(t1-t0)<<std::endl;
	}

	std::vector<int> rowDataBuffer(buffer_row_size);
	std::vector<int> colDataBuffer(buffer_col_size);
	std::vector<int> myRowData(buffer_row_size);
	std::vector<int> myColData(buffer_col_size);

	MPI_Barrier(MPI_COMM_WORLD);

	for(long k=0;k<nodes;k++)
	{
		if(rank==MASTER) {
			std::cout<<"k:"<<k<<std::endl;
		}

		switch(pe_position)
		{
			case INNER:
				if(pe_layout_row_offset == floor(k/buffer_col_size) && (pe_layout_col_offset == floor(k/buffer_row_size))) {
					// make copy of row to send
						memcpy(&rowDataBuffer[0], &matrixBuffer[(k%buffer_col_size)*buffer_row_size], buffer_row_size*sizeof(int));
					// make copy of column to send
					for(long i=0, j=(k%buffer_row_size); i<buffer_col_size; i++, j+=buffer_row_size) {
						colDataBuffer[i] = matrixBuffer[j];
					}

					//std::cout<<"k:"<<k<<" my inner "<<"rank:"<<rank<<std::endl;
					// Broadcast Col data to ROW_COMM
					MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, pe_layout_row_offset, ROW_COMM);
					// Broadcast Row data to COL_COMM
					
					MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, pe_layout_col_offset, COL_COMM);

				} 
				 else if(pe_layout_row_offset == floor(k/buffer_col_size)) {
					// receive ROW from  row_offset(th) RANK.
					//std::cout<<"k:"<<k<<" my row inner "<<"rank:"<<rank<<std::endl;

					MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, pe_layout_row_offset, ROW_COMM);
					// Store my column in ColDataBuffer
					memcpy(&myRowData[0], &matrixBuffer[(k%buffer_col_size)*buffer_row_size], buffer_row_size*sizeof(int));
					// Broadcast to my COL COMM.

					MPI_Bcast(myRowData.data(), myRowData.size(), MPI_INT, col_rank, COL_COMM);
				}
				 else if(pe_layout_col_offset == floor(k/n_by_rootp)) {
					// receive COL from col_offset(th) RANK.
					
					//std::cout<<"k:"<<k<<" my col inner "<<"rank:"<<rank<<std::endl;
					MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, pe_layout_col_offset, COL_COMM);

					// Store my Row in RowDataBuffer	
					memcpy(&rowDataBuffer[0], &matrixBuffer[(k%n_by_rootp)*n_by_rootp], n_by_rootp*sizeof(int));

					// Broadcast my Col data to my ROW COMM.
					for(long i=0, j=(k%n_by_rootp); i<n_by_rootp; i++, j+=n_by_rootp) {
						myColData[i] = matrixBuffer[j];
					}
					MPI_Bcast(myColData.data(), myColData.size(), MPI_INT, row_rank, ROW_COMM);

				}
				else {
					// Receive COL from ROW COMM.
					int pe_id;
					pe_id = k/n_by_rootp;
					pe_id = (pe_id > (root_pe_count-1))?(root_pe_count-1):pe_id;
					//std::cout<<"k:"<<k<<"rank:"<<rank<<" pe_id:"<<pe_id<<std::endl;
					MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, pe_id, ROW_COMM);

					pe_id = k/n_by_rootp;
					pe_id = (pe_id > (root_pe_count-1))?(root_pe_count-1):pe_id;
					//std::cout<<"k:"<<k<<"rank:"<<rank<<" col pe_id:"<<pe_id<<std::endl;
					// Receive ROW from COL COMM.
					MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, pe_id, COL_COMM);
				}			


				break;

			case CORNER:
				if(pe_layout_row_offset == floor((k+n_mod_rootp*n_mod_rootp)/buffer_col_size) && (pe_layout_col_offset == floor((k+n_mod_rootp*n_mod_rootp)/buffer_row_size))) {
						// make copy of row to send
						memcpy(&rowDataBuffer[0], &matrixBuffer[((k+(n_mod_rootp*n_mod_rootp))%buffer_col_size)*buffer_row_size], buffer_row_size*sizeof(int));
						// make copy of column to send
						for(long i=0, j=((k+(n_mod_rootp*n_mod_rootp))%buffer_row_size); i<buffer_col_size; i++, j+=buffer_row_size) {
							colDataBuffer[i] = matrixBuffer[j];
						}

						// Broadcast ROW data to ROW_COMM
					//						std::cout<<"k:"<<k<<" my corner "<<"rank:"<<rank<<std::endl;


						MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, pe_layout_row_offset, ROW_COMM);
											// Broadcast COL data to COL_COMM

						MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, pe_layout_col_offset, COL_COMM);
					} 
					else {
						// Receive COL from ROW COMM.
						int pe_id;
					pe_id = k/n_by_rootp;
						pe_id = (pe_id > (root_pe_count-1))?(root_pe_count-1):pe_id;
						MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, pe_id, ROW_COMM);

					pe_id = k/n_by_rootp;
						pe_id = (pe_id > (root_pe_count-1))?(root_pe_count-1):pe_id;

						// Receive ROW from COL COMM.
						MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, pe_id, COL_COMM);
					}			

				break;

			case RIGHT:
				if(pe_layout_row_offset == floor(k/buffer_col_size)) {
					// receive ROW from  row_offset(th) RANK.

					//												std::cout<<"k:"<<k<<" my row right "<<"rank:"<<rank<<std::endl;

					MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, pe_layout_row_offset, ROW_COMM);
					// Store my column in ColDataBuffer
					memcpy(&myRowData[0], &matrixBuffer[(k%buffer_col_size)*buffer_row_size], buffer_row_size*sizeof(int));
					// Broadcast to my COL COMM.
					MPI_Bcast(myRowData.data(), myRowData.size(), MPI_INT, col_rank, COL_COMM);
				}
				 else if(pe_layout_col_offset == floor((k+n_mod_rootp*n_mod_rootp)/buffer_row_size)) {
					// receive ROW from col_offset(th) RANK.
						//											std::cout<<"k:"<<k<<" my col right "<<"rank:"<<rank<<std::endl;

					
					MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, pe_layout_col_offset, COL_COMM);

					// Store my Row in RowDataBuffer	
					memcpy(&rowDataBuffer[0], &matrixBuffer[((k+(n_mod_rootp*n_mod_rootp))%buffer_col_size)*buffer_row_size], buffer_row_size*sizeof(int));

					// Broadcast my Col data to my ROW COMM.
					for(long i=0, j=(k%buffer_row_size); i<buffer_col_size; i++, j+=buffer_row_size) {
						colDataBuffer[i] = matrixBuffer[j];
					}
					MPI_Bcast(myColData.data(), myColData.size(), MPI_INT, row_rank, ROW_COMM);
					//std::cout<<"send right rank:"<<rank<<std::endl;
				}
				else {
					int pe_id;
					pe_id = k/n_by_rootp;
					pe_id = (pe_id > (root_pe_count-1))?(root_pe_count-1):pe_id;

					MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, pe_id, ROW_COMM);

					pe_id = k/n_by_rootp;
					pe_id = (pe_id > (root_pe_count-1))?(root_pe_count-1):pe_id;

					// Receive ROW from COL COMM.
					MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, pe_id, COL_COMM);
				}			

				break;

			case LOWER:
				if(pe_layout_row_offset == floor((k+n_mod_rootp*n_mod_rootp)/buffer_col_size)) {
						//											std::cout<<"k:"<<k<<" my row lower "<<"rank:"<<rank<<std::endl;


					// receive COL from  row_offset(th) RANK.
					MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, pe_layout_row_offset, ROW_COMM);
					// Store my column in ColDataBuffer
					memcpy(&myRowData[0], &matrixBuffer[(k%buffer_col_size)*buffer_row_size], buffer_row_size*sizeof(int));
					// Broadcast to my COL COMM.
					MPI_Bcast(myRowData.data(), myRowData.size(), MPI_INT, col_rank, COL_COMM);
				}
				 else if(pe_layout_col_offset == floor(k/buffer_row_size)) {
						//											std::cout<<"k:"<<k<<" my col lower "<<"rank:"<<rank<<std::endl;

					// receive Row from col_offset(th) RANK.

					MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, pe_layout_col_offset, COL_COMM);

					// Store my Row in RowDataBuffer	
					memcpy(&rowDataBuffer[0], &matrixBuffer[(k%n_by_rootp)*n_by_rootp], n_by_rootp*sizeof(int));

					// Broadcast my Col data to my ROW COMM.
					for(long i=0, j=(k%buffer_row_size); i<buffer_col_size; i++, j+=buffer_row_size) {
						colDataBuffer[i] = matrixBuffer[j];
					}
					MPI_Bcast(myColData.data(), myColData.size(), MPI_INT, row_rank, ROW_COMM);
				}
				else {
					int pe_id;
					pe_id = k/n_by_rootp;
					pe_id = (pe_id > (root_pe_count-1))?(root_pe_count-1):pe_id;

					MPI_Bcast(colDataBuffer.data(), colDataBuffer.size(), MPI_INT, pe_id, ROW_COMM);

					pe_id = k/n_by_rootp;
					pe_id = (pe_id > (root_pe_count-1))?(root_pe_count-1):pe_id;

					// Receive ROW from COL COMM.
					MPI_Bcast(rowDataBuffer.data(), rowDataBuffer.size(), MPI_INT, pe_id, COL_COMM);
				}			

			break;
		}
		//std::cout<<"rank:"<<rank<<" waiting at barrier\n";
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

//	MPI_File_open(COMM, outfile_path.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &outfile);
//	MPI_File_set_view(outfile, (row_offset + col_offset)*sizeof(int), MPI_INT, MATRIX_INT_VECTOR, "native", MPI_INFO_NULL);
//	if(rank==0) {
//		std::cout<<"matrix size:"<<matrixBuffer.size()<<std::endl;
//	}
//	MPI_File_write(outfile, matrixBuffer.data(), matrixBuffer.size(),  MPI_INT, &status);
//
//	double t3 = MPI_Wtime();
//	if(rank==MASTER) {
//		std::cout<<"File Writing time:"<<t3-t2<<std::endl;
//	}

	MPI_File_close(&infile);
//	MPI_File_close(&outfile);
	MPI_Comm_free(&ROW_COMM);
	MPI_Comm_free(&COL_COMM);
	switch(pe_position) {
		case INNER:
		MPI_Type_free(&MATRIX_INT_INNER_VECTOR);
		break;

		case RIGHT:
		MPI_Type_free(&MATRIX_INT_RIGHT_EDGE_VECTOR);
		break;

		case LOWER:
		MPI_Type_free(&MATRIX_INT_LOWER_EDGE_VECTOR);
		break;

		case CORNER:
		MPI_Type_free(&MATRIX_INT_LOWER_RIGHT_CORNER_VECTOR);
		break;
	}
}
