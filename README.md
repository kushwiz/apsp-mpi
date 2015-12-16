APSP MPI: All pairs shortest path algorithm using MPI
=======

**Goal:**

To parallelize the Floyd Warshall algorithm using MPI. The sequential version of the algorithm runs in O(n^3). Using 'p' processors, the running time reduces to O(n^3/p + n*log p) considering communication overheads.

**Implementation:**

 Here, we have considered the input graph to be unweighted and undirected. Adjacency matrix was converted to a dense matrix in binary format, which is given as input to the routine. 

To substantiate the correctness and performance of the implementation, tests were run with different combinations (compute nodes x processes). The running time of the algorithm is the factor determining the efficiency of this algorithm. Speedup and Efficiency was calculated for different the combinations and the observations have been analyzed using graphs.

By using a checkerboard version of Floyd-Warshall we visualize the placement of processors in the form of a chessboard and  the dense matrix also gets divided in the similar way. Given n nodes and p processors , each processor performs computations on n/√p * n/√p data.

This algorithm handles all the conditions where the data distribution is uneven ie., when √p doesn't equally divide n. In this case, our algorithm's communication pattern works differently in the different regions of the input dense matrix data.

Parallel version of this algorithm takes it's roots from a publication by Vipin Kumar and Vineet Singh (1991), "Scalability of parallel algorithms for All pairs shortest path problem".

**Experimental Section :**

Speedup and efficiency were calculated using the sequential and parallel running times. However for larger data sets where sequential running time is exponentially larger relative speedup was calculated.

Some interesting observations:
* Broadcasts happen in an exponentially lower order of time. On investigation we were bound to believe that it was because of the  extremely robust and fast IB(Infini-Band) network used in the data center.
* Broadcasts on different communicators appear to happen in parallel. We simulated this scenario using a program to broadcast 100,000 integers in row and column communicators simultaneously.


**Experimental Results:**

Below are the experimental results obtained for different sizes of input graph data.

The below graph shows the speedup of the algorithm for different number of cores using a 4000x4000 node input dense matrix.
![](graphs/4000node.png)

The below graph shows the speedup of the algorithm for different number of cores using a 6000x6000 node input dense matrix.
![](graphs/6000node.png)

The below graph shows the speedup of the algorithm for different number of cores using a 10000x10000 node input dense matrix.
![](graphs/10000node.png)





