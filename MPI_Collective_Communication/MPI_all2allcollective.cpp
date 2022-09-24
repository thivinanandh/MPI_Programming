#include <iostream>
#include <cstring>
#include <cmath>
#include "mpi.h"
#include "MPI_all2allcollective.h" 


using namespace std;

// Variable Naming Convention 
// All caps - HardCoded Value - Common across all processor
// m_variable name - MPI Variables 
// p_variable      - Processor Specific Value
// g_Variable      - Global Values


int main (int argc , char** argv)
{

    // --------------- SETUP PHASE FOR MPI ----------------------------- //
    // Obtain input variables from thr user 
    // 1 -  Message Size
    MPI_Init(&argc, &argv);
    int m_rank;
    int m_size;
    MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&m_size);

    // Create an array in each MPI Processor
    int p_chunkSize          =  stoi(argv[1]);
    int p_totalSizeinBytes   =  ( p_chunkSize   * m_size ) ; // Total Size in Bytes
    int p_sizeofArray        =  p_totalSizeinBytes / sizeof(int);    // Total Size for Declaration
    int numElemperChunk      =  p_chunkSize / sizeof(int);
    int* p_Array_Original    =  new int[p_sizeofArray]();
    int* p_Array_alltoall    =  new int[p_sizeofArray]();
    int* p_Array_t_all2all   =  new int[p_sizeofArray]();
    // cout << " Rank ; "<<m_rank <<  "  size of Array : " << p_sizeofArray <<endl;

    if(!m_rank)
        // cout << "p_totalSizeinBytes : " <<p_totalSizeinBytes<<" p_sizeofArray : "<<p_sizeofArray<<" p_chunkSize : " <<p_chunkSize<<endl;


    // Fill the Values in the Array
    for ( int i = 0  ; i < p_sizeofArray ; i++){
        p_Array_Original[i]   = m_rank*10 + i;
    }
    
    // --------------- SETUP PHASE FOR MPI  - END ----------------------------- //

    
    // Create a copy of the Array and Do the local arranging of the values

    double start_all2all  = MPI_Wtime();
    MPI_Alltoall(p_Array_Original,numElemperChunk,MPI_INT,p_Array_alltoall,numElemperChunk,MPI_INT,MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);
    double end_all2all  = MPI_Wtime();

    
    double start = MPI_Wtime();
    // thivin_MPI_all2all_bruks(p_Array_Original,numElemperChunk,MPI_INT,p_Array_t_all2all);
    thivin_MPI_all2all_pairwise(p_Array_Original,numElemperChunk,MPI_INT,p_Array_t_all2all);
    MPI_Barrier(MPI_COMM_WORLD);
    double end  = MPI_Wtime();


    int p_check = 0; 
    int totalCheck ;

    for ( int i = 0  ; i < numElemperChunk ; i++)
        if ( p_Array_alltoall[i] != p_Array_t_all2all[i])
            p_check += 1;
    
    MPI_Reduce(&p_check,&totalCheck,1,MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);



    // 
    // MPI_Barrier(MPI_COMM_WORLD);
    // // thivin_MPI_print(p_Array_t_all2all,p_sizeofArray);
    // MPI_Barrier(MPI_COMM_WORLD);

    if(!m_rank)
    if (!totalCheck)
    {
        cout << "size : "<<m_size << " chunk : " << p_chunkSize << " Pairwise" <<endl;
        cout << "" << end_all2all - start_all2all <<"\t";
        cout << "" << end - start <<"\t";
        cout << "" << double ((end_all2all - start_all2all)) /double(end - start) <<endl;  
    }

    // Delete the Dynamically Allocated Arrays
    delete[] p_Array_Original;
    delete[] p_Array_t_all2all;


    MPI_Finalize();

    return 0;
}

