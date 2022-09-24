#include <iostream>
#include <cstring>
#include <cmath>
#include "mpi.h"


using namespace std;

void thivin_MPI_all2all_bruks(const void *sendBuffer,int numElemperChunk,  MPI_Datatype sendtype , void *recvBuffer)
{
    // get the Parameters
    int m_rank;
    int m_size;
    int p_chunkSize;
    int p_sizeofArray;
    MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&m_size);

    int* p_Array; int* p_Array_temp; int* p_OutputArray;

    double start = MPI_Wtime();
    
    if (sendtype == MPI_INT){
        p_chunkSize     = sizeof(int)*numElemperChunk;
        p_sizeofArray   = m_size * numElemperChunk;
        p_Array         = (int *)sendBuffer;              // Has already Allocated Space
        p_OutputArray    = (int *)recvBuffer;              // Has already Allocated Space 
        p_Array_temp    =  new int[p_sizeofArray]();        // Have to Allocate Space
    }

    else     // Terminate the Program Gracefully 
    {
        cout << " Unsupported DataType in the function : 'thivin_MPI_all2all_bruks' " <<endl;
        MPI_Finalize();
        exit(0);             
    }
    
    // Create a copy of the Array and Do the local arranging of the values   

    double t1 = MPI_Wtime();
    for ( int chunk = 0  ; chunk < m_size ; chunk++)
    {
        int position   =  ( m_rank +  chunk );
        if(position >= m_size)    position -= m_size;
            int start      =  position ; 
            memcpy(p_Array_temp +  chunk*numElemperChunk, p_Array + position*numElemperChunk ,p_chunkSize );
    }
    double t2 = MPI_Wtime();

    // MPI_Barrier(MPI_COMM_WORLD);     // REMOVE - used for testing part only

    // Create Requests for isend and i receive 
    MPI_Request request[m_size*2];
    MPI_Status status[m_size*2];

    int stage = 0;
    int p_totalStages = ceil(log2(m_size));

    double t3 = MPI_Wtime();
    for ( int stage = 0 ; stage < p_totalStages ; stage++ )
    {
        int shift =  1 << stage ;   // 2^stage
        int counter = 0;
        for ( int chunkNo = 0 ;  chunkNo < m_size ; chunkNo++)
        {
            if ( chunkNo &  ( 1 << stage ) )
            {
                int toSendProcessor = m_rank + shift;
                int toRecvProcessor  = m_rank - shift;
                if (toSendProcessor >= m_size)      
                    toSendProcessor   -= m_size ; 
                if(toRecvProcessor < 0 )           
                    toRecvProcessor   +=  m_size ;
                
                MPI_Irecv(p_Array_temp + chunkNo*numElemperChunk,numElemperChunk,MPI_INT,toRecvProcessor,1000+m_rank,MPI_COMM_WORLD,
                    &request[counter+1]);

                MPI_Isend(p_Array_temp + chunkNo*numElemperChunk,p_chunkSize/sizeof(int),MPI_INT,toSendProcessor,1000+toSendProcessor,MPI_COMM_WORLD, &request[counter]);             
                // MPI_Recv(p_Array_temp + chunkNo*numElemperChunk, p_chunkSize/sizeof(int), MPI_INT, toRecvProcessor,1000 + m_rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                counter = counter + 2 ;
            }
        }
        MPI_Waitall(counter,request,status);
    }
    double t4 = MPI_Wtime();


    double t5 = MPI_Wtime();
    // Reverse Copying
    for ( int chunk = 0  ; chunk < m_size ; chunk++)
    {
        int position   =  ( m_rank -  chunk );
        if(position < 0)    position += m_size;
            int start      =  position ; 
            memcpy(p_OutputArray +  chunk*numElemperChunk, p_Array_temp + position*numElemperChunk ,p_chunkSize );
    }
    double t6 = MPI_Wtime();
    double end  = MPI_Wtime();

    double totalTime = end- start;


    // Delete the dynamically allocated arrays
    delete[] p_Array_temp;
    
    if(!m_rank){
        cout << "       Forward rotate : "<< t2 - t1 << "   " << ((t2 - t1)/ totalTime ) * 100 << " % " <<endl;
        cout << "       Calculation    : "<< t4 - t3 <<"   " << ((t4 - t3) / totalTime)  * 100 << " % " <<endl;
        cout << "       Reverse Rotate    : "<< t6 - t5<<"   " << ((t6 - t5) / totalTime) * 100 << " % " <<  endl;
    }
}


void thivin_MPI_all2all_pairwise(const void *sendBuffer,int numElemperChunk,  MPI_Datatype sendtype , void *recvBuffer)
{
    int m_rank;
    int m_size;
    int p_chunkSize;
    int p_sizeofArray;
    bool powerofTwo = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&m_size);

    int* p_Array; int* p_Array_temp; int* p_OutputArray;

    if (sendtype == MPI_INT){
        p_chunkSize     = sizeof(int)*numElemperChunk;
        p_sizeofArray   = m_size * numElemperChunk;
        p_Array         = (int *)sendBuffer;              // Has already Allocated Space
        p_OutputArray    = (int *)recvBuffer;              // Has already Allocated Space 
        p_Array_temp    =  new int[p_sizeofArray]();        // Have to Allocate Space
    }


    memcpy(p_OutputArray,p_Array,m_size*p_chunkSize);

    if(ceil(log2(m_size)) == floor(log2(m_size)) )  powerofTwo = 1;

    switch(powerofTwo)
    {
        case 1:
        {
             // Create Requests for isend and i receive 
            int p_totalStages = m_size - 1;
            MPI_Request request[2*p_totalStages];
            MPI_Status status[2*p_totalStages];
               
            for ( int stage = 0 ; stage < p_totalStages ; stage++ )
            {
                int toSendProcessor   = m_rank ^ (stage + 1);
                int revcFromProcessor = m_rank ^ (stage + 1);

                MPI_Irecv(p_OutputArray + toSendProcessor*numElemperChunk,numElemperChunk,MPI_INT,revcFromProcessor,1000+m_rank,MPI_COMM_WORLD,
                    &request[0 + stage*2]);
                MPI_Isend(p_Array + toSendProcessor*numElemperChunk,numElemperChunk,MPI_INT,
                            toSendProcessor,1000+toSendProcessor,MPI_COMM_WORLD, &request[1+  stage*2]);   
            }
            MPI_Waitall(2*p_totalStages,request,status);
            break;
        }
        

        case 0:
        {
            if(!m_rank) cout << "NON POWER TWO CASE " <<endl;
            // Create Requests for isend and i receive 
            int p_totalStages = m_size - 1;
            MPI_Request request[2*p_totalStages];
            MPI_Status status[2*p_totalStages];
                
            for ( int stage = 0 ; stage < p_totalStages ; stage++ )
            {
                int toSendProcessor             = m_rank + (stage + 1);
                int revcFromProcessor           = m_rank - (stage + 1);
                if (toSendProcessor >= m_size)          toSendProcessor -= m_size ; 
                if(revcFromProcessor < 0 )              revcFromProcessor        +=  m_size ;

                MPI_Irecv(p_OutputArray + revcFromProcessor*numElemperChunk,numElemperChunk,MPI_INT,revcFromProcessor,1000+m_rank,MPI_COMM_WORLD,
                    &request[0 + stage*2]);
                MPI_Isend(p_Array + toSendProcessor*numElemperChunk,numElemperChunk,MPI_INT,
                            toSendProcessor,1000+toSendProcessor,MPI_COMM_WORLD, &request[1+  stage*2]);   
            }

            MPI_Waitall(2*p_totalStages,request,status);
            break;

            break;
        }

        default:
        {
            cout << " This cannot happen "<<endl;
            MPI_Finalize();
            exit(0);
        }

    }

}


void thivin_MPI_print(int* p_Array_t_all2all, int  p_sizeofArray)
{
    int m_rank;
    int m_size;
    MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&m_size);
    for ( int proc = 0 ; proc < m_size ; proc ++ )
    {
        MPI_Barrier(MPI_COMM_WORLD);
       if(m_rank == proc)
       {
            // cout << " rank : " << m_rank <<endl;
            for ( int i = 0 ; i < p_sizeofArray ; i++)
            {
                cout << p_Array_t_all2all[i] << " \t" ;
            }
            cout <<endl;
       }
    }
}
