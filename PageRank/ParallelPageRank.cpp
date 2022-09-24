#include<iostream>
#include<fstream>
#include<cstring>
#include<sstream>
#include<algorithm>
#include<map>
#include<numeric>
#include "Sparse_Matrix.h"
#include "/opt/intel/impi/2019.6.166/intel64/include/mpi.h"


using namespace std;


struct CSRArray
{
    int CSRrow;
    int CSRcol;
    double CSRVal;
};

bool compareFunction(CSRArray const a,CSRArray const  b)
{

	if (a.CSRcol != b.CSRcol)
		return (a.CSRcol< b.CSRcol);
	else
        return (a.CSRrow < b.CSRrow);
}

bool compareFunction2(CSRArray const a,CSRArray const  b)
{
    if (a.CSRVal != b.CSRVal)
		return (a.CSRVal > b.CSRVal);
	else
        return (a.CSRrow < b.CSRrow);
}


int main()
{
   
    MPI_Init(NULL,NULL);
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    double totalStartTime , totalEndTime ;
    double startTime,endTime;


   // Create the Class for Adjacency matrix
    Sparse_Matrix* adjacencyGraph = new Sparse_Matrix(); // For Smaller Matrix in Each Processor
    const double dampingFactor = 0.85;


    int TotalrowSize = 0;
    int* VerticeTotal;
   
    // Matrix Data Structures in All processors
    int N_rows;
    int N_NNZ;
    Sparse_Matrix* PrimaryAdjacencyGraph; // To be Filled on root Processor only
    double* SolutionGlobal;   // Filled only in Root Processor
    double* SolutionLocal;    // Filled by All processors 
    double* SolutionLocal_old;

    if(rank == 0){
        totalStartTime =  MPI_Wtime();
        startTime =  MPI_Wtime();
    }
    
    if(rank == 0)   // write the File read in Processor one. 
    {
        PrimaryAdjacencyGraph = new Sparse_Matrix();
        // File read from input file 
        ifstream file;
        file.open("/home/thivin/Assignment1/graph.dat");
        cout <<" File Read " << endl;
        if(!file)
        {
            cout << " Cannot Open the file " <<endl;
        }

        std::string line;
        int linecount = 0;
        bool commentEnd = false;
        int N_V, N_E;
        while(!commentEnd )
        {
            getline(file,line);
        stringstream ss(line);
        string temp;
        ss >> temp ;
        linecount ++;
        if(temp == "p"){
                commentEnd = true;
                ss >> temp;
                ss >> N_V ;
                ss >> N_E ;
                PrimaryAdjacencyGraph->row = N_V;
                PrimaryAdjacencyGraph->col = N_V;
                PrimaryAdjacencyGraph->sizeRowPtr = N_V + 1;
                PrimaryAdjacencyGraph->NNZSize = N_E;
            }
        }
        int rowSize = N_V;
        TotalrowSize = rowSize;
        
        cout << " Vertices : " << PrimaryAdjacencyGraph->row <<endl;
        cout << " Edges : " << PrimaryAdjacencyGraph->NNZSize <<endl;

        PrimaryAdjacencyGraph->rowPtr.resize(PrimaryAdjacencyGraph->row + 1, 0);
        PrimaryAdjacencyGraph->colPtr.resize(PrimaryAdjacencyGraph->NNZSize,0);
        PrimaryAdjacencyGraph->values.resize(PrimaryAdjacencyGraph->NNZSize,0);

        PrimaryAdjacencyGraph->rowPtr[0] = 0;

        cout << " before Transpose " <<endl;        
        
        
        // Declare Size of the Struct
        CSRArray* CSRArrays_temp = new CSRArray[PrimaryAdjacencyGraph->NNZSize];

        VerticeTotal = new int[TotalrowSize]();

        int row;int col;double val;
        linecount = 0;
        // Get the Structures Ready
        while(linecount < PrimaryAdjacencyGraph->NNZSize)
        {
            getline(file,line);
            stringstream ss(line);
            string temp;
            ss >> temp;

            ss >>  row;
            row = row-1;
            ss >> col;
            col = col -1;
            ss >> val;

            CSRArrays_temp[linecount].CSRrow = row;
            CSRArrays_temp[linecount].CSRcol = col;
            CSRArrays_temp[linecount].CSRVal = val; 
            VerticeTotal[row] += val;
            linecount++;
        }   
        cout << " before Transpose " <<endl;
        //------------- TRANSPOSE THE MATRIX  -------------------------- //

        // Sort the CSR Array based on the Column Values
        std::sort(CSRArrays_temp,CSRArrays_temp + PrimaryAdjacencyGraph->NNZSize,compareFunction);

        int* RowPtr1 = PrimaryAdjacencyGraph->Get_RowPtr();
        int* ColPtr1 = PrimaryAdjacencyGraph->Get_ColPtr();
        double* Values1 = PrimaryAdjacencyGraph->Get_Values();


        for ( int index  = 0 ; index < PrimaryAdjacencyGraph->NNZSize ; index++)
        {
            ColPtr1[index] = CSRArrays_temp[index].CSRrow;
            Values1[index] = CSRArrays_temp[index].CSRVal;
            RowPtr1[CSRArrays_temp[index].CSRcol + 1] ++ ;
        }

        for (int k = 0 ; k < PrimaryAdjacencyGraph->row; k++)
            RowPtr1[k+1] += RowPtr1[k];
        
            // PrimaryAdjacencyGraph->Display_matrix();
        file.close();

         cout << " FILE READ COMPLETED " <<endl;
    } 

    
    MPI_Bcast(&TotalrowSize,1,MPI_INT,0,MPI_COMM_WORLD);
    if(rank) VerticeTotal = new int[TotalrowSize];
    MPI_Bcast(&VerticeTotal[0],TotalrowSize,MPI_INT,0,MPI_COMM_WORLD);
    
    //Create The Initial Solution Array and Vertice Sum in each Processor
    SolutionGlobal = new double[TotalrowSize]();


    for ( int  i = 0 ; i < TotalrowSize ; i++)
        SolutionGlobal[i] = double(1.0/TotalrowSize);

    // ----------------- Global Variables for each Processor  ------------------------------------ //
    int Proc_num_rows;


    int* GlobalRowPtr;
    int* GlobalColPtr;
    double* GlobalValues;

    int* rowEndEachProcessor  ;
    int* NNZatEachProcessor ;
    int* colStartEachProcessor;
    int* N_rowPointeratEP;
 
    int* N_rowatEachProcessor = new int[size]();    // Required at each processor for MPI_ALLGATHER of Solution array
    int* rowStartEachProcessor = new int[size]();

    // ----------------- Global Variables for each Processor  ------------------------------------ //

    if(rank ==0)
    {

        GlobalRowPtr = PrimaryAdjacencyGraph->rowPtr.data();
        GlobalColPtr = PrimaryAdjacencyGraph->colPtr.data();
        GlobalValues = PrimaryAdjacencyGraph->values.data();
        cout << " RANK ZERO ------- " <<endl;
        // Define the start and end row POinter for each Processor

        rowEndEachProcessor = new int[size]() ;
        NNZatEachProcessor = new int[size]();
        colStartEachProcessor = new int[size]();
        N_rowPointeratEP = new int[size]();


        int NNZSplitperProcessor = PrimaryAdjacencyGraph->NNZSize/size;
        const int NNZSplit = PrimaryAdjacencyGraph->NNZSize/size;
        cout << " NNZ- Split per processor : " << NNZSplitperProcessor <<endl;
        int PrevRowEnd = 0;
        for ( int k = 0 ; k < size - 1 ; k++)
        {
            // cout << " K : "<<k <<endl;
            bool complete = false;
            for ( int i = PrevRowEnd + 1 ; i < PrimaryAdjacencyGraph->sizeRowPtr && !complete; i++)
            {
                int diff = 0;
                // cout << "Row Ptr  = "<< PrimaryAdjacencyGraph->rowPtr[i] <<endl;
                if(PrimaryAdjacencyGraph->rowPtr[i] >= NNZSplitperProcessor)
                {
                    int current = PrimaryAdjacencyGraph->rowPtr[i];
                    int prev    = PrimaryAdjacencyGraph->rowPtr[i-1];
                    // cout << "(curr,prev) = ("<<current<<","<<prev<<") "<<endl;
                    if(fabs(prev - NNZSplitperProcessor) <= fabs(current - NNZSplitperProcessor))
                    {
                        // Current is Correct Value
                        // cout<< "i,prevRowEnd : " << i <<","<<PrevRowEnd<<endl;
                        rowStartEachProcessor[k] = PrevRowEnd ;
                        rowEndEachProcessor[k]   = i-1;
                        NNZatEachProcessor[k] = prev - PrimaryAdjacencyGraph->rowPtr[PrevRowEnd];
                        N_rowatEachProcessor[k] = rowEndEachProcessor[k] - PrevRowEnd;
                        PrevRowEnd = i-1;
                        diff = fabs(current - NNZSplitperProcessor);
                        NNZSplitperProcessor = (NNZSplitperProcessor + NNZSplit) ;
                        // cout << "** PrevRowEnd, rowStart,rowEnd,NNZEach : "<< rowStartEachProcessor[k] <<" , "
                                            // << rowEndEachProcessor[k] <<" , "<<NNZatEachProcessor[k] <<" , " << NNZSplitperProcessor <<endl;
                        complete = true;
                    }
                    else
                    {
                        rowStartEachProcessor[k] = PrevRowEnd;
                        rowEndEachProcessor[k]   = i;
                        NNZatEachProcessor[k] = current - PrimaryAdjacencyGraph->rowPtr[PrevRowEnd];
                        N_rowatEachProcessor[k] = rowEndEachProcessor[k] - PrevRowEnd;
                        PrevRowEnd = i;
                        diff = fabs(current - NNZSplitperProcessor);
                        NNZSplitperProcessor = (NNZSplitperProcessor + NNZSplit) ;
                        // cout << " PrevRowEnd, rowStart,rowEnd,NNZEach : "<< rowStartEachProcessor[k] <<","
                                            // << rowEndEachProcessor[k] <<","<<NNZatEachProcessor[k] <<"," << NNZSplitperProcessor <<endl;
                        complete = true;
                    }                   
                }

                
            }
        }

        // For last Processor 
        rowStartEachProcessor[size-1] = PrevRowEnd;
        rowEndEachProcessor[size - 1] = PrimaryAdjacencyGraph->row;
        NNZatEachProcessor[size - 1] =  PrimaryAdjacencyGraph->NNZSize - PrimaryAdjacencyGraph->rowPtr[PrevRowEnd];
        N_rowatEachProcessor[size -1] =  rowEndEachProcessor[size - 1]   - PrevRowEnd;
  


        for ( int i = 0 ; i < size ; i++){
            colStartEachProcessor[i] =  PrimaryAdjacencyGraph->rowPtr[rowStartEachProcessor[i]];
            N_rowPointeratEP[i] = N_rowatEachProcessor[i] +1; 
        }

        // COUT
        int NNZ_sum = 0;
        int row_sum = 0;
        for ( int k = 0 ; k < size ; k++)
        {
            // cout<<"k, "<<N_rowatEachProcessor[k]<<" , " << NNZatEachProcessor[k]<<endl;
            row_sum += N_rowatEachProcessor[k];
            NNZ_sum += NNZatEachProcessor[k];

        }

        cout << " NNZ SUM : " <<NNZ_sum <<endl;      cout << " Row SUM : " <<row_sum <<endl;

    }

    if(rank == 0)     {
        endTime =  MPI_Wtime();
        cout << " TIME :  File read & Transpose in Processor  0 : " << double(endTime - startTime) << " sec " <<endl;
    }


    if(rank == 0) startTime =  MPI_Wtime();
    MPI_Request null_request =  MPI_REQUEST_NULL;
    // Scatter the Row Start , Row End, N rows and N NNZ from root to all processors
    MPI_Scatter(&NNZatEachProcessor[0],1,MPI_INT,&N_NNZ,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(N_rowatEachProcessor,size,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(rowStartEachProcessor,size,MPI_INT,0,MPI_COMM_WORLD);
    

    // Adjuest the Send Count & send Buffer of the Root Processor to be zero, 
    // Since we dont want to send the same array to be sent to root processor agian


    // Declare sizes of arrays at each Processor
    N_rows =  N_rowatEachProcessor[rank];
    adjacencyGraph->row = N_rows;
    adjacencyGraph->col = TotalrowSize;
    adjacencyGraph->NNZSize = N_NNZ;
    adjacencyGraph->sizeRowPtr = N_rows + 1;
    adjacencyGraph->rowPtr.resize(N_rows + 1);
    adjacencyGraph->colPtr.resize(N_NNZ);
    adjacencyGraph->values.resize(N_NNZ);

    SolutionLocal = new double[N_rows]();              // b in our Ax = b   only has N_rows elements
    SolutionLocal_old =  new double[N_rows]();
    
    MPI_Barrier(MPI_COMM_WORLD);


 
    // ---------- MPI Scatter The matrix to every Processor ------------//
    // -- Send the Row Pointer --

    int tempRow  = N_rowatEachProcessor[rank] + 1;
    MPI_Scatterv(GlobalRowPtr,N_rowPointeratEP,rowStartEachProcessor,MPI_INT,
                    adjacencyGraph->rowPtr.data(), tempRow,MPI_INT,0,MPI_COMM_WORLD);
    // -- Send the Col Pointer

    MPI_Scatterv(GlobalColPtr,NNZatEachProcessor,colStartEachProcessor,MPI_INT,
                    adjacencyGraph->colPtr.data(), N_NNZ,MPI_INT,0,MPI_COMM_WORLD);


    // ---- Scatter the values Array
    MPI_Scatterv(GlobalValues,NNZatEachProcessor,colStartEachProcessor,MPI_DOUBLE,
                    adjacencyGraph->values.data(), N_NNZ,MPI_DOUBLE,0,MPI_COMM_WORLD);

    
    // --- Broadcast the Initial solution Array to entire Processors to Solution_old
    MPI_Ibcast(SolutionGlobal,TotalrowSize,MPI_DOUBLE,0,MPI_COMM_WORLD,&null_request);



    // Adjust The First Value of Row pointer to 0  in every processor except the root,
    if(rank)
    {
        if(adjacencyGraph->rowPtr[0] != 0)
        {
            int displace = adjacencyGraph->rowPtr[0];
            for ( int i = 0 ; i < adjacencyGraph->sizeRowPtr ; i++)
                adjacencyGraph->rowPtr[i] -= displace;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Delete the Arrays in

    if(rank == 0){
        endTime =  MPI_Wtime();
        cout << " TIME :  Matrix Setup in All Processors : " << double(endTime - startTime) << " sec " <<endl;
    }


    if(rank == 0) startTime =  MPI_Wtime();
    //------------------ Calculate Probability Matrix at each processor --------------------------------------- //
    for (int  i = 0 ; i < adjacencyGraph->row ; i++)
    {
        int start = adjacencyGraph->rowPtr[i];
        int end = adjacencyGraph->rowPtr[i+1];

        
        for ( int k = start ; k < end ; k++)
        {
            double VerticeSum = VerticeTotal[adjacencyGraph->colPtr[k]];
            if(VerticeSum)
                adjacencyGraph->values[k] /= VerticeSum;
        }
    }
    if(rank == 0) {
        endTime =  MPI_Wtime();
        cout << " TIME :  transpose  : " << double(endTime - startTime) << " sec " <<endl;
    }

    //----------------  --END-- Calculate Probability Matrix at each processor --------------------------------------- //


    // ---------------- POWER METHOD ITERATION ------------------------------------ //
   
   
   // ---- Local Variables in All the Processors -- . 
    double residual = 1000;  // Entry Criteria
    unsigned int iteration = 0;
    unsigned int rowMatrix = N_rows;
    const int* Rowptr = adjacencyGraph->rowPtr.data();
    const int* ColPtr = adjacencyGraph->colPtr.data();
    const double* Values = adjacencyGraph->values.data();
    double globalResidual = 0;
    MPI_Request request_all_gather;

    if(rank == 0) startTime =  MPI_Wtime();
    while(residual > 1e-9 && iteration < 70000)
    {
        // Matrix Vector Multiplication 
        for (int i = 0 ; i < rowMatrix ; i++)
        {
            int begin =  Rowptr[i];
            int end   = Rowptr[i+1];
            SolutionLocal[i] = 0;
            for ( int k = begin; k < end ; k++){
                SolutionLocal[i] += (Values[k])*SolutionGlobal[ColPtr[k]];
               // cout <<  "row : "<<i << " Val["<<k<<" ] : " <<  Values[k]  << "  colptr[k] : " <<ColPtr[k]<<  " Sol_old : " << SolutionGlobal[ColPtr[k]]<<endl;
            }
            // Apply Damping
            double temp = SolutionLocal[i] ;
            temp =  dampingFactor*(temp) + (1-dampingFactor)/TotalrowSize;
            SolutionLocal[i] = temp;
        }
        
        // residual 
        residual= 0;
        for ( int i = 0 ; i < rowMatrix; i++){
            residual += pow(fabs(SolutionLocal_old[i]-SolutionLocal[i]), 2)/TotalrowSize;
            SolutionLocal_old[i] = SolutionLocal[i];
            // if(PageRank[i] > 0.001 ) cout << PageRank[i] << "     ";
        }

        iteration ++;


        MPI_Allgatherv(SolutionLocal,N_rows,MPI_DOUBLE,SolutionGlobal,N_rowatEachProcessor,rowStartEachProcessor,MPI_DOUBLE,MPI_COMM_WORLD);
        // MPI_Iallgatherv(SolutionLocal,N_rows,MPI_DOUBLE,SolutionGlobal,N_rowatEachProcessor,rowStartEachProcessor,MPI_DOUBLE,MPI_COMM_WORLD,&request_all_gather);
        
        // Broadcsat Residual to all process
        MPI_Allreduce(&residual,&globalResidual,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        residual = sqrt(globalResidual);
    

   }

                // cout <<endl;
    
    if(rank == 0) {
        endTime =  MPI_Wtime();
        cout<< " Iteration --> : " << iteration <<endl;
        cout << " residual --> : "<< residual<<endl;
        cout << " TIME :  Iteration Time   : " << double(endTime - startTime) << " sec " <<endl;
    }

    if(rank == 0) {
        totalEndTime =  MPI_Wtime();
        cout << " TOTALLLLL  TIME :  total Time Taken    : " << double(totalEndTime - totalStartTime) << " sec " <<endl;
    }

    // Claculate the Rank of the Matrix
    if(rank == 0 ){
        CSRArray* result = new CSRArray[TotalrowSize];
        double min_val = 10000000;
        for(int i =0 ; i < TotalrowSize ; i++)
        {   
            double globsol  = SolutionGlobal[i];
            result[i].CSRrow = i+1;
            result[i].CSRVal = globsol;
            if(globsol < min_val && globsol < 1e-4)
                min_val = globsol;
        }
        if(fabs(min_val - 10000000) < 1e-8) min_val == 0;
        cout << " MIN  vale : " << min_val<<endl;
        sort(result ,result+TotalrowSize,compareFunction2);

        
        // Set Low Prob Value to '0' array
        for( int i = 0 ; i < TotalrowSize ; i++)
            if(fabs(result[i].CSRVal - min_val)  < 1e-10)
                result[i].CSRrow = 0 ;



        ofstream outfile;
        outfile.open("ans_par.txt");
        for ( int i = 0 ; i < TotalrowSize; i++ )
            outfile <<i+1 <<"      "  << result[i].CSRrow << endl;

        outfile.close();
    } 

    
    // --------- END -- POWER METHOD ITERATION  --------------------------------------------- //

    



    MPI_Finalize();
   
}