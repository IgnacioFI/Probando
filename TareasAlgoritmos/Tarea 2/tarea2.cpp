// Tarea 2: Ignacio Fullerton Infante
// Compilar: mpic++ -std=c++11 tarea2.cpp -o tarea2
// Ejecutar: mpirun -np 3 tarea2


#include <iostream>
#include <fstream>
#include <mpi.h>
using namespace std;

int main()
{
    int nrows, ncols;
    double *my_matrix;
    double tmp;

    ifstream file;

    file.open("matrix.txt");

    int* iter = NULL;
    if (file.is_open()){
        file >> nrows;
        file >> ncols;
        int* iter = (int*) calloc(ncols, sizeof(int));
        for (int i=0; i<ncols; i++) {
            iter[i] = 1;
        }
        
    }
    else{
        cout << "Unable to open file." << endl;
    }

    file.close();
    
    MPI_Init(NULL,NULL);

    int err;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int localResults;

    // ifstream file;
    
    // Guardar la matriz. Código del Profesor.

    file.open("matrix.txt");

    if (file.is_open())
    {
        file >> nrows;
        file >> ncols;
        
        if (world_rank == 0){
            cout << "Number of rows: " << nrows << endl;
            cout << "Number of columns: " << ncols << endl;
        }
        
        int firstIndex, localRows;
        localRows = nrows / world_size;
        firstIndex = world_rank * (nrows / world_size) + 1;
        if (world_rank == world_size - 1){
            localRows += nrows % world_size;
        }
        cout << "Rank: " << world_rank << ", first index: " << firstIndex << ", local size: " << localRows << endl;

        // Inicializar un vector de 1's
        int localVector[ncols];
        for (int n = 0; n < ncols; n++){
            localVector[n] = iter[n + world_rank];
            printf("%d ", localVector[n]);
            // cout << "Rank: " << world_rank << ", localVector[" << n << "] = " << localVector[n] << endl;
        }
        
        // Guardado del bloque de la matriz
        int my_firstrow = firstIndex;
        cout << "Read " << localRows << " rows starting from row " << my_firstrow << endl;

        my_matrix = new double [localRows * ncols];

        for (int i=0; i<(my_firstrow-1)*ncols; i++) {
            file >> tmp;
            // cout << "skipped: " << tmp << endl;
        }   

        cout << "Store matrix elements" << endl;
        for (int i=0; i<localRows*ncols; i++) {
            file >> my_matrix[i];
            cout << "Rank" << world_rank << ", " << i << " " << my_matrix[i] << endl;
        }

        file.close();
    }
    else
    {
        cout << "Unable to open file." << endl;
    }

    // Dejar al final del código para liberar memoria. 
    delete[] my_matrix;
    


// Paralelizar cálculos.

    MPI_Finalize();
    return 0;
}
