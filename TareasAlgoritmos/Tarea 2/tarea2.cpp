// Tarea 2: Ignacio Fullerton Infante
// Compilar: mpic++ -std=c++11 tarea2.cpp -o tarea2
// Ejecutar: mpirun -np 3 tarea2


#include <iostream>
#include <fstream>
#include <mpi.h>
using namespace std;

int main()
{   // While o for para iterar el código.
    ifstream file;

    int nrows, ncols;
    double *my_matrix;
    double tmp;
    

    file.open("matrix.txt");

    file >> nrows;
    file >> ncols;

    int b_0[ncols];
    for (int i = 0; i < ncols; i++){
        b_0[i] = 1;
    }
    // int b_k[ncols]

    file.close();

    MPI_Init(NULL,NULL);

    int err;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
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
        
        // Partir vector b_0
        int indiceVector, tamañoVector;
        tamañoVector = ncols / world_size;
        indiceVector = world_rank * (ncols / world_size);
        if (world_rank == world_size -1){
            tamañoVector += ncols % world_size;
        }
        int localVector[tamañoVector];
        for (int n = 0; n < tamañoVector; n++){
            localVector[n] = b_0[n + indiceVector];
            printf("%d ", localVector[n]);
            cout << "Rank: " << world_rank << ", localVector[" << n << "] = " << localVector[n] << endl;
        }

        int firstIndex, localRows;
        localRows = nrows / world_size;
        firstIndex = world_rank * (nrows / world_size) + 1;
        if (world_rank == world_size - 1){
            localRows += nrows % world_size;
        }
        cout << "Rank: " << world_rank << ", first index: " << firstIndex << ", local size: " << localRows << endl;

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

        // MatVec visto en ayudantía
        int* b_k = (int*) calloc(tamañoVector, sizeof(int));
        printf("Rank %i, empezando local mat vec\n", world_rank);
	    for (int i=0; i<localRows; i++) {
            for (int j=0; j<tamañoVector; j++) {
                b_k[j] += my_matrix[i * ncols + j] * localVector[j];
                cout << world_rank << "->" << b_k[j] << endl;
            }
	    }
        printf("Rank %i, terminó local mat vec\n", world_rank);

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
