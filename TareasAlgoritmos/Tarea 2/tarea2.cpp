// Tarea 2: Ignacio Fullerton Infante
// Compilar: mpic++ -std=c++11 tarea2.cpp -o tarea2
// Ejecutar: mpirun -np 3 tarea2


#include <iostream>
#include <fstream>
#include <mpi.h>
#include <cmath> //Se puede utilizar?
using namespace std;

void print_vector(double* vector, int n, int rank, const char* text) {
    printf("\nRank %i, %s:\n", rank, text);
    for (int i = 0; i < n; i++) {
        printf("%f ", vector[i]);
    }
    printf("\n\n");
} // Función extraída de ayudantía.

int main()
{   
    ifstream file;

    int nrows, ncols;
    double *my_matrix, *localVector, *b_k, *b_k_p, *resultado_parcial;
    double tmp, error, vp_p;
    int indiceVector, tamano_vector, firstIndex, localRows;
    double norma, norma_parcial;
    int itr = 0;
    file.open("matrix.txt");

    file >> nrows;
    file >> ncols;

    int b_0[ncols];
    for (int i = 0; i < ncols; i++){
        b_0[i] = 1;
    }

    file.close();

    MPI_Init(NULL,NULL);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    // Guardar la matriz. Código del Profesor.

    file.open("matrix.txt");

    if (file.is_open())
    {
        file >> nrows;
        file >> ncols;
        
        //if (world_rank == 0){
            //cout << "Number of rows: " << nrows << endl;
            //cout << "Number of columns: " << ncols << endl;
        //}
        
        // Partir vector b_0

        tamano_vector = ncols / world_size;
        indiceVector = world_rank * (ncols / world_size);
        if (world_rank == world_size -1){
            tamano_vector += ncols % world_size;
        }
        //cout << "Rank: " << world_rank << "\nTamaño de su vector: " << tamañoVector << endl;
        localVector = new double [tamano_vector];
        for (int n = 0; n < tamano_vector; n++){
            localVector[n] = b_0[n + indiceVector];
            //printf("%f ", localVector[n]);
            //cout << "Rank: " << world_rank << ", localVector[" << n << "] = " << localVector[n] << endl;
        }

        // Dimensiones matriz
        localRows = nrows / world_size;
        firstIndex = world_rank * (nrows / world_size) + 1;
        if (world_rank == world_size - 1){
            localRows += nrows % world_size;
        }
        //cout << "Rank: " << world_rank << ", first index: " << firstIndex << ", local size: " << localRows << endl;

        // Guardado del bloque de la matriz
        int my_firstrow = firstIndex;
        //cout << "Read " << localRows << " rows starting from row " << my_firstrow << endl;

        my_matrix = new double [localRows * ncols];

        for (int i=0; i<(my_firstrow-1)*ncols; i++) {
            file >> tmp;
            // cout << "skipped: " << tmp << endl;
        }   

        //cout << "Store matrix elements" << endl;
        for (int i=0; i<localRows*ncols; i++) {
            file >> my_matrix[i];
            //cout << "Rank" << world_rank << ", " << i << " " << my_matrix[i] << endl;
        }

        file.close();
    }
    else
    {
        cout << "Unable to open file." << endl;
    }
    double inicio = MPI_Wtime();
    while (true) {
        ++itr;
    // Comunicación para obtener el vector completo a utilizar. Código de ayudantía.
        b_k = new double [ncols];
        int recvcounts[world_size];
        int offsets[world_size];
        for (int i = 0; i < world_size; i++)
        {
            recvcounts[i] = nrows / world_size;
            if (i == world_size - 1){
                recvcounts[i] += ncols % world_size;
            }
            offsets[i] = i * (ncols / world_size);
            //if (i > 0)
            //{
            //    offsets[i] += offsets[i-1];
            //}
            //else { offsets[i] = 0; }
        }

        MPI_Allgatherv(localVector, tamano_vector, MPI_DOUBLE, b_k, recvcounts, offsets, MPI_DOUBLE, MPI_COMM_WORLD);
        // int MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
        //                void *recvbuf, const int *recvcounts, const int *displs,
        //                MPI_Datatype recvtype, MPI_Comm comm)
    
        //const char* text4 = "b_k";
        //print_vector(b_k, ncols, world_rank, text4);
    
    
    // MatVec visto en ayudantía
        //printf("Rank %i, empezando local mat vec\n", world_rank);
        b_k_p = new double[tamano_vector];
	    for (int i=0; i<localRows; i++) {
            b_k_p[i] = 0;
            for (int j=0; j<ncols; j++) {
                b_k_p[i] += my_matrix[i * ncols + j] * b_k[j];
                //cout << j << "->" << localVector[j] << endl;
            }
	    }
        //printf("Rank %i, terminó local mat vec\n", world_rank);
        //const char* parcial = "b_k+1 parcial";
        //print_vector(localVector, localRows, world_rank, parcial);
    
    // Falta la norma de la multiplicación.
        norma = 0;
        norma_parcial = 0;
        for (int i = 0; i < localRows; i++) {
            norma_parcial += pow(b_k_p[i], 2);
        }
 
        MPI_Allreduce(&norma_parcial, &norma, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        norma = sqrt(norma);
        for (int i = 0; i < localRows; i++) {
            b_k_p[i] = b_k_p[i] / norma;
        }
        //const char* normado = "b_k+1 parcial Normado";
        //print_vector(localVector, localRows, world_rank, normado);

    // Calcular valor propio actual y error
        
        MPI_Allgatherv(b_k_p, tamano_vector, MPI_DOUBLE, b_k, recvcounts, offsets, MPI_DOUBLE, MPI_COMM_WORLD);
        resultado_parcial = new double[localRows];
    // MatVec A_p * b_k_p
        for (int i=0; i<localRows; i++) {
            resultado_parcial[i] = 0.0;
            for (int j=0; j<ncols; j++) {
                resultado_parcial[i] += my_matrix[i * ncols + j] * b_k[j];
                //cout << j << "->" << resultado_parcial[j] << endl;
            }
	    }
    // (b_k)^T * resultado anterior
        vp_p = 0;
        for (int i = 0; i < localRows; i++) {
            vp_p += b_k_p[i] * resultado_parcial[i];
        }
    // Obtener la suma total
        double valor_propio = 0;
        MPI_Allreduce(&vp_p, &valor_propio, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        error = fabs(10 - valor_propio);
        //cout << "Rank: " << world_rank << ", global sum after allreduce: " << error << endl;
        localVector = b_k_p;
    // Juntar vector b_{k+1} si error < 10^{-5} o 1000 iteraciones -> Se convierte en b_k de la siguiente iteración.
        
        if (error < pow(10, -5) or itr == 1000) {
            if (world_rank == 0) {
                //const char* conf = "b_{k+1} completo";
                //#print_vector(b_k, ncols, world_rank, conf);
                cout << "Iteración: " << itr << "\nValor propio: " << valor_propio << "\nError: " << error << endl;  
            }
            char processor_name[MPI_MAX_PROCESSOR_NAME];
            int name_len;
            MPI_Get_processor_name(processor_name, &name_len);
            cout << "Procesador utilizado: " << processor_name << endl;
            break;
        }
    }
    double fin = MPI_Wtime();
    double t = fin - inicio;
    if (world_rank == 0) {
        cout << "Tiempo total: " << t << endl;
    }

    // Dejar al final del código para liberar memoria. 
    delete[] my_matrix;
    delete[] localVector;
    delete[] b_k, b_k_p;
    delete[] resultado_parcial;


    MPI_Finalize();

    return 0;
}
