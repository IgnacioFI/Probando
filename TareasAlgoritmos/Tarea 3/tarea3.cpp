// Tarea 3: Ignacio Fullerton Infante
// Compilar: g++ -o name_output tarea3.cpp -fopenmp
// Ejecutar: ./name_output


#include <vector> // Los gráficos fueron creados por Perplexity
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
using namespace std;

// Función creada con Perplexity
// f(x, y) = exp(-((x-cx)^2/(2*sx^2) + (y-cy)^2/(2*sy^2))) / (2*PI*sx*sy)
float func(float x, float y, float cx, float cy, float sx, float sy) {
    const float PI = 3.14159265358979323846;
    float exponent = -(
        ((x - cx) * (x - cx)) / (2.0 * sx * sx) +
        ((y - cy) * (y - cy)) / (2.0 * sy * sy)
    );
    float denominator = 2.0 * PI * sx * sy;
    return std::exp(exponent) / denominator;
}

float alpha(float x, float y) {
    return x * (x - 1) * y * (y - 1) + 1;
}

float calc_centro(int x, int y, float h_x, float h_y) {
    float primero = (float)(alpha((x - 0.5) * h_x, y * h_y) + alpha((x + 0.5) * h_x, y * h_y)) / (h_x * h_x);
    float segundo = (float)(alpha(x * h_x, (y - 0.5) * h_y) + alpha(x * h_x, (y + 0.5) * h_y)) / (h_y * h_y);
    return primero + segundo + 1;
}

float calc_norte(int x, int y, float h_x, float h_y) {
    return - (float)alpha(x * h_x, (y + 0.5) * h_y) / pow(h_y, 2);
}

float calc_sur(int x, int y, float h_x, float h_y) {
    return - (float)alpha(x * h_x, (y - 0.5) * h_y) / pow(h_y, 2);
}

float calc_este(int x, int y, float h_x, float h_y) {
    return - (float)alpha((x + 0.5) * h_x, y * h_y) / pow(h_x, 2);
}

float calc_oeste(int x, int y, float h_x, float h_y) {
    return - (float)alpha((x - 0.5) * h_x, y * h_y) / pow(h_x, 2);
}

// Función extraída de ayudantía y modificada para la tarea.
float* mat_vec_par(float *norte, float *sur, float *este, float *oeste, float *centro, float *vector, int n, int x, int y)
{
    float *result = (float*)calloc(n, sizeof(float));

    #pragma omp parallel for num_threads(10)
    for (int j = 1; j < y; j++)
    {
        for (int i = 1; i < x; i++)
        {
            int k = j * (x + 1) + i;
            result[k] += vector[k] * centro[k];
            result[k] += vector[k + x + 1] * norte[k];
            result[k] += vector[k - x - 1] * sur[k];
            result[k] += vector[k + 1] * este[k];
            result[k] += vector[k - 1] * oeste[k];
        }
    }
  return result;
}

// Función que imprime en consola un vector. Ayudantía
void print_vector(float* vector, int n)
{
    for (int i = 0; i < n; i++)
    {
        printf("%f\t", vector[i]);
    }
    printf("\n\n");
}

int main(){
    // Definiciones de las variables a utilizar

    int hilos = 0;

    int N_x = 600;
    int N_y = 500;
    int dim = (N_x + 1) * (N_y + 1);
    float h_x = (float)1 / N_x;
    float h_y = (float)1 / N_y;

    float* array_norte = (float*) calloc(dim, sizeof(float));
    float* array_oeste = (float*) calloc(dim, sizeof(float));
    float* array_este = (float*) calloc(dim, sizeof(float));
    float* array_sur = (float*) calloc(dim, sizeof(float));
    float* array_centro = (float*) calloc(dim, sizeof(float));

    float* array_b = (float*) calloc(dim, sizeof(float));
    float* array_x = (float*) calloc(dim, sizeof(float));

    float* array_r = (float*) calloc(dim, sizeof(float));
    float* array_z = (float*) calloc(dim, sizeof(float));
    float* array_p = (float*) calloc(dim, sizeof(float));
    float* array_q = (float*) calloc(dim, sizeof(float));

    float norma_r;
    float ro_0;
    float ro_1 = 0;
    float beta;
    float delta;

    float cx = 0.4;
    float cy = 0.8;
    float sx = 0.2;
    float sy = 0.1;
    
    // Vectores para almacenar el número de iteración y el error/residuo correspondiente
    std::vector<int> iter_history;
    std::vector<double> error_history;


    // Inicializar vector b y matriz A
    #pragma omp parallel for num_threads(hilos)
    for (int j = 1; j < N_y; j++){
        for (int i = 1; i < N_x; i++){
            int k = j * (N_x + 1) + i;
            float x = i * h_x;
            float y = j * h_y;

            array_b[k] += func(x, y, cx, cy, sx, sy);
            array_r[k] -= array_b[k];
            array_centro[k] += calc_centro(i, j, h_x, h_y);
            array_norte[k] += calc_norte(i, j, h_x, h_y);
            array_este[k] += calc_este(i, j, h_x, h_y);
            array_oeste[k] += calc_oeste(i, j, h_x, h_y);
            array_sur[k] += calc_sur(i, j, h_x, h_y);
        }
    }

    for (int iter = 0; iter < 2000; iter++) {
        // cout << "Iteración: " << iter << endl;
        memcpy(array_z, array_r, dim * sizeof(float)); // Función entregada por Perplexity para copiar arreglos.

        // Producto interno paralelizable
        ro_0 = ro_1;
        ro_1 = 0;
        #pragma omp parallel for num_threads(hilos) reduction(+:ro_1) // Extraído de ayudantía
        for (int k = 0; k < dim; k++){
            ro_1 += array_r[k] * array_z[k];
        }

        if (iter == 0) {
            memcpy(array_p, array_z, dim * sizeof(float));
        }
        else {
            beta = ro_1 / ro_0;
            // Iteración paralelizable
            #pragma omp parallel for num_threads(hilos) schedule(static) // Extraído de ayudantía
            for (int k = 0; k < dim; k++) {
                array_p[k] = array_z[k] + beta * array_p[k];
            }
        }
        // cout << "Beta: " << beta <<endl;
        // Matvec paralelo
        array_q = mat_vec_par(array_norte, array_sur, array_este, array_oeste, array_centro, array_p, dim, N_x, N_y);
        
        // Paralelizable
        float denominador = 0;
        #pragma omp parallel for num_threads(hilos) reduction(+:denominador) // Extraído de ayudantía
        for (int k = 0; k < dim; k++) {
            denominador += (array_p[k] * array_q[k]);
        }
        delta = ro_1 / denominador;
        // cout << "Delta: " << delta <<endl;

        // Paralelizable
        #pragma omp parallel for num_threads(hilos) schedule(static) // Extraído de ayudantía
        for (int k = 0; k < dim; k++) {
            if (k % N_x != 0 && k % N_y != 0) {
                array_x[k] -= delta * array_p[k];
            }
        }

        // Paralelizable
        #pragma omp parallel for num_threads(hilos) schedule(static) // Extraído de ayudantía
        for (int k = 0; k < dim; k++) {
            array_r[k] -= delta * array_q[k];
        }

        // Paralelizable
        norma_r = 0;
        #pragma omp parallel for num_threads(hilos) reduction(+:norma_r) // Extraído de ayudantía
        for (int k = 0; k < dim; k++) {
            norma_r += array_r[k] * array_r[k];
        }
        norma_r = sqrt(norma_r);
        // cout << "\nVector x:\n" << endl;
        // print_vector(array_x, dim);
        // cout << "Norma:" << norma_r << endl;


        // Guarda historial
        iter_history.push_back(iter);
        error_history.push_back(norma_r);


        if (norma_r < pow(10, -6)) {
            // cout << "\nVector x:\n" << endl;
            // print_vector(array_x, dim);
            break;
        }
    }
    

    // Al finalizar, exporta los vectores a archivo para graficar (por ejemplo, en formato CSV)
    std::ofstream outfile("hist_convergencia_10.csv");
    if (!outfile) {
    std::cerr << "No se puede abrir el archivo de salida.\n";
    }
    else {
        outfile << "Iteracion,Error\n";
        for (size_t i = 0; i < iter_history.size(); ++i) {
            outfile << iter_history[i] << "," << error_history[i] << "\n";
        }
        outfile.close();
        std::cout << "Guardado el historial en hist_convergencia_10.csv\n";
    }

    cout << "Norma:" << norma_r << endl;

    free(array_centro);
    free(array_este);
    free(array_norte);
    free(array_sur);
    free(array_oeste);

    free(array_b);
    free(array_x);

    free(array_r);
    free(array_z);
    free(array_p);
    free(array_q);

    return 0;
}

