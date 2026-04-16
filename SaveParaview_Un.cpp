//------------------------------------------------------------
// g++ -std=c++17 -O2 SaveParaview_Un.cpp -o save
// ./save
//------------------------------------------------------------

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

//------------------------------------------------------------
// Escribir VTK
//------------------------------------------------------------
void write_vtk(const string& filename,
               const vector<double>& x,
               const vector<double>& y,
               const vector<vector<int>>& cells,
               const vector<double>& U)
{
    ofstream file(filename);

    int Np = x.size();
    int Nc = cells.size();

    file << "# vtk DataFile Version 2.0\n";
    file << "Triangular mesh solution\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    file << "POINTS " << Np << " double\n";
    for (int i = 0; i < Np; i++)
        file << x[i] << " " << y[i] << " 0.0\n";

    int total_size = 0;
    for (auto& c : cells)
        total_size += c.size() + 1;

    file << "CELLS " << Nc << " " << total_size << "\n";
    for (auto& c : cells) {
        file << "3 " << c[0] << " " << c[1] << " " << c[2] << "\n";
    }

    file << "CELL_TYPES " << Nc << "\n";
    for (int i = 0; i < Nc; i++)
        file << "5\n";

    file << "POINT_DATA " << Np << "\n";
    file << "SCALARS U double\n";
    file << "LOOKUP_TABLE default\n";

    for (int i = 0; i < Np; i++)
        file << scientific << setprecision(10) << U[i] << "\n";

    file.close();
}

//------------------------------------------------------------
// MAIN
//------------------------------------------------------------
int main()
{
    const double pi = M_PI;
    //--------------------------------------------------------
    // NODOS 
    //--------------------------------------------------------
    vector<double> x = {0.0, 1.5, 1.5, 0.0, 0.75, 1.5,  
                        0.75, 0.0, 1.059, 0.937, 0.528, 0.421};

    vector<double> y = {0.0, 0.0, 1.5, 1.5,  0.0, 0.75, 
                        1.5, 0.75, 0.440, 0.937, 0.534, 1.078};

    //--------------------------------------------------------
    // TRIÁNGULOS 
    //--------------------------------------------------------
    vector<vector<int>> cells = {
        {7, 0, 10},
        {5, 2, 9},
        {2, 6, 9},
        {0, 4, 10},
        {6, 3, 11},
        {3, 7, 11},
        {4, 1, 8},
        {1, 5, 8},
        {7, 10, 11},
        {8, 5, 9},
        {8, 9, 10},
        {9, 6, 11},
        {10, 9, 11},
        {4, 8, 10}};

    //--------------------------------------------------------
    // SOLUCIÓN
    //--------------------------------------------------------
    int N = x.size();
    vector<double> U(N);

    double t = 1;

    for (int i = 0; i < N; i++) {
        U[i] = exp(-t) * sin(pi*x[i]) * sin(pi*y[i]);
    }

    //--------------------------------------------------------
    // GUARDAR
    //--------------------------------------------------------
    write_vtk("solucion_triangulos.vtk", x, y, cells, U);

    return 0;
}