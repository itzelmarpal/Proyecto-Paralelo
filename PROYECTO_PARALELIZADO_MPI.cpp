#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>
#include <math.h>
#include <mpi.h>
using namespace std;

// LECTURA DE DATOS GMSH

// Funcion que separa un renglon en palabras
vector<string> split_sentence(string sen) {

    stringstream ss(sen);
    string word;
    vector<string> words;

    while (ss >> word) {
        words.push_back(word);
    }
    return words;
}

struct point {
    double x;
    double y;
    point& operator+=(const point& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }
};

struct trian{
    long long p1;
    long long p2;
    long long p3;
};

double f(point p) {
    double x = p.x, y = p.y;
    return 20 * (sin(4 * x) * sin(3 * y) + 0.3 * cos(6 * x) * sin(5 * y) + 1.0);
}

point grad_f(point p) {
    double x = p.x, y = p.y;
    point g;
    g.x = 20.0 * (4.0 * cos(4*x) * sin(3*y) - 1.8 * sin(6*x) * sin(5*y));
    g.y = 20.0 * (3.0 * sin(4*x) * cos(3*y) + 1.5 * cos(6*x) * cos(5*y));
    return g;
}

bool areNeighbours(const trian& T1, const trian& T2) {
    long long T1_vertices[3] = {T1.p1, T1.p2, T1.p3};
    long long T2_vertices[3] = {T2.p1, T2.p2, T2.p3};
    int shared_vert = 0;
    for (long long p : T1_vertices) 
        for (long long q : T2_vertices)
            if (p == q) shared_vert++;
    if (shared_vert == 2) 
        return true;
    return false;
}

void point_sum(void* inputBuffer, void* outputBuffer, int* len, MPI_Datatype* datatype) {
    point* input = (point*) inputBuffer;
    point* output = (point*) outputBuffer;
    for (long long i = 0; i < *len; i++)
        output[i] += input[i];
    return;
}

void write_vtk(const string& filename,
               const vector<point>& points,
               const vector<trian>& triangles,
               const vector<double>& real_point_val,
               const vector<point>& grad_analitico,
               const vector<point>& point_grad)
{
    ofstream file(filename);

    int Np = points.size();
    int Nt = triangles.size();

    file << "# vtk DataFile Version 2.0\n";
    file << "Triangular mesh solution\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    file << "POINTS " << Np << " double\n";
    for (int i = 0; i < Np; i++)
        file << points[i].x << " " << points[i].y << " 0.0\n";

    int total_size = Nt * 4;

    file << "CELLS " <<  Nt << " " << total_size << "\n";
    for (auto& t : triangles) {
        file << "3 " << t.p1 << " " << t.p2 << " " << t.p3 << "\n";
    }

    file << "CELL_TYPES " << Nt << "\n";
    for (int i = 0; i < Nt; i++)
        file << "5\n";

    file << "POINT_DATA " << Np << "\n";
    file << "SCALARS FuncionAnalitica double 1\n";
    file << "LOOKUP_TABLE default\n";

    for (int i = 0; i < Np; i++)
        file << scientific << setprecision(10) << real_point_val[i] << "\n";
    file << "VECTORS GradienteAproximado double\n";
    for (int i = 0; i < Np; i++)
        file << scientific << setprecision(10) << point_grad[i].x
        << " " << point_grad[i].y << " 0.0" << "\n";

    vector<point> error_grad(Np);

    for(int i = 0; i < Np; i++){
        point pnt;
        pnt.x = fabs(point_grad[i].x - grad_analitico[i].x);
        pnt.y = fabs(point_grad[i].y - grad_analitico[i].y);
        error_grad.push_back(pnt);
    }
    
    file << "SCALARS Error_Absoluto_X double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < Np; i++) file << error_grad[i].x << "\n";

    file << "SCALARS Error_Absoluto_Y double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < Np; i++)
        file << error_grad[i].y << "\n";
    
    file.close();
}

int main(){
    int n_tasks, rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &n_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Datatype TRIAN;
    MPI_Type_vector(1, 3, 3, MPI_LONG_LONG, &TRIAN);
    MPI_Type_commit(&TRIAN);

    MPI_Datatype POINT;
    MPI_Type_vector(1, 2, 2, MPI_DOUBLE, &POINT);
    MPI_Type_commit(&POINT);

    MPI_Op POINT_SUM;
    MPI_Op_create(&point_sum, 1, &POINT_SUM);

    vector<point> points;
    vector<trian> triangles;
    long long n_points = 0;
    long long n_triangles = 0;
    
    if (rank == 0) {
        long long counter_triangles = 0;
        // Abrimos el archivo .vtk en modo lectura
        ifstream gmsh_file;
        gmsh_file.open("Malla_Ejemplo.vtk");

        string line;

        // Leemos cada renglon del archivo
        while(getline(gmsh_file, line)){

            // Llegamos a la lista de coordenadas de todos los nodos
            if ((line.length() >= 6) && (line.substr(0, 6) == "POINTS")){
                vector<string> words = split_sentence(line);
                n_points = stoi(words[1].substr()); 

                // Ajustamos el tamaño al vector de puntos
                points.resize(n_points);
                for(int i = 0; i < n_points; i++){
                    getline(gmsh_file, line);
                    stringstream ss(line);

                    double val1, val2;
                    if (ss >> val1 >> val2) {
                        points[i].x = val1;
                        points[i].y = val2;
                    }
                }
            }

            if ((line.length() >= 5) && (line.substr(0, 5) == "CELLS")){
                vector<string> words = split_sentence(line);
                n_triangles = stoi(words[1].substr());
                bool triangles_on_the_way = true;

                while(triangles_on_the_way){

                    getline(gmsh_file, line);
                    stringstream ss(line);
                    double n_nodes;
                    ss >> n_nodes;

                    if(n_nodes == 3){
                        n_triangles = n_triangles - counter_triangles;
                        triangles.resize(n_triangles);

                        long long p1, p2, p3;
                        ss >> p1 >> p2 >> p3;
                        trian new_triangle;

                        triangles[0].p1 = p1;
                        triangles[0].p2 = p2;
                        triangles[0].p3 = p3;
                        triangles_on_the_way = false;
                    }
                    else{
                        counter_triangles++;
                    }
                }
                for(int i = 1; i < n_triangles; i++){
                    getline(gmsh_file, line);
                    stringstream ss(line);

                    long long p1, p2, p3, trash;
                    ss >> trash >> p1 >> p2 >> p3;

                    trian new_triangle;
                    triangles[i].p1 = p1;
                    triangles[i].p2 = p2;
                    triangles[i].p3 = p3;
                }
            }
        }
        cout << "\n\tInput .vtk succesfully read." << endl;
        gmsh_file.close();
    }

    // Broadcast to all processes n_points, n_triangles, points 
    // and triangles
    MPI_Bcast(&n_points, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_triangles, 1, MPI_LONG_LONG,  0, MPI_COMM_WORLD);
    if (rank != 0) {
        points.resize(n_points);
        triangles.resize(n_triangles);
    }
    MPI_Bcast(points.data(), n_points, POINT, 0, MPI_COMM_WORLD);
    MPI_Bcast(triangles.data(), n_triangles, TRIAN, 0, MPI_COMM_WORLD);

    // Iniciamos la particion
    long long partition_np = n_points / n_tasks;
    long long partition_nt = n_triangles / n_tasks;

    long long residuo_np = n_points % n_tasks;
    long long residuo_nt = n_triangles % n_tasks;

    long long local_nt  = partition_nt + (rank == n_tasks - 1 ? residuo_nt : 0);
    long long local_np  = partition_np + (rank == n_tasks - 1 ? residuo_np : 0);

    long long triang_it = rank * partition_nt;
    long long points_it = rank * partition_np;

    vector<int> counts_nt(n_tasks), displs_nt(n_tasks);
    vector<int> counts_np(n_tasks), displs_np(n_tasks);

    for (int i = 0; i < n_tasks; i++) {
        counts_nt[i] = partition_nt + (i == n_tasks - 1 ? residuo_nt : 0);
        displs_nt[i] = i * partition_nt;
        counts_np[i] = partition_np + (i == n_tasks - 1 ? residuo_np : 0);
        displs_np[i] = i * partition_np;
    }

    // Start timer
    double start = MPI_Wtime();

    // neighbours[t]: indexes of triangles neighbouring with triangle t.
    // Adjacency lists of dual graph of triangulation.
    // 1 - crear local_neighbours
    vector<vector<long long>> local_neighbours(local_nt);
    

    // 2 - llenar local_neighbours para los triangulos del proceso actual
    // Compute adjacency lists.
    for(long long t = 0; t < local_nt; t++){
        long long global_t = t + triang_it;
        for(long long s = 0; s < n_triangles; s++){
            if (areNeighbours(triangles[global_t], triangles[s]))
                local_neighbours[t].push_back(s);
        }
    }

    // 1 - CREATE local_centroids
    vector<point> local_centroids(local_nt, {0.0, 0.0});

    // 2 - calcular centroides de triangulos correspondientes a ese proceso

    for(long long t = 0; t < local_nt; t++){
        long long global_t = t + triang_it;
        local_centroids[t] += points[triangles[global_t].p1];
        local_centroids[t] += points[triangles[global_t].p2];
        local_centroids[t] += points[triangles[global_t].p3];
        local_centroids[t].x /= 3.0;
        local_centroids[t].y /= 3.0;
    }

    // 3 - pegar todos los local_centroids en un vector centroids
    // FALTA ARREGLAR PARA QUE EL ULTIMO PROCESO NO SALGA DEL VECTOR CENTROIDS
    vector<point> centroids(n_triangles);
    MPI_Allgatherv(local_centroids.data(), local_nt, POINT, centroids.data(), counts_nt.data(),
                   displs_nt.data(), POINT, MPI_COMM_WORLD);

    // ---------- GRADIENTS
    vector<point> local_gradients(local_nt);
    vector<point> gradients;
    for (long long t = 0; t < local_nt; t++) {
        double S_xx = 0.0;
        double S_yy = 0.0;
        double S_xy = 0.0;
        double S_xf = 0.0;
        double S_yf = 0.0;
        // t is a local index. t_global is a global index.
        long long t_global = t + triang_it;
        // s is a global index
        for (long long s : local_neighbours[t]) {
            double delta_x = centroids[s].x - centroids[t_global].x;
            double delta_y = centroids[s].y - centroids[t_global].y;
            double delta_f = f(centroids[s]) - f(centroids[t_global]);
            S_xx += delta_x * delta_x;
            S_yy += delta_y * delta_y;
            S_xy += delta_x * delta_y;
            S_xf += delta_x * delta_f;
            S_yf += delta_y * delta_f;
        } 
        double D = S_xx * S_yy - S_xy * S_xy;
        local_gradients[t].x = (S_yy * S_xf - S_xy * S_yf) / D;
        local_gradients[t].y = (S_xx * S_yf - S_xy * S_xf) / D;
    }
    // Gather local_gradients in process 0
    if (rank == 0)
        gradients.resize(n_triangles);
    MPI_Gatherv(local_gradients.data(), local_nt, POINT,
            rank == 0 ? gradients.data() : nullptr, counts_nt.data(), displs_nt.data(),
            POINT, 0, MPI_COMM_WORLD);
        
    // ---------- ANALYTIC F VALUE
    // local real_point_val
    vector<double> local_rpv(local_np);
    vector<double> real_point_val;

    for (long long i = 0; i < local_np; i++) {
        long long global_i = i + points_it;
        local_rpv[i] = f(points[global_i]);
    }
    if (rank == 0) real_point_val.resize(n_points);
    MPI_Gatherv(local_rpv.data(), local_np, MPI_DOUBLE,
            rank == 0 ? real_point_val.data() : nullptr, counts_np.data(), displs_np.data(),
            MPI_DOUBLE, 0, MPI_COMM_WORLD);


    // ----- INTERPOLATION -----
    // These three use global indexing
    vector<double> local_pval(n_points, 0.0);
    vector<point> local_pgrad(n_points, {0.0, 0.0});
    vector<int> local_trian_cnt(n_points, 0);
    vector<double> point_val;
    vector<point> point_grad;
    vector<int> trian_cnt;
    
    for (long long t = 0; t < local_nt; t++) {  // t: local index
        long long t_global = t + triang_it;
        long long vertices[3] = {   
                                    triangles[t_global].p1, 
                                    triangles[t_global].p2, 
                                    triangles[t_global].p3
                                };
        for (long long p : vertices) {  // p global index
            local_trian_cnt[p]++;
            local_pval[p] += f(centroids[t_global]);
            local_pgrad[p] += local_gradients[t];
        }
    }
    if (rank == 0) {
        point_val.assign(n_points, 0.0);
        point_grad.assign(n_points, {0.0, 0.0});
        trian_cnt.assign(n_points, 0);
    }
    MPI_Reduce(local_trian_cnt.data(), rank == 0 ? trian_cnt.data()  : nullptr,
               n_points, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_pval.data(),      rank == 0 ? point_val.data()  : nullptr,
               n_points, MPI_DOUBLE,  MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_pgrad.data(),     rank == 0 ? point_grad.data() : nullptr,
               n_points, POINT,       POINT_SUM, 0, MPI_COMM_WORLD);

    vector<double> my_pval(local_np);
    vector<point>  my_pgrad(local_np);
    vector<int>    my_cnt(local_np);

    MPI_Scatterv(rank == 0 ? point_val.data()  : nullptr, counts_np.data(), displs_np.data(),
                 MPI_DOUBLE, my_pval.data(),  local_np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(rank == 0 ? point_grad.data() : nullptr, counts_np.data(), displs_np.data(),
                 POINT, my_pgrad.data(), local_np, POINT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(rank == 0 ? trian_cnt.data()  : nullptr, counts_np.data(), displs_np.data(),
                 MPI_INT, my_cnt.data(), local_np, MPI_INT, 0, MPI_COMM_WORLD);

    for (long long i = 0; i < local_np; i++) {
        my_pval[i] /= my_cnt[i];
        my_pgrad[i].x /= my_cnt[i];
        my_pgrad[i].y /= my_cnt[i];
    }

    MPI_Gatherv(my_pval.data(),  local_np, MPI_DOUBLE,
                rank == 0 ? point_val.data()  : nullptr,
                counts_np.data(), displs_np.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(my_pgrad.data(), local_np, POINT,
                rank == 0 ? point_grad.data() : nullptr,
                counts_np.data(), displs_np.data(), POINT, 0, MPI_COMM_WORLD);
    // Stop timer
    double stop = MPI_Wtime();
    // Compute execution time
    double duration = stop - start;

    // ENCERRAMOS TODO EL CÁLCULO FINAL Y ESCRITURA EN EL RANK 0
    if (rank == 0) {

        vector<point> grad_analitico(n_points);
        for (long long i = 0; i < n_points; i++)
            grad_analitico[i] = grad_f(points[i]);

        // Compute max error
        double max_error = 0.0;
        for (long long i = 0; i < n_points; i++) {
            // Usamos fabs() en lugar de abs() para preservar los decimales
            max_error = max(max_error, fabs(point_val[i] - real_point_val[i]));
        }
        
        // Print results
        cout << "\n\tN_POINTS:\t" << n_points << endl;
        cout << "\tN_TRIANGLES:\t" << n_triangles << endl;
        cout << "\tExecution time (microseconds):\t" << duration << endl;
        cout << "\tMax error:\t" << scientific << setprecision(10) << max_error << endl;

        write_vtk("solucion_triangulos.vtk", points, triangles, real_point_val, grad_analitico, point_grad);
    }

    MPI_Op_free(&POINT_SUM);
    MPI_Type_free(&TRIAN);
    MPI_Type_free(&POINT);
    MPI_Finalize();

    return 0;
}
