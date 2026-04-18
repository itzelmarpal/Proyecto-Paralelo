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
    return 20 * (sin(4 * x) * sin(3 * y) + 0.3 * cos(6 * x) * sin(5 * y) + 1);
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

int main(){
    int n_tasks, rank;
    MPI_INIT(NULL, NULL);
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

                        double p1, p2, p3;
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

                    double p1, p2, p3, trash;
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

        // Broadcast to all processes n_points, n_triangles, points 
        // and triangles
        MPI_Bcast(&n_points, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&n_triangles, 1, MPI_LONG_LONG,  0, MPI_COMM_WORLD);
        MPI_Bcast(points.data(), n_points, POINT, 0, MPI_COMM_WORLD);
        MPI_Bcast(triangles.data(), n_triangles, TRIAN, 0, MPI_COMM_WORLD);
    }

    // Start timer
    double start = MPI_Wtime();

        // neighbours[t]: indexes of triangles neighbouring with triangle t.
    // Adjacency lists of dual graph of triangulation.
    
    local_np = n_points / n_tasks;
    local_nt = n_triangles / n_tasks;

    // 1 - crear local_neighbours
    vector<vector<long long>> local_neighbours(local_nt);
    long long triang_it = rank * local_nt;

    // 2 - llenar local_neighbours para los triangulos del proceso actual
    // Compute adjacency lists.
    for(long long t = triang_it; t < triang_it + local_nt; t++){
        for(long long s = 0; s < n_triangles; s++){
            if (areNeighbours(triangles[t], triangles[s]))
                local_neighbours[triang_it - t].push_back(s);
        }
    }

    // 1 - CREATE local_centroids
    vector<point> local_centroids(local_nt, {0.0, 0.0});

    // 2 - calcular centroides de triangulos correspondientes a ese proceso

    for(long long t = 0; t < local_nt; t++){
        if(triang_it + t >= n_triangles)
            break;
        local_centroids[t] += points[triangles[triang_it + t].p1];
        local_centroids[t] += points[triangles[triang_it + t].p2];
        local_centroids[t] += points[triangles[triang_it + t].p3];
        local_centroids[t].x /= 3;
        local_centroids[t].y /= 3;
    }

    // 3 - pegar todos los local_centroids en un vector centroids
    // FALTA ARREGLAR PARA QUE EL ULTIMO PROCESO NO SALGA DEL VECTOR CENTROIDS
    vector<point> centroids;
    if (rank == 0)
        centroids.resize(n_triangles);
    if (rank != 0) {
        MPI_Send(local_centroids.data(), local_nt, MPI_POINT, 0, 0, MPI_COMM_WORLD);
    } 
    else {
        for (long long t = 0; t < local_nt; t++)
            centroids[triang_it + t] = local_centroids[t];
        for (int p = 1; p < n_tasks; p++)
            MPI_Recv(&centroids[local_nt * rank], local_nt, MPI_POINT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Bcast(&centroids, n_triangles, POINT, 0, MPI_COMM_WORLD);
    }

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
        long long t_global = t + local_nt * rank;
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
    if (rank == 0) {
        gradients.resize(n_triangles);
        MPI_Gather(local_gradients.data(), local_nt, POINT, gradients.data(), 
                    local_nt, POINT, 0, MPI_COMM_WORLD);
        // Then we write in .vtk file
    } 
        
    // ---------- ANALYTIC F VALUE
    // local real_point_val
    vector<double> local_rpv(local_np);
    vector<double> real_point_val;
    for (long long i = 0; i < local_nt; i++) {
        long long global_i = i + rank * local_np;
        local_rpv[i] = f(points[global_i]);
    }
    // Gather local_rpv's in process 0
    if (rank == 0) {
        real_point_val.resize(n_points);
        MPI_Gather(local_rpv.data(), local_np, MPI_DOUBLE, 
                    real_point_val.data(), local_np, MPI_DOUBLE, 
                    0, MPI_COMM_WORLD);
        // again... we write this in .vtk file
    }

    // ----- INTERPOLATION -----
    // These three use global indexing
    vector<double> local_pval(n_points, 0.0);
    vector<point> local_pgrad(n_points, {0.0, 0.0});
    vector<int> local_trian_cnt(n_points, 0);
    vector<double> point_val;
    vector<point> point_grad;
    vector<int> trian_cnt;
    for (long long t = 0; t < local_nt; t++) {  // t: local index
        long long t_global = t + rank * local_nt;
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
        MPI_Reduce(local_trian_cnt.data(), trian_cnt.data(), n_points, MPI_INT, 
                    MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(local_pval.data(), point_val.data(), n_points, MPI_DOUBLE, 
                    MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(local_pgrad.data(), point_grad.data(), n_points, POINT, 
                    POINT_SUM, 0, MPI_COMM_WORLD);
    }
    
    MPI_Scatter(point_val.data(), local_np, MPI_DOUBLE, local_pval.data(), 
                local_np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(point_grad.data(), local_np, POINT, local_pgrad.data(), 
                local_np, POINT, 0, MPI_COMM_WORLD);
    MPI_Scatter(trian_cnt.data(), local_np, MPI_INT, local_trian_cnt.data(),
                local_np, MPI_INT, 0, MPI_COMM_WORLD);
    for (long long i = 0; i < local_np; i++) {
        local_pval[i] /= local_trian_cnt[i];
        local_pgrad[i].x /= local_trian_cnt[i];
        local_pgrad[i].y /= local_trian_cnt[i];
    }
    MPI_Gather(local_pval.data(), local_np, MPI_DOUBLE, point_val.data(), 
                local_np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(local_pgrad.data(), local_np, POINT, point_grad.data(),
                local_np, POINT, 0, MPI_COMM_WORLD);
    MPI_Gather(local_trian_cnt.data(), local_np, MPI_INT, trian_cnt.data(),
                local_np, MPI_INT, 0, MPI_COMM_WORLD);
    // Stop timer
    double stop = MPI_Wtime();
    // Compute execution time
    double duration = stop - start;

    // Compute max error
    double max_error = 0.0;
    for (long long i = 0; i < n_points; i++) 
        max_error = max(max_error, abs(point_val[i] - real_point_val[i]));
    
    // Print results
    cout << "\n\tN_POINTS:\t" << n_points << endl;
    cout << "\tN_TRIANGLES:\t" << n_triangles << endl;
    cout << "\tExecution time (microseconds):\t" << duration.count() << endl;
    cout << "\tMax error:\t" << max_error << endl;

    MPI_Type_free(&TRIAN);
    MPI_Type_free(&POINT);
    MPI_Finalize();

    return 0;
}
