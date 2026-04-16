#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>
#include <math.h>
#include <chrono>
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

int main(){
    vector<point> points;
    vector<trian> triangles;
    vector<vector<int>> neighbours; // HACE FALTA LLENAR ESTA!!
    int n_points = 0;
    int n_triangles = 0;
    int counter_triangles = 0;

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

    // Start timer
    auto start = chrono::steady_clock::now();

    // centroids[i]: centroid of triangle i.
    vector<point> centroids(n_triangles, {0.0, 0.0});

    for (int t = 0; t < n_triangles; t++) {
        centroids[t] += points[triangles[t].p1];
        centroids[t] += points[triangles[t].p2];
        centroids[t] += points[triangles[t].p3];
        centroids[t].x /= 3;
        centroids[t].y /= 3;
    }

    // gradients[i]: approximation of grad f evaluated at centroids[i].
    vector<point> gradients(n_triangles);

    // Compute gradient evaluated at centroids
    for (int t = 0; t < n_triangles; t++) {
        double S_xx = 0.0;
        double S_yy = 0.0;
        double S_xy = 0.0;
        double S_xf = 0.0;
        double S_yf = 0.0;

        for (int s : neighbours[t]) {
            double delta_x = centroids[s].x - centroids[t].x;
            double delta_y = centroids[s].y - centroids[t].y;
            double delta_f = f(centroids[s]) - f(centroids[t]);
            S_xx += delta_x * delta_x;
            S_yy += delta_y * delta_y;
            S_xy += delta_x * delta_y;
            S_xf += delta_x * delta_f;
            S_yf += delta_y * delta_f;
        }

        double D = S_xx * S_yy - S_xy * S_xy;
        gradients[t].x = (S_yy * S_xf - S_xy * S_yf) / D;
        gradients[t].y = (S_xx * S_yf - S_xy * S_xf) / D;
    }

    // ----- INTERPOLATION -----

    // point_val[i]: mean value of f(centroids[j]) for all triangles j incident
    // on vertex i.
    vector<double> point_val(n_points, 0);
    // point_grad[i]: mean value of gradients[j] for all triangles j incident
    // on vertex i.
    vector<point> point_grad(n_points, {0.0, 0.0});
    // incident_triangles_cnt[i]: # of triangles incident on vertex i.
    vector<int> incident_triangles_cnt(n_points, 0);
    
    for (int t = 0; t < n_triangles; t++) {
        int triangle_vertices[3] = {triangles[t].p1, triangles[t].p2, triangles[t].p3};
        for (int p : triangle_vertices) {
            incident_triangles_cnt[p]++;
            point_val[p] += f(centroids[t]);
            point_grad[p] += gradients[t];
        }
    }

    for (int i = 0; i < n_points; i++) {
        point_val[i] /= incident_triangles_cnt[i];
        point_grad[i].x /= incident_triangles_cnt[i];
        point_grad[i].y /= incident_triangles_cnt[i];
    }

    // real_point_val[i]: stores f(points[i]), the analytic value of f
    // evaluated at point i.
    vector<double> real_point_val(n_points);
    for (int i = 0; i < n_points; i++)
        real_point_val[i] = f(points[i]);

    // Stop timer
    auto stop = chrono::steady_clock::now();
    // Compute execution time
    auto duration = chrono::duration_cast<std::chrono::microseconds>(stop - start);

    // Compute max error
    double max_error = 0.0;
    for (int i = 0; i < n_points; i++) 
        max_error = max(max_error, abs(point_val[i] - real_point_val[i]));
    
    // Print results
    cout << "\n\tN_POINTS:\t" << n_points << endl;
    cout << "\tN_TRIANGLES:\t" << n_triangles << endl;
    cout << "\tExecution time (microseconds):\t" << duration.count() << endl;
    cout << "\tMax error:\t" << max_error << endl;

    return 0;
}
