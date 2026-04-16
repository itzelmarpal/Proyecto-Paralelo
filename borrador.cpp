#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>
#include <math.h>
#include <cmath>
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

double df_x(point p){
    double x = p.x, y = p.y;
    return 20.0 * (4.0 * cos(4.0 * x) * sin(3.0 * y) - 1.8 * sin(6.0 * x) * sin(5.0 * y));
}

double df_y(point p){
    double x = p.x, y = p.y;
    return 20.0 * (3.0 * sin(4.0 * x) * cos(3.0 * y) + 1.5 * cos(6.0 * x) * cos(5.0 * y));
}

bool vtxInCommon(trian t1, trian t2){
    int count = 0;
    if((t1.p1 == t2.p1) || (t1.p1 == t2.p2) || (t1.p1 == t2.p3))
        count++;
    if((t1.p2 == t2.p1) || (t1.p2 == t2.p2) || (t1.p2 == t2.p3))
        count++;
    if((t1.p1 == t2.p1) || (t1.p1 == t2.p2) || (t1.p1 == t2.p3))
        count++;
    return (count == 2);
}

void write_vtk(const string& filename,
               const vector<point> points,
               const vector<trian> triangles,
               vector<double> real_point_val,
               vector<point> grad_analitico,
               vector<point> point_grad)
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

    file << "SCALARS Error_X double 1\n";
    file << "LOOKUP_TABLE default\n";

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
    for (int i = 0; i < Np; i++) file << error_grad[i].y << "\n";
    

    file.close();
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
    gmsh_file.close();

    for(int i = 0; i < n_triangles; i++){
        for(int j = i; j < n_triangles; j++){
            if (vtxInCommon(triangles[i], triangles[j])){
                neighbours[i].push_back(j);
                neighbours[j].push_back(i);
            }
        }
    }

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
        long long triangle_vertices[3] = {triangles[t].p1, triangles[t].p2, triangles[t].p3};
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
    
    vector<point> grad_analitico(n_points);
    for(int pnt = 0; pnt < n_points; pnt++){
        grad_analitico[pnt].x = df_x(points[pnt]);
        grad_analitico[pnt].y = df_y(points[pnt]);
    }

    write_vtk("solucion_triangulos.vtk", points, triangles, real_point_val, grad_analitico, point_grad);
    return 0;
}
