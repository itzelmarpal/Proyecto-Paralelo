#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>
#include <math.h>
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
    point operator+(const point& other) const {
        return {x + other.x, y + other.y};
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

int main(){
    vector<point> points;
    vector<trian> triangles;
    vector<vector<int>> neighbours; 
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

    neighbours.resize(n_triangles);


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

    for (int i = 0; i < n_triangles; i++) {
        centroids[i] = centroids[i] + points[triangles[i].p1];
        centroids[i] = centroids[i] + points[triangles[i].p2];
        centroids[i] = centroids[i] + points[triangles[i].p3];
        centroids[i].x /= 3;
        centroids[i].y /= 3;
    }

    // gradients[i]: approximation of grad f evaluated at centroids[i].
    vector<point> gradients(n_triangles);

    // Compute gradient evaluated at centroids
    for (int i = 0; i < n_triangles; i++) {
        double S_xx = 0.0;
        double S_yy = 0.0;
        double S_xy = 0.0;
        double S_xf = 0.0;
        double S_yf = 0.0;

        for (int j : neighbours[i]) {
            double delta_x = centroids[j].x - centroids[i].x;
            double delta_y = centroids[j].y - centroids[i].y;
            double delta_f = f(centroids[j]) - f(centroids[i]);
            S_xx += delta_x * delta_x;
            S_yy += delta_y * delta_y;
            S_xy += delta_x * delta_y;
            S_xf += delta_x * delta_f;
            S_yf += delta_y * delta_f;
        }

        double D = S_xx * S_yy - S_xy * S_xy;
        gradients[i].x = (S_yy * S_xf - S_xy * S_yf) / D;
        gradients[i].y = (S_xx * S_yf - S_xy * S_xf) / D;
    }

    return 0;
}#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>
#include <math.h>
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
    point operator+(const point& other) const {
        return {x + other.x, y + other.y};
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

int main(){
    vector<point> points;
    vector<trian> triangles;
    vector<vector<int>> neighbours; 
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

    neighbours.resize(n_triangles);


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

    for (int i = 0; i < n_triangles; i++) {
        centroids[i] = centroids[i] + points[triangles[i].p1];
        centroids[i] = centroids[i] + points[triangles[i].p2];
        centroids[i] = centroids[i] + points[triangles[i].p3];
        centroids[i].x /= 3;
        centroids[i].y /= 3;
    }

    // gradients[i]: approximation of grad f evaluated at centroids[i].
    vector<point> gradients(n_triangles);

    // Compute gradient evaluated at centroids
    for (int i = 0; i < n_triangles; i++) {
        double S_xx = 0.0;
        double S_yy = 0.0;
        double S_xy = 0.0;
        double S_xf = 0.0;
        double S_yf = 0.0;

        for (int j : neighbours[i]) {
            double delta_x = centroids[j].x - centroids[i].x;
            double delta_y = centroids[j].y - centroids[i].y;
            double delta_f = f(centroids[j]) - f(centroids[i]);
            S_xx += delta_x * delta_x;
            S_yy += delta_y * delta_y;
            S_xy += delta_x * delta_y;
            S_xf += delta_x * delta_f;
            S_yf += delta_y * delta_f;
        }

        double D = S_xx * S_yy - S_xy * S_xy;
        gradients[i].x = (S_yy * S_xf - S_xy * S_yf) / D;
        gradients[i].y = (S_xx * S_yf - S_xy * S_xf) / D;
    }

    return 0;
}
