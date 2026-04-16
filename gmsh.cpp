#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>
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

struct point{
    double x;
    double y;
};

struct trian{
    long long p1;
    long long p2;
    long long p3;
};

int main(){

    vector <point> points;
    vector <trian> triangles;
    int n_points = 0;
    int n_triangles = 0;
    int counter = 0;

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

            for(int i = 0; i < n_triangles; i++){
                getline(gmsh_file, line);
                stringstream ss(line);

                double n_nodes;
                ss >> n_nodes;
                if(n_nodes == 3){
                    n
                    double p1, p2, p3;
                    ss >> p1 >> p2 >> p3;
                    trian new_triangle;
                        new_triangle.p1 = p1;
                        new_triangle.p2 = p2;
                        new_triangle.p3 = p3;
                    triangles.append(new_triangle);
                }
                else{
                    counter++;
                }
            }
        }
    }
    /*
    cout << fixed << setprecision(16);
    for(int i = 0; i < n_points; i++)
        cout << points[i].x << " " << points[i].y << endl;

    gmsh_file.close();
    */
    return 0;
}
