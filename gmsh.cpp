#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
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

struct point {
    double x;
    double y;
};

int main(){

    // Abrimos el archivo .vtk en modo lectura
    ifstream gmsh_file;
    gmsh_file.open("Malla_Ejemplo.vtk");

    string line;

    // Leemos cada renglon del archivo
    while(getline(gmsh_file, line)){

        // Llegamos a la lista de coordenadas de todos los nodos
        if ((line.length() >= 6) && (line.substr(0, 6) == "POINTS")){
            vector<string> words = split_sentence(line);
            int n_points = stoi(words[1].substr()); 

            // Creamos un vector de puntos para guardar la informacion
            vector <point> points(n_points);
            for(int i = 0; i < n_points; i++){
                
            }
        }
    }

    

    gmsh_file.close();

    return 0;
}
