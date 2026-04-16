#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
using namespace std;

// LECTURA DE DATOS GMSH

vector<string> split_sentence(string sen) {
  
    // Create a stringstream object
    stringstream ss(sen);
    
    // Variable to hold each word
    string word;
    
    // Vector to store the words
    vector<string> words;
    
    // Extract words from the sentence
    while (ss >> word) {
      
        // Add the word to the vector
        words.push_back(word);
    }
    
    return words;
}

int main(){

    // Abrimos el archivo .vtk en modo lectura
    ifstream gmsh_file;
    gmsh_file.open("Malla_Ejemplo.vtk");

    string line;
    while(getline(gmsh_file, line)){
        if ((line.length() >= 6) && (line.substr(0, 6) == "POINTS")){
            vector<string> words = split_sentence(line);
            int n_points = stoi(words[1].substr()); 
            cout << n_points << endl;
        }
    }

    

    gmsh_file.close();

    return 0;
}