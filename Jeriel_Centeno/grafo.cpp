#include <iostream>
#include <cstdlib>
using namespace std;

int main(int argc, char* argv[]) 
{
    // Se requieren al menos 9 argumentos
    if (argc < 9) 
    {
        cerr << "Not enough arguments provided.\n";
        return 1;
    }

    // Directorio de entrada y nombre base del archivo
    string pr_dir = argv[1];
    string k = argv[2];

    // Variables para ejes x, y, z
    string x_axis, y_axis, z_axis;
    int i = 3;

    // Leer valores para el eje x hasta el delimitador "|"
    while (i < argc && string(argv[i]) != "|") 
    {
        x_axis += argv[i];
        x_axis += " ";
        i++;
    }
    i++; // Saltar "|"

    // Leer valores para el eje y hasta el siguiente delimitador "|"
    while (i < argc && string(argv[i]) != "|") 
    {
        y_axis += argv[i];
        y_axis += " ";
        i++;
    }
    i++; // Saltar "|"

    // Leer valores para el eje z hasta que se encuentre "true" o "false"
    while (i < argc && string(argv[i]) != "true" && string(argv[i]) != "false") 
    {
        z_axis += argv[i];
        z_axis += " ";
        i++;
    }

    // Leer los flags de entrada: 3D, apilable, y nombre del archivo de salida
    string three_d = argv[i++];
    string stackable = argv[i++];
    string output = argv[i++];

    // Construir el comando para ejecutar el script Python
    string cmd = "python3 read_csv.py -i \"" + pr_dir + "/" + k + "mers_F1Scores.csv\"";
    cmd += " -x " + x_axis;
    cmd += " -y " + y_axis;

    // Si es gráfico 3D, agregar eje z y el flag -d (3D)
    if (three_d == "true") 
    {
        cmd += " -z " + z_axis;
        cmd += " -d";
        if (stackable == "true") 
            cmd += " -s"; // stack en gráfico 3D
    } 
    else if (stackable == "true") 
        cmd += " -s"; // stack en gráfico 2D

    // Agregar el nombre de archivo de salida
    cmd += " -o " + output;

    // Mostrar y ejecutar comando
    cout << "Ejecutando: " << cmd << endl;
    system(cmd.c_str());

    return 0;
}
