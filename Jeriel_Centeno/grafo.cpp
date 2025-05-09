#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstring>
using namespace std;

// Une los elementos del vector separados por espacio
string join(const vector<string>& vec) {
    string result;
    for (size_t i = 0; i < vec.size(); ++i) {
        result += vec[i];
        if (i < vec.size() - 1) result += " ";
    }
    return result;
}

int main(int argc, char* argv[]) {
    if (argc < 9) { //validacion
        cerr << "Not enough arguments provided.\n";
        return 1;
    }

    string pr_dir = argv[1];
    string k = argv[2];
    vector<string> x_axis, y_axis, z_axis;

    // Reunir x_axis hasta el "|"
    int i = 3;
    while (i < argc && string(argv[i]) != "|") {
        x_axis.push_back(argv[i]);
        i++;
    }
    i++; // skip "|"

    while (i < argc && string(argv[i]) != "|") {
        y_axis.push_back(argv[i]);
        i++;
    }
    i++; // skip "|"

    while (i < argc && string(argv[i]) != "true" && string(argv[i]) != "false") {
        z_axis.push_back(argv[i]);
        i++;
    }

    string three_d = argv[i++];
    string stackable = argv[i++];
    string output = argv[i++];

    string cmd = "python3 read_csv.py -i \"" + pr_dir + "/" + k + "mers_F1Scores.csv\"";
    cmd += " -x " + join(x_axis);
    cmd += " -y " + join(y_axis);

    if (three_d == "true") {
        cmd += " -z " + join(z_axis);
        cmd += " -d";
        if (stackable == "true") cmd += " -s";
    } else if (stackable == "true") {
        cmd += " -s";
    }

    cmd += " -o " + output;

    cout << "Ejecutando: " << cmd << endl;
    system(cmd.c_str());

    return 0;
}

//################################ SECCION DE JERIEL END #######################