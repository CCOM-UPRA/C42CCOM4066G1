#include <limits.h>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <unistd.h>
#include <cstdio>
#include "includes/cxxopts.hpp"

using std::string, std::vector, std::cout, std::endl;

/*
Precision & Recall

Authors: 
    Juan O. Lopez (juano.lopez@upr. edu)
    Emanuel Martinez (emanuel.martinez8@upr.edu)
    Javier Quinones (javier.quinones3@upr.edu)

License: Creative Commons Attribution-ShareAlike 4.0
http://creativecommons.org/licenses/by-sa/4.0/1egalcode
*/

string get_slurm_dir() {
    const char* slurm_job_id = std::getenv("SLURM_JOB_ID");
    if (slurm_job_id) {
        string command = "scontrol show job " + string(slurm_job_id) + " | awk -F= '/WorkDir=/{print $2}'";
        FILE* pipe = popen(command.c_str(), "r");
        if (!pipe) return "";
        char buffer[256];
        string result;
        if (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            result = buffer;
            result.erase(result.find_last_not_of(" \n\r\t") + 1);
        }
        pclose(pipe);
        return result;
    } else {
        char result[PATH_MAX];
        ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
        if (count != -1) {
            string path(result, count);
            size_t last_slash = path.find_last_of('/');
            if (last_slash != string::npos) {
                return path.substr(0, last_slash);
            }
            return path;
        } else {
            return "";
        }
    }
}

struct PRConfig {
    string genome;
    string fasta;
    string sam;
    string meta;
    int l1count;
    string printdir;
    bool graph;
    bool three_d;
    bool stackable;
    string output;
    vector<string> x_axis;
    vector<string> y_axis;
    vector<string> z_axis;
    int EditStart;
    int EditEnd;
    int EditStep;
    int TresholdStart;
    int TresholdEnd;
    int TresholdStep;
    string DIR;
};

PRConfig run_precision_recall(int argc, char* argv[]) {
    // Define command line options with descriptions matching the original Bash script flags
    cxxopts::Options options("precision_recall", "Precision and Recall computation script");

    options.add_options()
        // -g: Genome file
        ("g,genome", "Genome file.", cxxopts::value<string>())
        // -f: Fasta file
        ("f,fasta", "Fasta file with all the probes.", cxxopts::value<string>())
        // -s: SAM directory
        ("s,sam", "Directory where SAM files will be generated.", cxxopts::value<string>())
        // -c: Meta directory
        ("c,meta", "Directory with LINE-1 meta data file in .csv format.", cxxopts::value<string>())
        // -l: Entries count
        ("l,entries", "Entries for the entire genome.", cxxopts::value<int>())
        // -p: Print directory, default "result"
        ("p,print", "Directory where results and raw files will be saved in.", cxxopts::value<string>()->default_value("result"))
        // -e: Edit distance range, default "5, 20, 5"
        ("e,edit", "Vary the edit distance for MRFAST.", cxxopts::value<string>()->default_value("5,20,5"))
        // -t: Threshold range (note typo 'treshold' kept for exact match), default "25, 500, 25"
        ("t,treshold", "Vary threshold for Precision and Recall.", cxxopts::value<string>()->default_value("25,500,25"))
        // -m: Minimum kmers, default 2
        ("m,MinKmers", "Minimum amount of kmers to be used at once.", cxxopts::value<int>()->default_value("2"))
        // -G: Generate graph, default false
        ("G,graph", "Generates a graph from the results of precision and recall.", cxxopts::value<bool>()->default_value("false"))
        // -x: x-axis headers
        ("x,x_axis", "Specifies the header name for the x axis of the graph.", cxxopts::value<vector<string>>())
        // -y: y-axis headers
        ("y,y_axis", "Specifies the header name for the y axis of the graph.", cxxopts::value<vector<string>>())
        // -z: z-axis headers
        ("z,z_axis", "Specifies the header name for the z axis of the graph.", cxxopts::value<vector<string>>())
        // -d: 3D graph flag, default false
        ("d,three_d", "Graph will be three dimensional.", cxxopts::value<bool>()->default_value("false"))
        // -s: Stackable graphs, default false. 
        // Note: duplicate short flag 's' with sam directory, short flag removed.
        ("stackable", "Graphs will be displayed on the same plain.", cxxopts::value<bool>()->default_value("false"))
        // -o: Output graph filename, default "graph.png".
        ("o,output", "Name of the output graph", cxxopts::value<string>()->default_value("graph.png"))
        // -h: Print help menu.
        ("h,help", "Print usage");

    // Parse command line arguments
    auto result = options.parse(argc, argv);

    // If help requested, print usage and exit program
    if (result.count("help")) {
        cout << options.help() << endl;
        exit(0);
    }
    
    // Create the PRConfig instance that this function will return
    PRConfig config;

    // Extract values from parsed options
    // Required arguments (no default values)
    config.genome = result["genome"].as<string>();
    config.fasta = result["fasta"].as<string>();
    config.sam = result["sam"].as<string>();
    config.meta = result["meta"].as<string>();
    config.l1count = result["entries"].as<int>();

    // Optional arguments with defaults
    config.printdir = result["print"].as<string>();
    config.graph = result["graph"].as<bool>();
    config.three_d = result["three_d"].as<bool>();
    config.stackable = result["stackable"].as<bool>();
    config.output = result["output"].as<string>();

    // Optional vector arguments, default to single headers if not provided
    config.x_axis = result.count("x_axis") ? result["x_axis"].as<vector<string>>() : vector<string>{"Edit dist"};
    config.y_axis = result.count("y_axis") ? result["y_axis"].as<vector<string>>() : vector<string>{"Threshold"};
    config.z_axis = result.count("z_axis") ? result["z_axis"].as<vector<string>>() : vector<string>{"F1 score"};

    // Parse edit distance range string into integers
    string edit_range = result["edit"].as<string>();
    std::sscanf(edit_range.c_str(), "%d,%d,%d", &config.EditStart, &config.EditEnd, &config.EditStep);

    // Parse threshold range string into integers
    string treshold_range = result["treshold"].as<string>(); // 'treshold' typo kept intentionally for exact match with original Bash script
    std::sscanf(treshold_range.c_str(), "%d,%d,%d", &config.TresholdStart, &config.TresholdEnd, &config.TresholdStep);

    // Get SLURM working directory if applicable
    config.DIR = get_slurm_dir();

    return config;
}

void configOptTester(PRConfig config) {
    /* Testing */
    cout << "SLURM or current executable directory: " << config.DIR << endl;
    cout << "Genome file: " << config.genome << endl;
    cout << "Fasta file: " << config.fasta << endl;
    cout << "SAM directory: " << config.sam << endl;
    cout << "Meta directory: " << config.meta << endl;
    cout << "Number of L1 entries: " << config.l1count << endl;
    cout << "Print directory: " << config.printdir << endl;
    cout << "Generate graph: " << (config.graph ? "Yes" : "No") << endl;
    cout << "3D Graph: " << (config.three_d ? "Yes" : "No") << endl;
    cout << "Stackable graphs: " << (config.stackable ? "Yes" : "No") << endl;
    cout << "Output file: " << config.output << endl;
    cout << "Edit distance range: " << config.EditStart << " to " << config.EditEnd << " step " << config.EditStep << endl;
    cout << "Threshold range: " << config.TresholdStart << " to " << config.TresholdEnd << " step " << config.TresholdStep << endl;
    cout << "X axis headers: ";
    for (const auto& header : config.x_axis) cout << header << " ";
        cout << endl;
        cout << "Y axis headers: ";
    for (const auto& header : config.y_axis) cout << header << " ";
        cout << endl;
        cout << "Z axis headers: ";
    for (const auto& header : config.z_axis) cout << header << " ";
}

int main(int argc, char* argv[]) {
    PRConfig config = run_precision_recall(argc, argv);

    // ##########################################
    // # Simplify access to parameter variables #
    // ##########################################

    string GENOME=config.genome;
    string FASTAfile=config.fasta;
    string SAMfile=config.sam;
    string DATA_DIR=config.meta;
    int L1Count=config.l1count;

    configOptTester(config);

    return 0;
}
