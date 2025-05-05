#include <limits.h>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <unistd.h>
#include <cstdio>
#include "includes/cxxopts.hpp"
#include <omp.h>
// Possibly unnecessary
#include <fstream>

using std::string, std::vector, std::cout, std::endl;
// This is the consolidated script so far work is still being done on the individual parts 

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

void executePythonScript(const string script_path, const string args, bool time)
{
    string command = "python3 " + script_path + " " + args;
    if (time)
        command = "time " + command;

    int result = system(command.c_str());
    if (result != 0)
        cerr << "Error executing: " << command << endl;
}

int executeCommand(const string command, string &output)
{
    FILE *p;
    int ch;
    output.clear();

    p = popen(command.c_str(), "r");
    if (p == NULL)
    {
        puts("Unable to open process");
        return 1;
    }

    while ((ch = fgetc(p)) != EOF)
        output += ch;

    return pclose(p);
}



int main(int argc, char* argv[])
 {
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

    // If help requested, print usage and exit
    if (result.count("help")) {
        cout << options.help() << endl;
        return 0;
    }

    // Extract values from parsed options
    // Required arguments (no default values)
    string genome = result["genome"].as<string>();
    string fasta = result["fasta"].as<string>();
    string sam = result["sam"].as<string>();
    string meta = result["meta"].as<string>();
    int l1count = result["entries"].as<int>();

    // Optional arguments with defaults
    string printdir = result["print"].as<string>();
    bool graph = result["graph"].as<bool>();
    bool three_d = result["three_d"].as<bool>();
    bool stackable = result["stackable"].as<bool>();
    string output = result["output"].as<string>();

    // Optional vector arguments, default to single headers if not provided
    vector<string> x_axis = result.count("x_axis") ? result["x_axis"].as<vector<string>>() : vector<string>{"Edit dist"};
    vector<string> y_axis = result.count("y_axis") ? result["y_axis"].as<vector<string>>() : vector<string>{"Threshold"};
    vector<string> z_axis = result.count("z_axis") ? result["z_axis"].as<vector<string>>() : vector<string>{"F1 score"};

    // Parse edit distance range string into integers
    string edit_range = result["edit"].as<string>();
    int EditStart, EditEnd, EditStep;
    std::sscanf(edit_range.c_str(), "%d,%d,%d", &EditStart, &EditEnd, &EditStep);

    // Parse threshold range string into integers
    string treshold_range = result["treshold"].as<string>(); // 'treshold' typo kept intentionally for exact match with original Bash script
    int TresholdStart, TresholdEnd, TresholdStep;
    std::sscanf(treshold_range.c_str(), "%d,%d,%d", &TresholdStart, &TresholdEnd, &TresholdStep);

    // Get SLURM working directory if applicable
    string DIR = get_slurm_dir();

    // Simulated outer loop to continue consolidation 
    for (e = $EditStart; e <= $EditEnd; e += $EditStep)
    {
        // Seccion de Giovanni

        int runCommand(const string command, string &output)

        FILE *p;
        int ch;

        p = popen(command.c_str(), "r");
        if (p == NULL)
        {
            puts("Unable to open process");
            return 1;
        }

        while ((ch = fgetc(p)) != EOF)
            output += ch;

        return pclose(p);


        // Using predetermined threshold values to see the effect on the P&R
        for(int t=ThresholdStart; t<ThresholdEnd; t+=ThresholdStep)
        {
            // Generate the CSV file with the established values
            // A file is created with "ofstream outFile" using the values from the for loops
            ofstream outFile(pr_dir+"/results_e"+to_string(e)+"_t"+to_string(t)+".csv");

            // Checking if the file was created successfully
            if(outFile.is_open())
            {
                // Writing the header of the file
                outFile<<"K-mers min\tPatterns\tTrue +\tFalse +\tFalse -\tPrecision\tRecall\tF1 Score\n";

        // #######################################################################################################
		// #                                       MinKmer Runtime Reducer                                       #
		// # As a specific m value provides best F1 Score in the mayority of cases with any threshold in a       #
		// # given edit distance, do a test run using whole range of m values (50% to 100% of total probes)      #
		// # for each edit distance with first treshold value. After seeing which m value gave highest F1 value  #
		// # for current e value, repeat first and next treshold values with one m value for each edit distance. #
		// #######################################################################################################
        for (int m = MinKmers; m <= MaxKmers; m += 1)
        {
            string strE = to_string(e);
            string strT = to_string(t);
            string strM = to_string(m);
    
            cout << "edit " << strE << ", threshold " << strT << ", MinKmers " << strM << endl;
    
            string rawfilename = "L1s_raw_e" + strE + "_t" + strT + "_m" + strM + ".csv";
            string filteredfilename = "L1s_filtered_e" + strE + "_t" + strT + "_m" + strM + ".csv";
    
            // Generate a CSV with possible L1s
            executePythonScript(DIR + "/L1PD_files/L1PD.py", SAMfile + " " + FASTAfile + " -t " + strT + " -m " + strM + " --data_dir " + DATA_DIR + " --csvoutput > " + pr_dir + "/" + rawfilename, true);
    
            // Compare possible L1s from L1PD against L1s from L1Base2 to filter L1s with no matches
            executePythonScript(DIR + "/precision_recall_files/filter_possible_L1s_CSV.py", pr_dir + "/" + rawfilename + " " + FASTAfile + " -t " + strT + " --data_dir " + DATA_DIR + " > " + pr_dir + "/" + filteredfilename, true);
    
            // Use wc t get line count, with sed's help to remove trailing file name
            // Line counts correspond to Positives and True Positives.
            // cambiar por grep (maybe)
            string output = "";
            executeCommand("wc -l " + pr_dir + "/" + rawfilename + " | sed 's/ .*//'", output);
            int PatCount = stoi(output);
    
            if (PatCount > 0)
            {
                // True positives, false positives, and false negatives
                executeCommand("grep -v '^DUPLICATE' " + pr_dir + "/" + filteredfilename + " | wc -l | sed 's/ .*//'", output);
    
                int TPCount = stoi(output);
                int FPCount = PatCount - TPCount;
                int FNCount = L1Count - TPCount;
    
                // We assume the file with full-length intact L1s has 'FLI-L1'
                // in its name (upper or lower case).
                executeCommand("grep -Fci 'FLI-L1' " + pr_dir + "/" + filteredfilename, output);
                int FLIL1sFound = stoi(output);
    
                // Calculate precision, recall, and F1 Score.
                float precision = TPCount / PatCount;
                float recall = TPCount / L1Count;
                float F1Score = 2 * precision * recall / (precision + recall);
    
                string strPatCount = to_string(PatCount);
                string strTPCount = to_string(TPCount);
                string strFPCount = to_string(FPCount);
                string strFNCount = to_string(FNCount);
                string strPrecision = to_string(precision);
                string strRecall = to_string(recall);
                string strF1Score = to_string(F1Score);
    
                ofstream resultsOutfile(pr_dir + "/results_e" + strE + "_t" + strT + "_m" + strM + ".txt");
    
                if (resultsOutfile.is_open())
                {
                    resultsOutfile << "FLI-L1s:   " + to_string(FLIL1sFound) + "\n";
                    resultsOutfile << "L1Count:   " + to_string(L1Count) + "\n";
                    resultsOutfile << "PatCount:  " + strPatCount + "\n";
                    resultsOutfile << "TPCount:   " + strTPCount + "\n";
                    resultsOutfile << "Precision: " + strPrecision + "\n";
                    resultsOutfile << "Recall:    " + strRecall + "\n";
                    resultsOutfile << "F1 Score:  " + strF1Score + "\n";
                }
    
                // Save data for best m score for e and t value
                ofstream bestmOutfile(pr_dir + "/results_e" + strE + "_t" + strT + ".csv");
    
                if (bestmOutfile.is_open())
                    bestmOutfile << strM + "\t" + strPatCount + "\t" + strTPCount + "\t" + strFPCount + "\t" + strFNCount + "\t" + strPrecision + "\t" + strRecall + "\t" + strF1Score + "\n";
    
                // If current run is part of real runs or only 2 probes are availible (no need to run test run)
                if (MinKmers == MaxKmers)
                {
                    // Store all F1 scores with their respectively used parameters
                    ofstream bestmOutfile(pr_dir + "/" + k + "mers_F1Scores.csv");
    
                    if (bestmOutfile.is_open())
                        bestmOutfile << k + "\t" + strE + "\t" + strT + "\t" + strM + "\t" + strPrecision + "\t" + strRecall + "\t" + strF1Score + "\n";
    
                    if (F1Score > maxF1)
                    {
                        // Stores only F1 scores that show improvement in a given kmer value
                        ofstream bestF1Outfile(pr_dir + "/" + k + "mers_best_F1_history.csv");
    
                        if (bestF1Outfile.is_open())
                        {
                            bestF1Outfile << k + "\t" + strE + "\t" + strT + "\t" + strM + "\t" + strPrecision + "\t" + strRecall + "\t" + strF1Score + "\n";
                        }
    
                        // Stores current best F1 score to be used in if statmement comparison
                        maxF1 = F1Score;
                    }
    
                    // As best m value has been selected, skip process of selecting best m value below
                    break;
                }
                // Check if current m value gave higher F1 score than previous m values in current edit distance during test runs
                if (F1Score > BestMF1)
                {
                    // Stores best m value for a given edit distance
                    BestM = m;
                    // Stores F1 Score for BestM to be compared in if statement
                    BestMF1 = F1Score;
                }
            }
    
            else
            {
                ofstream emptyPatCountOutfile(pr_dir + "/results_e" + strE + "_t" + strT + "_m" + strM + ".csv");
    
                if (emptyPatCountOutfile.is_open())
                    emptyPatCountOutfile << "PatCount = 0";
            }
        }
            }
        }






    }
    return 0;
}