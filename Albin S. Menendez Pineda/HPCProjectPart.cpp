#include <iostream>
#include <string>
#include <omp.h>
// Possibly unnecessary
#include <fstream>
using namespace std;
// string PR_DIR = "";
// string RESULTS_FILE_NAME = "";
// string DATA_DIR = "";

// Functions from Joel's part of the program are used here.  REDUNDANT AF
void runPythonScript(const string script_path, const string args, bool time)
{
    string command = "python3 " + script_path + " " + args;
    if (time)
        command = "time " + command;

    int result = system(command.c_str());
    if (result != 0)
        cerr << "Error executing: " << command << endl;
}

int runCommand(const string command, string &output)
{
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
}

// vvvvv Actual program part. vvvvv
int main()
{
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

            // I believe this is where Joel's code would sit in... SEND HEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEELP oh nevermind, I just wanted to be seen :)

            // Checking if the run done was a test run
            if(MinKmers!=MaxKmers)
            {
                // Resetting threshold value back to prepare for real run
                t-=ThresholdStep;

                // Prepare future runs of this edit distance to use -m=BestM instead of the range of m values
                MinKmers=BestM
                MaxKmers=BestM
            }

            // Close the file to avoid unwanted editing
            outFile.close();
        }
        // Just in case of error.  Feel free to remove or change (veeeeery unsure if it works or not... never used "cerr" before)
        else
        {
            cerr<<"Unable to open file for writing."<<endl;
        }
    }
}