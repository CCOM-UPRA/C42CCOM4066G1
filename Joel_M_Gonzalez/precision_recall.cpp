#include <iostream>
#include <string>
#include <omp.h>
#include <string>
#include <fstream>
using namespace std;

string DIR = "/ccsopen/home/joelgonzalez35/C42CCOM4066G1";
string SAMfile = "orf12_sam/ORF12_75mers.sam";
string FASTAfile = "/ccsopen/home/joelgonzalez35/C42CCOM4066G1/orf12_probes/orf12_75mers.fasta";
string pr_dir = "precision_recall_results_75mers";
string DATA_DIR = "/ccsopen/home/joelgonzalez35/C42CCOM4066G1/hum_csv/";
int L1Count = 13671;
string k = "75";
int maxF1 = 0;
int BestM = 0;
int BestMF1 = 0;

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

int main()
{
    int MinKmers = 11;
    int MaxKmers = 12;
    int e = 15;
    int t = 100;

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

    return 0;
}
