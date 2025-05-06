#include <iostream>
#include <fstream>
#include <filesystem>
#include <cstdlib>
#include <string>
#include <omp.h>

int main() {
    std::string printdir = "output_dir"; // directorio de salida
    std::string pr_dir = printdir;
    std::string k = "21"; // k-mer size
    std::string genome = "ATCG..."; // secuencia
    std::string FASTAfile = "reads.fasta";
    std::string GENOME = "reference_genome";

    int EditStart = 1, EditEnd = 5, EditStep = 1;
    int OriginalMaxKmers = 100, OriginalMinKmers = 10;

    std::filesystem::create_directories(pr_dir);

    std::ofstream bestF1File(pr_dir + "/" + k + "mers_best_F1_history.csv");
    std::ofstream allF1File(pr_dir + "/" + k + "mers_F1Scores.csv");

    std::string header = "K-mer size\tEdit dist\tThreshold\tK-mers min\tPrecision\tRecall\tF1 Score\n";
    bestF1File << header;
    allF1File << header;

    bestF1File.close();
    allF1File.close();

    std::ofstream testFile("test.txt", std::ios::app);
    testFile << genome << "\n";
    testFile.close();

    // Paralelizamos el bucle con OpenMP
    #pragma omp parallel for
    for (int e = EditStart; e <= EditEnd; e += EditStep) {
        // para evitar conflictos genera su propio archivo de salida 
        std::string samFileName = pr_dir + "/output_e" + std::to_string(e) + ".sam";

        // Proteger el chequeo del índice con una sección crítica
        #pragma omp critical
        {
            if (!std::filesystem::exists(GENOME + ".index")) {
                std::string indexCmd = "mrfast --index " + GENOME;
                std::system(indexCmd.c_str());
            }
        }

        // Ejecutar mrfast para este valor de e
        std::string searchCmd = "mrfast --search " + GENOME +
                                " --seq " + FASTAfile +
                                " -e " + std::to_string(e) +
                                " -o " + samFileName;

        std::cout << "[Hilo " << omp_get_thread_num() << "] Ejecutando: " << searchCmd << std::endl;
        std::system(searchCmd.c_str());

        // Reiniciar variables
        int MaxKmers = OriginalMaxKmers;
        int MinKmers = OriginalMinKmers;
        int BestM = 0;
        double BestMF1 = 0.0;

    }

    return 0;
}
