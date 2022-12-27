#include <highfive/H5File.hpp>
#include <vector>
#include <random>
#include <list>
#include <chrono>
#include <algorithm>
#include <map>
#include <getopt.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sys/stat.h>
#include <omp.h>

using namespace std;

void PrintHelp()
{
    cout << "--tssLen, -k:        define the mean of pause sites across cells [default: 50]\n"
            "--kSd:               define the standard deviation of pause sites across cells [default: 0]\n"
            "--kMin:              upper bound of pause site allowed [default: 17]\n"
            "--kMax:              lower bound of pause site allowed [default: 200]\n"
            "--geneLen:           define the length of the whole gene [default: 2e3]\n"
            "--alpha, -a:         initiation rate [default: 1 event per min]\n"
            "--beta, -b:          pause release rate [default: 1 event per min]\n"
            "--zeta, -z:          the mean of elongation rates across sites [default: 2000 per min]\n"
            "--zetaSd:            the standard deviation of elongation rates across sites [default: 1000]\n"
            "--zetaMax:           the maximum of elongation rates allowed [default: 2500 per min]\n"
            "--zetaMin:           the minimum of elongation rates allowed [default: 1500 per min]\n"
            "--cellNum, -n:       number of cells being simulated [default: 10]\n"
            "--polSize, -s:       Polymerase II size [default: 33]\n"
            "--addSpace:          Additional space in addition to RNAP size [default: 17]\n"
            "--time, -t:          Total time of simulating data in a cell [default: 0.1 min]\n"
            "--hdf5:              Record position matrix to HDF5 file for remaining number of steps specified [default: 0 steps]\n"
            "--csv:               Record position matrix to csv file for remaining number of steps specified [default: 1 step]"
            "--help:              Show help\n";
    exit(1);
}

template <typename T>
void PrintMatrixToCSV(vector<T> &matrix, int nrows, int ncols, string file_name)
{
    ofstream out(file_name);

    for (int i = 0; i < nrows; i++)
    {
        vector<int> *sites = &matrix[i];
        int site_idx = 0;
        for (int j = 0; j < ncols; j++)
        {
            if (j == (*sites)[site_idx])
            {
                out << "1,";
                site_idx++;
            }
            else
            {
                out << "0,";
            }
        }
        out << '\n';
    }
}

template <typename T>
void PrintVectorToCSV(const vector<T> &input, string file_name, string type)
{
    ofstream out(file_name);

    for (size_t i = 0; i < input.size(); i++)
    {
        out << type << " " << i << ": " << input[i] << ',' << '\n';
    }
}

void ConvertSiteDataToMatrix(vector<vector<int>> &input, vector<vector<int>> &output)
{
    for (size_t i = 0; i < output.size(); i++)
    {
        vector<int> *sites = &input[i];
        for (size_t j = 0; j < sites->size(); j++)
        {
            output[i][(*sites)[j]] = 1;
        }
    }
}

class Generator
{
    default_random_engine generator;
    normal_distribution<double> distribution;
    double min;
    double max;

public:
    Generator(double mean, double stddev, double min, double max) : distribution(mean, stddev), min(min), max(max)
    {
    }

    double operator()()
    {
        while (true)
        {
            double number = this->distribution(generator);
            if (number >= this->min && number <= this->max)
                return number;
        }
    }
};

int main(int argc, char **argv)
{
    /* set defaults for parameters */
    int k = 50;                    // mean pause site across cells
    double ksd = 0;                // standard deviation of pause sites across cells
    static int k_min = 17;         // upper bound of pause sites allowed
    static int k_max = 200;        // lower bound of pause sites allowed
    size_t gene_len = 2e3;         // length of the entire gene
    double alpha = 1;              // initiation rate
    double beta = 1;               // pause release rate
    double zeta = 2000;            // mean of elongation rates across sites
    double zeta_sd = 1000;         // standard deviation of elongation rates across sites
    static double zeta_max = 2500; // Max elongation rates allowed
    static double zeta_min = 1500; // Min elongation rates allowed
    int total_cells = 10;
    int s = 33;        // polymerase II size
    int h = 17;        // Additional space in addition to RNAP size
    double time = 0.1; // Total time of simulating data in a cell in minutes
    double delta_t = 1e-4;
    int hdf5_steps_to_record = 0;
    int csv_steps_to_record = 1;

    const char *const short_opts = "k:a:b:z:n:s:t:h";
    const option long_opts[] = {
        {"tssLen", required_argument, 0, 'k'},
        {"kSd", required_argument, 0, 0},
        {"kMin", required_argument, 0, 0},
        {"kMax", required_argument, 0, 0},
        {"geneLen", required_argument, 0, 0},
        {"alpha", required_argument, 0, 'a'},
        {"beta", required_argument, 0, 'b'},
        {"zeta", required_argument, 0, 'z'},
        {"zetaSd", required_argument, 0, 0},
        {"zetaMax", required_argument, 0, 0},
        {"zetaMin", required_argument, 0, 0},
        {"cellNum", required_argument, 0, 'n'},
        {"polSize", required_argument, 0, 's'},
        {"addSpace", required_argument, 0, 0},
        {"time", required_argument, 0, 't'},
        {"hdf5", required_argument, 0, 0},
        {"csv", required_argument, 0, 0},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, no_argument, nullptr, 0}};
    int option_index = 0;

    while (true)
    {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, &option_index);

        if (-1 == opt)
            break;

        switch (opt)
        {
        case 0:
            if (long_opts[option_index].flag != 0)
                break;
            if (strcmp(long_opts[option_index].name, "kSd") == 0)
                ksd = stod(optarg);
            if (strcmp(long_opts[option_index].name, "kMin") == 0)
                k_min = stoi(optarg);
            if (strcmp(long_opts[option_index].name, "kMax") == 0)
                k_max = stoi(optarg);
            if (strcmp(long_opts[option_index].name, "geneLen") == 0)
                gene_len = stoi(optarg);
            if (strcmp(long_opts[option_index].name, "alpha") == 0)
                alpha = stod(optarg);
            if (strcmp(long_opts[option_index].name, "beta") == 0)
                beta = stod(optarg);
            if (strcmp(long_opts[option_index].name, "zeta") == 0)
                zeta = stod(optarg);
            if (strcmp(long_opts[option_index].name, "zetaSd") == 0)
                zeta_sd = stod(optarg);
            if (strcmp(long_opts[option_index].name, "zetaMax") == 0)
                zeta_max = stod(optarg);
            if (strcmp(long_opts[option_index].name, "zetaMin") == 0)
                zeta_min = stod(optarg);
            if (strcmp(long_opts[option_index].name, "cellNum") == 0)
                total_cells = stoi(optarg);
            if (strcmp(long_opts[option_index].name, "polSize") == 0)
                s = stoi(optarg);
            if (strcmp(long_opts[option_index].name, "addSpace") == 0)
                h = stoi(optarg);
            if (strcmp(long_opts[option_index].name, "time") == 0)
                time = stod(optarg);
            if (strcmp(long_opts[option_index].name, "hdf5") == 0)
                hdf5_steps_to_record = stoi(optarg);
            if (strcmp(long_opts[option_index].name, "csv") == 0)
                csv_steps_to_record = stoi(optarg);
            break;
        case 'k':
            k = stoi(optarg);
            break;
        case 'a':
            alpha = stod(optarg);
            break;
        case 'b':
            beta = stod(optarg);
            break;
        case 'z':
            zeta = stod(optarg);
            break;
        case 'n':
            total_cells = stoi(optarg);
            break;
        case 's':
            s = stoi(optarg);
            break;
        case 't':
            time = stod(optarg);
            break;
        case 'h':
        case '?':
        default:
            PrintHelp();
            break;
        }
    }

    int steric_hindrance = s + h;
    const int total_sites = gene_len + 1;
    double steps = time / delta_t;

    /* Create output directory */
    mkdir("results", 0755);
    if(csv_steps_to_record > 0)
    {
        mkdir("results/positions", 0755);
    }

    /* Initialize an array to hold Pol II presence and absence*/
    vector<vector<int>> pos_matrix;
    for (int i = 0; i < total_cells; i++)
    {
        vector<int> sites;
        sites.push_back(0);
        pos_matrix.insert(pos_matrix.begin() + i, sites);
    }

    /* Construct a probability matrix to control RNAP movement
     * Generate pause sites located from kmin to kmax with sd = ksd
     */
    vector<double> y;
    Generator y_distribution(k, ksd, k_min, k_max);
    for (int i = 0; i < total_cells; i++)
    {
        y.push_back(round(y_distribution()));
    }
    /* Output pause sites in csv format */
    PrintVectorToCSV(y, "results/pause_sites.csv", "cell");

    /* A matrix of probabilities to control transition from state to state
     * cols are cells, rows are positions
     */
    vector<double> zv;
    Generator zv_distribution(zeta, zeta_sd, zeta_min, zeta_max);
    for (int i = 0; i < total_sites; i++)
    {
        zv.push_back(zv_distribution() * delta_t);
    }
    /* Output probability values per site in csv format */
    PrintVectorToCSV(zv, "results/probability_vector.csv", "site");

    random_device rd; // Get seed for random number generator
    mt19937 gen(rd());
    uniform_real_distribution<double> distrib(0.0, 1.0);

    auto start = chrono::high_resolution_clock::now();
    vector<vector<vector<int>>> pos_matrices_csv_record; 
    if(hdf5_steps_to_record > 0)
    {
        HighFive::File file("results/position_matrices.h5", HighFive::File::Create);
    }
    for (int step = 0; step < steps; step++)
    {
#pragma omp parallel for
        for (int cell = 0; cell < total_cells; cell++)
        {
            vector<int> *sites = &pos_matrix[cell];
            for (size_t i = 0; i < sites->size(); i++)
            {
                /* Determine whether polymerase can move or not
                 * criteria 1, probability larger than random draw
                 * criteria 2, enough space ahead to let polymerase advance
                 */
                int idx = (*sites)[i];
                double prob = idx == 0            ? alpha * delta_t
                              : idx == y.at(cell) ? beta * delta_t
                                                  : zv.at(idx);
                double draw = distrib(gen);
                if (prob > draw)
                {
                    size_t last_polymerase = sites->size() - 1;
                    /* Check if space ahead is larger than polymerase size */
                    if (i != last_polymerase && (*sites)[i + 1] - (*sites)[i] > steric_hindrance)
                    {
                        (*sites)[i]++;
                    }
                    /* Always allow the polymerase at the end to move */
                    else if (i == last_polymerase)
                    {
                        if ((*sites)[i] + 1 < total_sites)
                        {
                            (*sites)[i]++;
                        }
                        else
                        {
                            /* Remove polymerase if past final site */
                            sites->pop_back();
                            break; // to prevent iterating past the end of the linked list
                        }
                    }
                }
            }

            /* Ensure there are always polymerases waiting to be initialized (i.e., first row equals 1) */
            if (sites->size() == 0 || (*sites)[0] != 0)
            {
                sites->insert(sites->begin(), 0);
            }
        }
        /* Record info for studying steric hindrance */
        bool record_to_hdf5 = step >= (steps - hdf5_steps_to_record);
        bool record_to_csv = step >= (steps - csv_steps_to_record);
        if (record_to_hdf5)
        {
            HighFive::File file("results/position_matrices.h5", HighFive::File::ReadWrite);
            vector<size_t> dims(total_cells, total_sites);
            vector<vector<int>> pos_matrix_to_record(total_cells, vector<int>(total_sites));
            ConvertSiteDataToMatrix(pos_matrix, pos_matrix_to_record);
            HighFive::DataSet dataset = file.createDataSet<int>("/group/dataset_" + to_string(step), HighFive::DataSpace(dims));
            dataset.write(pos_matrix_to_record);
        }
        if (record_to_csv)
        {
            pos_matrices_csv_record.push_back(pos_matrix);
        }
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    printf("Total time used for the simulation is %.2f mins.\n", duration.count() / 60.0);

    /* Calculate the # of polymerase at each site across all cells and output the results in csv format*/
    vector<int> res_all(total_sites, 0);
    for (int i = 0; i < total_cells; i++)
    {
        vector<int> *sites = &pos_matrix[i];
        for (size_t j = 0; j < sites->size(); j++)
        {
            res_all[(*sites)[j]]++;
        }
    }
    PrintVectorToCSV(res_all, "results/combined_cell_data.csv", "site");

    if(csv_steps_to_record > 0)
    {
        #pragma omp parallel for
        for(int i = 0; i < csv_steps_to_record; i++)
        {
            /* Output position matrix in csv format */
            PrintMatrixToCSV(pos_matrices_csv_record[i], total_cells, total_sites, "results/positions/position_matrix_" + to_string((int)steps - csv_steps_to_record + i) + ".csv");       
        }
    }

    return 0;
}