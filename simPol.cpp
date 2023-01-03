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
            "--zetaVec            a file contains vector to scale elongation rates. All cells share the same set of parameters [default: ""]"
            "--cellNum, -n:       number of cells being simulated [default: 10]\n"
            "--polSize, -s:       Polymerase II size [default: 33]\n"
            "--addSpace:          Additional space in addition to RNAP size [default: 17]\n"
            "--time, -t:          Total time of simulating data in a cell [default: 0.1 min]\n"
            "--hdf5:              Record position matrix to HDF5 file for remaining number of steps specified [default: 0 steps]\n"
            "--csv:               Record position matrix to csv file for remaining number of steps specified [default: 1 step]\n"
            "--outputDir, -d:     Directory for saving results [default: 'results']"
            "--help:              Show help\n";
    exit(1);
}

template <typename T>
void PrintPositionMatrixToCSV(vector<T> &matrix, int nrows, int ncols, string file_name)
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

vector<double> NormalDistrubtionGenerator(double mean, double stddev, double min, double max, size_t length, bool round_result, double multiplication_factor=1)
{
    random_device rd; // Get seed for random number generator
    default_random_engine generator;
    normal_distribution<double> distribution(mean, stddev);
    generator.seed(rd());
    vector<double> random_values;
    while(random_values.size() < length)
    {
        double number = distribution(generator);
        if (number >= min && number <= max) {
            if(round_result)
            {
                number = round(number);
            }
            random_values.push_back(number * multiplication_factor);
        }
    }

    return random_values;
}

int main(int argc, char **argv)
{
    /* set defaults for parameters */
    int k = 50, k_min = 17, k_max = 200;     // mean, max and min of pause sites
    double ksd = 0;                          // std dev of pause sites across cells
    size_t gene_len = 2e3;                   // length of the entire gene
    double alpha = 1, beta = 1;              // initiation rate, pause release rate
    double zeta = 2000, zeta_sd = 1000, zeta_max = 2500, zeta_min = 1500; // mean, std dev, max, and min of elongation rates across sites
    string zeta_vec = "";                  // file containing vector to scale elongation rates
    int total_cells = 10;
    int s = 33, h = 17;                      // polymerase II size, Additional space in addition to RNAP size
    double time = 0.1, delta_t = 1e-4;       // total time of simulating data in a cell in minutes
    int hdf5_steps_to_record = 0, csv_steps_to_record = 1;
    string output_dir = "results";

    const char *const short_opts = "k:a:b:z:n:s:t:d:h";
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
        {"zetaVec", required_argument, 0, 0},
        {"cellNum", required_argument, 0, 'n'},
        {"polSize", required_argument, 0, 's'},
        {"addSpace", required_argument, 0, 0},
        {"time", required_argument, 0, 't'},
        {"hdf5", required_argument, 0, 0},
        {"csv", required_argument, 0, 0},
        {"outputDir", required_argument, 0, 'd'},
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
            if (strcmp(long_opts[option_index].name, "zetaVec") == 0)
                zeta_vec = optarg;
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
            if (strcmp(long_opts[option_index].name, "outputDir") == 0)
                output_dir = optarg;
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
        case 'd':
            output_dir = optarg;
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

    /* Create output directories */
    string positions_dir = output_dir + "/positions";
    string pause_sites_file_name = output_dir + "/pause_sites.csv";
    string probability_file_name = output_dir + "/probability_vector.csv";
    string positions_hdf5_file_name = output_dir + "/position_matrices.h5";
    string combined_cells_file_name = output_dir + "/combined_cell_data.csv";
    string positions_file_name = positions_dir + "/position_matrix_";
    mkdir(output_dir.c_str(), 0755);
    if(csv_steps_to_record > 0)
    {
        mkdir(positions_dir.c_str(), 0755);
    }
    ofstream out;

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
    vector<double> y = NormalDistrubtionGenerator(k, ksd, k_min, k_max, total_cells, true);

    /* Output pause sites per cell in csv format */
    out.open(pause_sites_file_name);
    for (size_t i = 0; i < y.size(); i++)
    {
        out << y[i] << '\n';
    }
    out.close();

    /* A matrix of probabilities to control transition from state to state
     * cols are cells, rows are positions
     */
    vector<double> zv;
    if(zeta_vec == "")
    {
        zv = NormalDistrubtionGenerator(zeta, zeta_sd, zeta_min, zeta_max, total_sites, false, delta_t);
    }
    else {
        ifstream  data(zeta_vec);
        if(data.is_open())
        {
            string line;
            while(getline(data, line, '\n')) {
                zv.emplace_back(stod(line));
            }
        }
        if((int)zv.size() >= total_sites)
        {
            zv.resize(total_sites);
        }
        else if((int)zv.size() == total_sites - 1)
        {
            double mean = accumulate(zv.begin(), zv.end(), 0.0) / zv.size();
            zv.insert(zv.begin(), mean);
        }
        else {
            printf("Vector for scaling zeta is too short, check total length of the vector!");
            return -1;
        }
        double transform_val = zeta * delta_t;
        transform(zv.begin(), zv.end(), zv.begin(), [&transform_val](auto& c){return c*transform_val;});
    }
    /* Output probability values per site in csv format */
    out.open(probability_file_name);
    for (size_t i = 0; i < zv.size(); i++)
    {
        out << zv[i] << '\n';
    }
    out.close();

    random_device rd; // Get seed for random number generator
    mt19937 gen(rd());
    uniform_real_distribution<double> distrib(0.0, 1.0);

    auto start = chrono::high_resolution_clock::now();
    vector<vector<vector<int>>> pos_matrices_csv_record; 
    if(hdf5_steps_to_record > 0)
    {
        HighFive::File file(positions_hdf5_file_name, HighFive::File::Create | HighFive::File::Overwrite);
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
                int site_idx = (*sites)[i];
                double prob = site_idx == 0            ? alpha * delta_t
                              : site_idx == y.at(cell) ? beta * delta_t
                                                  : zv.at(site_idx);
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

            /* Ensure there are always polymerases waiting to be initialized (i.e., first col is always 1) */
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
            HighFive::File file(positions_hdf5_file_name, HighFive::File::ReadWrite);
            vector<size_t> dims(2);
            dims[0] = total_cells;
            dims[1] = total_sites;
            vector<vector<int>> pos_matrix_record_to_hdf5(total_cells, vector<int>(total_sites));
            ConvertSiteDataToMatrix(pos_matrix, pos_matrix_record_to_hdf5);
            HighFive::DataSetCreateProps dsprops;
            hsize_t num_chunk_rows = total_cells / 10;
            hsize_t num_chunk_cols = total_sites / 10;
            dsprops.add(HighFive::Chunking(vector<hsize_t>{num_chunk_rows >= 1 ? num_chunk_rows : 1, num_chunk_cols >= 1 ? num_chunk_cols : 1}));
            dsprops.add(HighFive::Deflate(9));
            HighFive::DataSet dataset = file.createDataSet<int>("/group/dataset_" + to_string(step + 1), HighFive::DataSpace(dims), dsprops);
            dataset.write(pos_matrix_record_to_hdf5);
        }
        if (record_to_csv)
        {
            pos_matrices_csv_record.push_back(pos_matrix);
        }
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    printf("Total time used for the simulation is %.2f mins.\n", duration.count() / 60.0);

    /* Calculate the # of polymerase at each site across all cells and output the results in csv format */
    vector<int> res_all(total_sites, 0);
    for (int cell = 0; cell < total_cells; cell++)
    {
        vector<int> *sites = &pos_matrix[cell];
        for (size_t j = 0; j < sites->size(); j++)
        {
            res_all[(*sites)[j]]++;
        }
    }
    out.open(combined_cells_file_name);
    for (size_t i = 0; i < res_all.size(); i++)
    {
        out << res_all[i] << "\n";
    }
    out.close();

    if(csv_steps_to_record > 0)
    {
        int total_steps_to_record = csv_steps_to_record > (int)pos_matrices_csv_record.size() ? (int)pos_matrices_csv_record.size() : csv_steps_to_record;
        #pragma omp parallel for
        for(int i = 0; i < total_steps_to_record; i++)
        {
            int step_idx = csv_steps_to_record > steps ? i + 1 : steps - csv_steps_to_record + i + 1; 
            PrintPositionMatrixToCSV(pos_matrices_csv_record[i], total_cells, total_sites, positions_file_name + to_string(step_idx) + ".csv");       
        }
    }

    return 0;
}