#include <highfive/H5File.hpp>
#include <vector>
#include <random>
#include <list>
#include <chrono>
#include <map>
#include <getopt.h>
#include <fstream>
#include <sys/stat.h>
#include <omp.h>

using namespace std;

void PrintHelp()
{
    cout << "--verbose, -v:       print messages [default]\n"
            "--quietly, -q:       print no messages\n"
            "--tssLen, -k:        define the mean of pause sites across cells [default: 50]\n"
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
            "--help:              Show help\n";
    exit(1);
}

template <typename T>
void PrintMatrixToCSV(const vector<T> &matrix, string file_name)
{
    ofstream out(file_name);

    for (auto &row : matrix)
    {
        for (auto col : row)
            out << col << ',';
        out << '\n';
    }
}

void ConvertListDataToMatrix(vector<list<int>> &input, vector<vector<int>> &output)
{
    for (int cell = 0; cell < output.size(); cell++)
    {
        list<int> *sites = &input[cell];
        for (list<int>::iterator itr = sites->begin(); itr != sites->end(); itr++)
        {
            output[cell][*itr] = 1;
        }
    }
}

int main(int argc, char **argv)
{
    /* set defaults for parameters */
    bool verbose = true;
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
    int total_cells = 1000;
    int s = 33;       // polymerase II size
    int h = 17;       // Additional space in addition to RNAP size
    double time = 10; // Total time of simulating data in a cell in minutes
    int length = gene_len - k;
    int steric_hindrance = s + h;
    double delta_t = 1e-4;
    double steps = time / delta_t;
    const int total_sites = 2e3 + 1;
    int steps_to_record = 100;

    int c;

    const char *const short_opts = "k:a:b:z:n:s:t:d:";
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
        {"addSpace", required_argument, 0, 's'},
        {"time", required_argument, 0, 't'},
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

    /* Create output directory */
    mkdir("results", 0755);

    /* Initialize an array to hold Pol II presence and absence*/
    vector<list<int>> pos_matrix;
    for (int i = 0; i < total_cells; i++)
    {
        list<int> sites;
        sites.push_front(0);
        pos_matrix.insert(pos_matrix.begin() + i, sites);
    }

    /* Construct a probability matrix to control RNAP movement
     * Generate pause sites located from kmin to kmax with sd = ksd
     */
    vector<double> y;
    static default_random_engine generator;
    while (y.size() < total_cells)
    {
        static normal_distribution<double> distribution(k, ksd); // mean, std dev
        vector<double> x(total_cells * 2);
        generate(x.begin(), x.end(), []()
                 { return distribution(generator); });
        x.erase(remove_if(x.begin(), x.end(), [](const double &v)
                          { return v < k_min + 1 || v > k_max + 1; }),
                x.end());
        y.insert(y.end(), x.begin(), x.end());
    }
    for_each(y.begin(), y.end(), [](double &v)
             { round(v); });

    /* A matrix of probabilities to control transition from state to state
     * cols are cells, rows are positions
     * ignore gamma for now
     */
    vector<double> zv;
    while (zv.size() < total_sites * total_cells)
    {
        static normal_distribution<double> distribution(zeta, zeta_sd); // mean, std dev
        vector<double> x(total_sites * total_cells * 2);
        generate(x.begin(), x.end(), []()
                 { return distribution(generator); });
        x.erase(remove_if(x.begin(), x.end(), [](const double &v)
                          { return v < zeta_min || v > zeta_max; }),
                x.end());
        zv.insert(zv.end(), x.begin(), x.end());
    }
    zv.resize(total_sites * total_cells);
    transform(zv.begin(), zv.end(), zv.begin(), [&delta_t](auto &c)
              { return c * delta_t; });

    vector<vector<double>> prob_matrix(total_cells, vector<double>(total_sites, 0));
    int zv_idx = 0;
    for (int i = 0; i < total_cells; i++)
    {
        for (int j = 0; j < total_sites; j++)
        {
            prob_matrix[i][j] = j == 0         ? alpha * delta_t
                                : j == y.at(i) ? beta * delta_t
                                               : zv.at(zv_idx);
            zv_idx++;
        }
    }

    int total_length = total_sites * total_cells;

    random_device rd; // Get seed for random number generator
    mt19937 gen(rd());
    uniform_real_distribution<double> distrib(0.0, 1.0);

    auto start = chrono::high_resolution_clock::now();
    HighFive::File file("results/position_matrices.h5", HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    vector<size_t> dims(2);
    dims[0] = total_cells;
    dims[1] = total_sites;
    for (int step = 0; step < steps; step++)
    {
        #pragma omp parallel for
        for (int cell = 0; cell < total_cells; cell++)
        {
            list<int> *sites = &pos_matrix[cell];
            list<int>::iterator itr = sites->begin();
            while (itr != sites->end())
            {
                /* Determine whether polymerase can move or not
                 * criteria 1, probability larger than random draw
                 * criteria 2, enough space ahead to let polymerase advance
                 */
                double prob = prob_matrix[cell][*itr];
                double draw = distrib(gen);
                if (prob > draw)
                {
                    /* Check if space ahead is larger than polymerase size */
                    if (next(itr) != sites->end() && *next(itr) - *itr > steric_hindrance)
                    {
                        *itr = *itr + 1;
                    }
                    /* Always allow the polymerase at the end to move */
                    else if (next(itr) == sites->end())
                    {
                        int new_index = *itr + 1;
                        if (new_index < total_sites)
                        {
                            *itr = new_index;
                        }
                        else
                        {
                            /* Remove polymerase if past final site */
                            sites->pop_back();
                            break; // to prevent iterating past the end of the linked list
                        }
                    }
                }
                itr++;
            }

            /* Ensure there are always polymerases waiting to be initialized (i.e., first row equals 1) */
            if (sites->size() == 0 || sites->front() != 0)
            {
                sites->push_front(0);
            }
        }
        /* Record info for studying steric hindrance */
        bool record_to_hdf5 = step >= (steps - steps_to_record);
        bool final_step = step == steps - 1;
        if (record_to_hdf5 || final_step)
        {
            vector<vector<int>> pos_matrix_to_record(total_cells, vector<int>(total_sites));
            ConvertListDataToMatrix(pos_matrix, pos_matrix_to_record);
            if(record_to_hdf5)
            {
                HighFive::DataSet dataset = file.createDataSet<int>("/group/dataset" + to_string(step), HighFive::DataSpace(dims));
                dataset.write(pos_matrix_to_record);
            }
            if(final_step)
            {
                /* Output final position matrix in csv format */
                PrintMatrixToCSV(pos_matrix_to_record, "results/final_position_matrix.csv");
            }
        }
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    printf("Total time used for the simulation is %.2f mins.\n", duration.count() / 60.0);

    /* Calculate the # of polymerase at each site across all cells*/
    map<int, int> res_all;
    for (int j = 0; j < total_sites; j++)
    {
        res_all[j] = 0;
    }
    for (int i = 0; i < total_cells; i++)
    {
        list<int> *sites = &pos_matrix[i];
        list<int>::iterator itr = sites->begin();
        for (itr = sites->begin(); itr != sites->end(); ++itr)
        {
            res_all[*itr]++;
        }
    }

    /* Output initial probability matrix in csv format */
    PrintMatrixToCSV(prob_matrix, "results/probability_matrix.csv");

    return 0;
}