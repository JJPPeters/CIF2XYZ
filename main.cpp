#include <iostream>
#include <istream>
#include <limits>
#include <algorithm>
#include <chrono>
#include <iomanip>

#include <getopt.h>

#include "reader.h"
#include "supercell.h"

#include <Eigen/Dense>

static int verbose_flag;

bool endsWith(std::string str, std::string suffix)
{
    // https://stackoverflow.com/questions/25777655/a-c-function-that-tests-if-the-c-string-ends-with-a-suffix
    if (str.length() < suffix.length())
        return false;
    return str.substr(str.length() - suffix.length()) == suffix;
}

std::string removeEnd(std::string str, int num)
{
    return str.substr(0, str.length() - num);
}

void makeXYZ(CIF::Reader cif, std::string output_file, double u, double v, double w, double a, double b, double c, double x_width, double y_width, double z_width, double alpha, double beta, double gamma);

/// Prints the help message
void printHelp()
{
    std::cout << "usage: cif2xyz cif_file [options]\n"
                 "  options:\n"
                 "    -h : (--help) print this help message and exit\n"
                 "    -v : (--version) print the clTEM command line version number and exit\n"
                 "    -z : (--zone) REQUIRED comma separated list for the crystal zone axis\n"
                 "    -d : (--dimensions) REQUIRED the size of the supercell in Ångströms\n"
                 "    -o : (--output) the output file name\n"
                 "    -n : (--normal) comma separated list for the axis to be aligned with the x direction\n"
                 "    -t : (--tilt) comma separated list for tilts about the x,y,z axes in degrees\n"
                 "    -f : (--fix) attempt to resolve common mistakes in the .cif file\n"
                 "    --verbose : show full output" << std::endl;
}

void printVersion() { std::cout << "cif2xyz version: 0.2a" << std::endl; }

/// Parse a string with 3 numbers separated by commas
template <typename T>
void parseThreeCommaList(std::string lst, T& a, T& b, T& c)
{
    std::istringstream ss(lst);
    std::string part;
    T v;
    std::vector<T> vec;
    while(ss.good()) {
        getline( ss, part, ',' );
//        vec.emplace_back(std::stoi(part));
        std::istringstream ssp(part);
        ssp >> v;
        vec.emplace_back(v);
    }

    if (vec.size() != 3)
        throw std::runtime_error("Error: Input vector was " + std::to_string(vec.size()) + ", expected 3.");

    a = vec[0];
    b = vec[1];
    c = vec[2];
}

void makeXYZ(CIF::Reader cif, std::string output_file, CIF::SuperCellInfo si)
{
    std::vector<std::string> Atoms;
    std::vector<double> pos_x, pos_y, pos_z, occ, ux, uy, uz;

    CIF::makeSuperCell(cif, si, Atoms, pos_x, pos_y, pos_z, occ, ux, uy, uz);

    std::cout << "Opening file: " << output_file << " for writing" << std::endl;

    std::ofstream myfile;
    myfile.open(output_file);

    myfile << std::fixed << std::showpoint << std::setprecision(5);
    myfile << Atoms.size() << "\n" << "occ nm\n";

    for(int i = 0; i < Atoms.size(); ++i)
            myfile << Atoms[i] << " " << pos_x[i] / 10 << " " << pos_y[i] / 10 << " " << pos_z[i] / 10 << " " << occ[i] << "\n";

    myfile.close();
    std::cout << "All done!" << std::endl;
}



int main(int argc, char *argv[])
{
    verbose_flag = 0; // set the here to be 100% that it is zero
    bool fix = false;
    int cc;

    std::string output_xyz = "";
    std::string input_cif = "";
    double u, v, w, a, b, c, x, y, z, alpha, beta, gamma;
    u = v = w = a = b = c = x = y = z = alpha = beta  = gamma = 0;

    bool valid = true;

    while (true)
    {
        static struct option long_options[] =
                {
                        {"help",       no_argument,       nullptr,       'h'},
                        {"version",    no_argument,       nullptr,       'v'},
                        {"zone",       required_argument, nullptr,       'z'},
                        {"dimensions", required_argument, nullptr,       'd'},
                        {"output",     optional_argument, nullptr,       'o'},
                        {"normal",     optional_argument, nullptr,       'n'},
                        {"tilt",       optional_argument, nullptr,       't'},
                        {"fix",       no_argument, nullptr,       'f'},
                        {"verbose",    no_argument,       &verbose_flag, 1},
                        {nullptr,      0,                 nullptr,       0}
                };
        // getopt_long stores the option index here.
        int option_index = 0;
        cc = getopt_long (argc, argv, "fhvz:d:o:n:t:", long_options, &option_index);

        // Detect the end of the options.
        if (cc == -1)
            break;

        switch (cc)
        {
            case 0:
                break;
            case 'h':
                printHelp();
                return 0;
            case 'v':
                printVersion();
                return 0;
            case 'o':
                output_xyz = optarg;
                break;
            case 'z':
                try {
                    parseThreeCommaList<double>(optarg, u, v, w);
                } catch (const std::exception& ex) {
                    std::cerr << "Error: Cannot parse zone axis:" << "\n\t" << ex.what() << std::endl;
                    valid = false;
                }
                break;
            case 'n':
                try {
                    parseThreeCommaList<double>(optarg, a, b, c);
                } catch (const std::exception& ex) {
                    std::cerr << "Error: Cannot parse normal vector:" << "\n\t" << ex.what() << std::endl;
                    valid = false;
                }
                break;
            case 't':
                try {
                    parseThreeCommaList<double>(optarg, alpha, beta, gamma);
                } catch (const std::exception& ex) {
                    std::cerr << "Error: Cannot parse tilt angles:" << "\n\t" << ex.what() << std::endl;
                    valid = false;
                }
                break;
            case 'd':
                try {
                    parseThreeCommaList<double>(optarg, x, y, z);
                } catch (const std::exception& ex) {
                    std::cerr << "Error: Cannot parse superstructure dimensions:" << "\n\t" << ex.what() << std::endl;
                    valid = false;
                }
                break;
            case 'f':
                fix = true;
                break;
            case '?':
                // getopt_long already printed an error message.
                break;

            default:
                return 1;
        }
    }

    // do the input file (it's the only non option arg)

    // get the non option args
    std::vector<std::string> non_option_args;
    while (optind < argc)
        non_option_args.emplace_back(argv[optind++]);

    if (non_option_args.size() < 1) {
        std::cerr << "Error: Need one non-option argument (the input .cif)" << std::endl;
        valid = false;
    }


    if (non_option_args.size() > 1) {
        std::cerr << "Error: Only expecting one non-option argument (the input .cif). Instead got:" << std::endl;
        for (std::string& s : non_option_args)
            std::cerr << "\t" << s << std::endl;
        valid = false;
    }

    // this is where the input is actually set
    if (non_option_args.size() == 1)
        input_cif = non_option_args[0];

    if (u == 0 && v == 0 && w == 0)  {
        std::cerr << "Error: Did not get valid zone axis" << std::endl;
        valid = false;
    }

    if (x <= 0 || y <= 0 || z <= 0)  {
        std::cerr << "Error: Did not get valid superstructure dimensions" << std::endl;
        valid = false;
    }

    std::ifstream test_stream(input_cif);

    std::stringstream ss;
    ss << test_stream.rdbuf();
    test_stream.close();
    std::string file = ss.str();

    if (file.empty()) {
        std::cerr << "Error: Invalid input .cif file" << std::endl;
        valid = false;
    }

    if (!valid)
        return 1;
    // we now have everything we need in some form. Need to verify the file/directory now.s

    // output file:
    // if we were not give an output file, use the input file (but change the extension)
    if (output_xyz.empty())
    {
        output_xyz = input_cif;
        if (endsWith(output_xyz, ".cif"))
            output_xyz = removeEnd(output_xyz, 4) + ".xyz";
    }

    // else just go with it, but append xyz if it doesnt have it already.
    if (!endsWith(output_xyz, ".xyz"))
        output_xyz += ".xyz";


    if (verbose_flag) {
        std::cout << "Verbose flag is set" << std::endl << std::endl;
    }

    std::shared_ptr<CIF::Reader> cif;

    try {
        cif = std::make_shared<CIF::Reader>(input_cif, fix);
    } catch (const std::exception& ex) {
        std::cerr << "Error: could not read cif file:\n\t" << ex.what() << std::endl;
        return 1;
    }

    CIF::SuperCellInfo si;
    si.setUVW(u, v, w);
    si.setABC(a, b, c);
    si.setWidths(x, y, z);
    si.setTilts(alpha, beta, gamma);

    try {
        makeXYZ(*cif, output_xyz, si);
    } catch (const std::exception& ex) {
        std::cerr << "Error: could not create superstructure:\n\t" << ex.what() << std::endl;
        return 1;
    }

    return 0;
}
