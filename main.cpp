#include <iostream>
#include <istream>
#include <limits>
#include <algorithm>
#include <chrono>
#include <iomanip>

#include <getopt.h>

#include "CIFReader.h"

#include <Eigen/Dense>

static int verbose_flag;



std::string double2string(double dbl, int precision = 2)
{
    std::ostringstream strs;
    strs << std::setprecision(precision) << dbl;
    std::string str = strs.str();
    return str;
}

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

void makeXYZ(CIFReader cif, std::string output_file, double u, double v, double w, double a, double b, double c, double x_width, double y_width, double z_width, double alpha, double beta, double gamma);

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
                 "    --verbose : show full output" << std::endl;
}

void printVersion() { std::cout << "cif2xyz version: 0.1a" << std::endl; }

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

void makeXYZ(CIFReader cif, std::string output_file, double u, double v, double w, double a, double b, double c, double x_width, double y_width, double z_width, double alpha, double beta, double gamma)
{
    UnitCell  cell = cif.getUnitCell();

    // TODO: rename this method (get cell axes?)
    auto geom  = cell.getCellGeometry();
    auto a_vector = Eigen::Vector3d(geom.getAVector().data());
    auto b_vector = Eigen::Vector3d(geom.getBVector().data());
    auto c_vector = Eigen::Vector3d(geom.getCVector().data());

    if (verbose_flag)
        std::cout << "Unit cell with bases:\na:\n" << a_vector << "\nb:\n" << b_vector << "\nc:\n" << c_vector << std::endl;

    // this gives us the element type, occupancy and FRACTIONAL coords as well as the basis vectors
    // now we just need to tile them to fill a space?

    // create vector from the zone axis we want
    Eigen::Vector3d uvw = u * a_vector + v * b_vector + w * c_vector;
    // create matrix with direction we want to map onto (here is is the z direction)
    Eigen::Vector3d z_direction(std::vector<double>({0.0, 0.0, 1.0}).data());
    // create the rotation matrix
    auto za_rotation = Utilities::generateNormalisedRotationMatrix<double>(uvw, z_direction);

    if (za_rotation.array().isNaN().sum() != 0) {
        throw std::runtime_error("Zone axis rotation matrix contains NaNs");
    }
    if (verbose_flag)
        std::cout << "Zone axis rotation:\n" << za_rotation << std::endl;

    // similarly for rotation in x,y plane, though here we default to 0 rotation
    Eigen::Matrix3d xy_rotation = Eigen::Matrix3d::Identity();
    if (a != 0 || b != 0 || c != 0)
    {
        Eigen::Vector3d abc = a * a_vector + b * b_vector + c * c_vector;
        abc = za_rotation * abc;
        abc(2) = 0.0; // no rotate the z-axis!
        abc.normalize();
        Eigen::Vector3d x_direction(std::vector<double>({1.0, 0.0, 0.0}).data());
        // Need the negative angle here? so give inputs in opposite order?
        xy_rotation = Utilities::generateNormalisedRotationMatrix<double>(x_direction, abc);
    }

    if (xy_rotation.array().isNaN().sum() != 0) {
        throw std::runtime_error("Normal rotation matrix contains NaNs");
    }
    if (verbose_flag)
        std::cout << "Normal rotation:\n" << xy_rotation << std::endl;

    // convert angles to radians
    alpha = alpha * PI / 180;
    beta = beta * PI / 180;
    gamma = gamma * PI / 180;

    // generate the small rotation matrices
    auto x_rotation = Utilities::generateRotationMatrix<double>({1.0, 0.0, 0.0}, alpha);
    auto y_rotation = Utilities::generateRotationMatrix<double>({0.0, 1.0, 0.0}, beta);
    auto z_rotation = Utilities::generateRotationMatrix<double>({0.0, 0.0, 1.0}, gamma);

    if (verbose_flag)
        std::cout << "Tilt adjustments:\nalpha:\n" << x_rotation << "\nbeta:\n" << y_rotation << "\ngamma:\n" << z_rotation << std::endl;

    // for easy loopage
    std::vector<Eigen::Vector3d> basis = {a_vector, b_vector, c_vector};

    // apply all rotations in order
    for (int i = 0; i < basis.size(); ++i)
    {
        basis[i] = za_rotation * basis[i];
        basis[i] = xy_rotation * basis[i];
        basis[i] = x_rotation * basis[i];
        basis[i] = y_rotation * basis[i];
        basis[i] = z_rotation * basis[i];
    }

    if (verbose_flag)
        std::cout << "Rotated bases:\na:\n" << basis[0] << "\nb:\n" << basis[1] << "\nc:\n" << basis[2] << std::endl;

    // AFFINE TRANSFORM STUFF

    // make the basis matrices
    Eigen::MatrixXd affine_basis(3, 4);
    affine_basis << 0.0, 1.0, 0.0, 0.0,
                    0.0, 0.0, 1.0, 0.0,
                    0.0, 0.0, 0.0, 1.0;
    Eigen::Matrix4d real_basis;
    real_basis << 1.0, 1.0, 1.0, 1.0,
                  0.0, basis[0](0), basis[1](0), basis[2](0),
                  0.0, basis[0](1), basis[1](1), basis[2](1),
                  0.0, basis[0](2), basis[1](2), basis[2](2);

    auto affine_transform_matrix = affine_basis * real_basis.inverse();

    // make matrix representing the corners of our cuboid
    Eigen::MatrixXd cuboid(8, 4);
    cuboid << 1.0, 0.0, 0.0, 0.0,
              1.0, x_width, 0.0, 0.0,
              1.0, 0.0, y_width, 0.0,
              1.0, 0.0, 0.0, z_width,
              1.0, x_width, y_width, 0,
              1.0, x_width, 0.0, z_width,
              1.0, 0.0, y_width, z_width,
              1.0, x_width, y_width, z_width;

    // affine transform the cuboid
    Eigen::MatrixXd affine_cube(8, 3);
    for (int i = 0; i < cuboid.rows(); ++i)
    {
        Eigen::Vector4d temp_row = cuboid.row(i);
//        affine_cube.setRow(Mat::Multiply(affine_transform_matrix, temp_row), i);
        affine_cube.row(i) = affine_transform_matrix * temp_row;
    }

    // get the limits of the cube in affine space
    double x_max = std::numeric_limits<double>::min();
    double x_min = std::numeric_limits<double>::max();
    double y_max = std::numeric_limits<double>::min();
    double y_min = std::numeric_limits<double>::max();
    double z_max = std::numeric_limits<double>::min();
    double z_min = std::numeric_limits<double>::max();

    for (int i = 0; i < affine_cube.rows(); ++i)
    {
        if (affine_cube(i, 0) > x_max)
            x_max = affine_cube(i, 0);
        if (affine_cube(i, 0) < x_min)
            x_min = affine_cube(i, 0);
        if (affine_cube(i, 1) > y_max)
            y_max = affine_cube(i, 1);
        if (affine_cube(i, 1) < y_min)
            y_min = affine_cube(i, 1);
        if (affine_cube(i, 2) > z_max)
            z_max = affine_cube(i, 2);
        if (affine_cube(i, 2) < z_min)
            z_min = affine_cube(i, 2);
    }

    // get integer cell range (rounded up to nearest even number)
    int x_range = static_cast<int>(std::ceil(x_max - x_min));
    if (x_range % 2 != 0)
        x_range += 1;
    int y_range = static_cast<int>(std::ceil(y_max - y_min));
    if (y_range % 2 != 0)
        y_range += 1;
    int z_range = static_cast<int>(std::ceil(z_max - z_min));
    if (z_range % 2 != 0)
        z_range += 1;

    // get total number of atoms in our xyz file
    int count = 0;
    for (auto at : cell.getAtoms())
        count += at.getElements().size() * at.getPositions().size();
    count = count * x_range * y_range * z_range;

    if (verbose_flag)
        std::cout << "Calculated tiling needed to span cell:\nx:\n" << x_range << "\ny:\n" << y_range << "\nz:\n" << z_range << std::endl;
    if (x_range == 0 || y_range == 0 || z_range == 0 )
        throw std::runtime_error("Did not find any atoms inside limits");

    std::vector<std::string> element_list(count);
    std::vector<double> x_list(count);
    std::vector<double> y_list(count);
    std::vector<double> z_list(count);
    std::vector<double> occupancy_list(count);
    int it = 0;


    for (auto at : cell.getAtoms())
    {
        for (auto pos : at.getPositions())
        {
            // convert from fractional coordinates
            auto p = pos[0] * basis[0] + pos[1] * basis[1] + pos[2] * basis[2];
            for (int k = 0; k < z_range; ++k)
            {
                //std::cout << "k: " << k << std::endl;
                auto k_factor = basis[2] * k;
                for (int j = 0; j < y_range; ++j)
                {
                    //std::cout << "j: " << j << std::endl;
                    auto j_factor = basis[1] * j;
                    for (int i = 0; i < x_range; ++i)
                    {
                        //std::cout << "i: " << i << std::endl;
                        auto i_factor = basis[0] * i;
                        auto new_pos = p + i_factor + j_factor + k_factor;

                        for (int ind = 0; ind < at.getElements().size(); ++ind)
                        {
                            element_list[it] = at.getElements()[ind];
                            x_list[it] = new_pos(0);
                            y_list[it] = new_pos(1);
                            z_list[it] = new_pos(2);
                            occupancy_list[it] = at.getOccupancies()[ind];
                            ++it;
                        }
                    }
                }
            }
        }
    }

    if (x_list.size() == 0 || y_list.size() == 0 || z_list.size() == 0 )
        throw std::runtime_error("Did not find any atoms inside limits");

    // TODO: all this can be put into a function
    // re use the min and max variables here
    x_min = *std::min_element(x_list.begin(), x_list.end());
    x_max = *std::max_element(x_list.begin(), x_list.end());
    y_min = *std::min_element(y_list.begin(), y_list.end());
    y_max = *std::max_element(y_list.begin(), y_list.end());
    z_min = *std::min_element(z_list.begin(), z_list.end());
    z_max = *std::max_element(z_list.begin(), z_list.end());
    // TODO: find mid point of new crystal
    // TODO: find crop limits

    double x_mid = (x_min + x_max) / 2;
    double y_mid = (y_min + y_max) / 2;
    double z_mid = (z_min + z_max) / 2;

    double x_low = x_mid - (x_width / 2);
    double y_low = y_mid - (y_width / 2);
    double z_low = z_mid - (z_width / 2);

    double x_high = x_mid + (x_width / 2);
    double y_high = y_mid + (y_width / 2);
    double z_high = z_mid + (z_width / 2);

    std::vector<bool> position_valid(count, false);
    int valid_count = 0;

    // find valid positions
    for (int i = 0; i < count; ++i)
        if (x_list[i] > x_low && x_list[i] < x_high)
            if (y_list[i] > y_low && y_list[i] < y_high)
                if (z_list[i] > z_low && z_list[i] < z_high)
                {
                    position_valid[i] = true;
                    valid_count++;
                }

    std::cout << "Opening file: " << output_file << " for writing" << std::endl;

    std::ofstream myfile;
    myfile.open(output_file);

    myfile << std::fixed << std::showpoint << std::setprecision(5);
    myfile << valid_count << "\n" << "occ nm\n";

    for(int i = 0; i < count; ++i)
        if (position_valid[i])
            myfile << element_list[i] << " " << (x_list[i] - x_low)/10 << " " << (y_list[i] - y_low)/10 << " " << (z_list[i] - z_low)/10 << " " << occupancy_list[i] << "\n";

    myfile.close();
    std::cout << "All done!" << std::endl;
}



int main(int argc, char *argv[])
{
    verbose_flag = 0; // set the here to be 100% that it is zero

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
                        {"verbose",    no_argument,       &verbose_flag, 1},
                        {nullptr,      0,                 nullptr,       0}
                };
        // getopt_long stores the option index here.
        int option_index = 0;
        cc = getopt_long (argc, argv, "hvz:d:o:n:t:", long_options, &option_index);

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

    std::shared_ptr<CIFReader> cif;

    try {
        cif = std::make_shared<CIFReader>(input_cif);
    } catch (const std::exception& ex) {
        std::cerr << "Error: could not read cif file:\n\t" << ex.what() << std::endl;
        return 1;
    }

    try {
        makeXYZ(*cif, output_xyz, u, v, w, a, b, c, x, y, z, alpha, beta, gamma);
    } catch (const std::exception& ex) {
        std::cerr << "Error: could not create superstructure:\n\t" << ex.what() << std::endl;
        return 1;
    }

    return 0;
}
