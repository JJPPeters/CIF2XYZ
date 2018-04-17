#include <iostream>
#include <limits>
#include <algorithm>
#include <chrono>
#include <iomanip>

#include "CIFReader.h"

#include <Eigen/Dense>

std::string double2string(double dbl, int precision = 2)
{
    std::ostringstream strs;
    strs << std::setprecision(precision) << dbl;
    std::string str = strs.str();
    return str;
}

void makeXYZ(CIFReader cif, std::string output_file, double u, double v, double w, double a, double b, double c, double x_width, double y_width, double z_width, double alpha, double beta, double gamma);

int main()
{
    // get file to open
//    std::string folder = "D:\\Dropbox\\Documents\\PhD physics\\Analysis\\Geanina PbTiO\\Structures\\LSMO structure\\";
//    std::string file = "PTO_ferro";
    std::string folder = "D:\\Dropbox\\Documents\\PhD physics\\Data\\CIF\\";
    std::string file = "znblende";


//    std::string out_folder = folder + file;
//    CreateDirectory(out_folder.c_str(), NULL);

    double u = 0, v = 0, w = 1; // zone axis of crystal (aka will be z direction)
    double a = 1, b = 0, c = 0;
//    double u = 1, v = -1, w = -1; // zone axis of crystal (aka will be z direction)
//    double a = 1, b = 1, c = 0;
    double x_width = 50, y_width = 50, z_width = 200;

//    double alpha, beta, gamma;
//    gamma = 0.0;

    auto cif = CIFReader(folder + file + ".cif");

    std::string out_file = folder + file + ".xyz";
    makeXYZ(cif, out_file, u, v, w, a, b, c, x_width, y_width, z_width, 0.0, 0.0, 0.0);

//    for (alpha = -2.0; alpha < 2.1; alpha += 0.2)
//        for (beta = -2.0; beta < 2.1; beta += 0.2)
//        {
//            std::string out = folder + file + "\\" + file + "_alpha_" + double2string(alpha) + "_beta_" + double2string(beta) + ".xyz";
//            makeXYZ(cif, out, u, v, w, a, b, c, x_width, y_width, z_width, alpha, beta, gamma);
//        }

    return 0;
}

void makeXYZ(CIFReader cif, std::string output_file, double u, double v, double w, double a, double b, double c, double x_width, double y_width, double z_width, double alpha, double beta, double gamma)
{
//    std::chrono::time_point<std::chrono::system_clock> start, end;
//    start = std::chrono::system_clock::now();
//    // USER SET VARIABLES
//    // these will eventually be (command line?) arguments to function so they cna be user set
//    double u = 1, v = 0, w = 0; // zone axis of crystal (aka will be z direction)
//    double a = 0, b = 1, c = 0; // crystal direction of x in image (projected)
//    double alpha = 0, beta = 0, gamma = 0; // tilt of crystal in degrees
//    double x_width = 50, y_width = 50, z_width = 200;


    // here just make call to function to create xyz file
//    auto cif = CIFReader("D:\\GoogleDrive\\PhDPhysics\\Reports Presentations\\Domain walls as new 2D functional materials poster\\Images\\STO sim\\STO.cif");

    UnitCell  cell = cif.getUnitCell();

    // small test to print out basic unit cell information
//    for (auto at : cell.getAtoms())
//    {
//        std::cout << "Element: ";
//        for ( auto n : at.getElements())
//            std::cout << n << " ";
//        std::cout << std::endl;
//
//        std::cout << "Occupancy: ";
//        for ( auto o : at.getOccupancies())
//            std::cout << o << " ";
//        std::cout << std::endl;
//
//        std::cout << "Fractional positions: " << std::endl;
//        for (auto p : at.getPositions())
//            std::cout << p(0) << ", " << p(1) << ", " << p(2) << std::endl;
//    }

    // TODO: rename this method (get cell axes?)
    auto geom  = cell.getCellGeometry();
    auto a_vector = Eigen::Vector3d(geom.getAVector().data());
    //std::cout << "A vector: " << a_vector(0) << ", " << a_vector(1) << ", " << a_vector(2) << std::endl;
    auto b_vector = Eigen::Vector3d(geom.getBVector().data());
    //std::cout << "B vector: " << b_vector(0) << ", " << b_vector(1) << ", " << b_vector(2) << std::endl;
    auto c_vector = Eigen::Vector3d(geom.getCVector().data());
    //std::cout << "C vector: " << c_vector(0) << ", " << c_vector(1) << ", " << c_vector(2) << std::endl;

    // this gives us the element type, occupancy and FRACTIONAL coords as well as the basis vectors
    // now we just need to tile them to fill a space?

    // create vector from the zone axis we want
    Eigen::Vector3d uvw = u * a_vector + v * b_vector + w * c_vector;
    // create matrix with direction we want to map onto (here is is the z direction)
    Eigen::Vector3d z_direction(std::vector<double>({0.0, 0.0, 1.0}).data());
    // create the rotation matrix
    auto za_rotation = Utilities::generateNormalisedRotationMatrix<double>(uvw, z_direction);

    // similarly for rotation in x,y plane
    Eigen::Vector3d abc = a * a_vector + b * b_vector + c * c_vector;
    abc = za_rotation * abc;
    abc(2) = 0.0; // no rotate the z-axis!
    abc.normalize();
    Eigen::Vector3d x_direction(std::vector<double>({1.0, 0.0, 0.0}).data());
    // Need the negative angle here? so give inputs in opposite order?
    auto xy_rotation = Utilities::generateNormalisedRotationMatrix<double>(x_direction, abc);

//    double theta = -abc.dot(x_direction) / (abc.norm() * x_direction.norm());
//
//    double t_cos = std::cos(theta);
//    double t_sin = std::sin(theta);
//
//    std::vector<double> data = {t_cos, -t_sin, 0.0, t_sin, t_cos, 0.0, 0.0, 0.0, 1.0};
//    Eigen::Matrix3d xy_rotation(data.data());
//
//    std::vector<double> test1(za_rotation.size());
//    Eigen::Matrix3d::Map(&test1[0], test1.size()) = za_rotation;
//
//    std::vector<double> test2(xy_rotation.size());
//    Eigen::Matrix3d::Map(&test2[0], test2.size()) = xy_rotation;

    // convert angles to radians
    alpha = alpha * PI / 180;
    beta = beta * PI / 180;
    gamma = gamma * PI / 180;

    // generate the small rotation matrices
    auto x_rotation = Utilities::generateRotationMatrix<double>({1.0, 0.0, 0.0}, alpha);
    auto y_rotation = Utilities::generateRotationMatrix<double>({0.0, 1.0, 0.0}, beta);
    auto z_rotation = Utilities::generateRotationMatrix<double>({0.0, 0.0, 1.0}, gamma);

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

    std::ofstream myfile;
    myfile.open (output_file);

    myfile << std::fixed << std::showpoint << std::setprecision(5);
    myfile << valid_count << "\n" << "occ\n";

    for(int i = 0; i < count; ++i)
        if (position_valid[i])
            myfile << element_list[i] << " " << (x_list[i] - x_low)/10 << " " << (y_list[i] - y_low)/10 << " " << (z_list[i] - z_low)/10 << " " << occupancy_list[i] << "\n";

    myfile.close();

//    end = std::chrono::system_clock::now();
//    std::chrono::duration<double> elapsed_seconds = end - start;

//    std::cout << "Created .xyz in " << elapsed_seconds.count() << "s" << std::endl;

//    return 0;
}