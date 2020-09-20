#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <chrono>
#include <unordered_set>

#include "SimulationDriver.h"

int main(int argc, char* argv[])
{
    using T = float;
    constexpr int dim = 3;
    using TV = Eigen::Matrix<T,dim,1>;

    SimulationDriver<T,dim> driver;

    // set up mass spring system
    T youngs_modulus = 0;
    T damping_coeff = 0; 
    T dt = 0;

    // node data
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;
    std::vector<bool> node_is_fixed;

    // segment data
    std::vector<Eigen::Matrix<int,2,1> > segments;
    std::vector<T> rest_length;

    if (argc < 2) 
    {
        std::cout << "Please indicate test case number: 0 (cloth) or 1 (volumetric bunny)" << std::endl;
        exit(0);
    }

    if (strcmp(argv[1], "0") == 0) // cloth case
    {
        // TODO
        /* 
            1. Create node data: position, mass, velocity
            2. Fill segments and rest_length, including struct springs, shearing springs and bending springs.
            3. Choose proper youngs_modulus, damping_coeff and dt.
            4. Set boundary condition (node_is_fixed) and helper function (to achieve moving boundary condition).
            5. Generate quad mesh for rendering.
        */
        driver.helper = [&](T t, T dt) {
            // TODO
        };
        driver.test="cloth";
    }

    else if (strcmp(argv[1], "1") == 0) // volumetric bunny case
    { 
        // TODO
        /* 
            1. Create node data from data/points: The first line indicates the number of points and dimension (which is 3). 
            2. Fill segments and rest_length from data/cells: The first line indicates the number of tetrahedra and the number of vertices of each tet (which is 6). Each edge in this tetrahedral mesh will be a segment. Be careful not to create duplicate edges. 
            3. Choose proper youngs_modulus, damping_coeff, dt; 
            4. Set boundary condition (node_is_fixed) and helper function (to achieve moving boundary condition).
        */

        driver.helper = [&](T t, T dt) {
            // TODO
        };
        driver.test="bunny";
    }

    else {
        std::cout << "Wrong case number!" << std::endl;
        exit(0);
    }

    // simulate
    
    driver.dt = dt;
    driver.ms.segments = segments;
    driver.ms.m = m;
    driver.ms.v = v;
    driver.ms.x = x;
    driver.ms.youngs_modulus = youngs_modulus;
    driver.ms.damping_coeff = damping_coeff;
    driver.ms.node_is_fixed = node_is_fixed;
    driver.ms.rest_length = rest_length;

    driver.run(120);

    return 0;
}
