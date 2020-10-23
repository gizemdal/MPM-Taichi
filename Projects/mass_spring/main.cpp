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
#include <iomanip>

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
    T dt = T(1/24.);

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

    if (strcmp(argv[1], "1") == 0) { // bunny case
        youngs_modulus = 0.1; // TODO: iterate in [0.1, 1, 10, 100, 1000]
        damping_coeff = 2;
        // TODO: 
        /* 
            1. Copy the loading codes from your hw1. Fix two ears (2140, 2346) only, and you don't need helper function here.
            2. Set the initial velocities of non_fixed_nodes to be (10, 0, 0)
            
            The output folder will automatically renamed by bunny_[youngs_modulus], don't worry about overwriting.
        */
        std::stringstream ss;
        ss << std::fixed << std::setprecision(2) << youngs_modulus;
        driver.test="bunny_"+ss.str();
    }

    else if (strcmp(argv[1], "2") == 0) { //brush case
        driver.gravity.setZero();
        driver.collision_stiffness = 0.1;
        youngs_modulus = 10000;
        damping_coeff = 100; // 0
        int N = 32; // z direction
        int M = 4; // y direction
        int L = 32; // x direction
        int N_points = N*M*L;
        T dx = (T)0.1/(N-1);
        m.resize(N_points);
        x.resize(N_points);
        v.resize(N_points);
        node_is_fixed = std::vector<bool>(N_points, false);
        for(int i=0; i<N; i++){ // z
            for(int j=0; j<M; j++) { // y
                for (int k=0; k<L; k++) { // x
                    int id = i * M * L + j * L + k;
                    m[id] = (T)0.001/N_points;
                    x[id](2) = i*dx;
                    x[id](1) = j*dx;
                    x[id](0) = k*dx;
                    v[id] = TV::Zero();
                    if (k <= 2) node_is_fixed[id] = true;
                    // struct spring
                    if (k > 0) { 
                        segments.push_back(Eigen::Matrix<int,2,1>(id, id-1));
                        rest_length.push_back((x[id]-x[id-1]).norm());
                    }
                    // bending spring
                    if (k > 1) {
                        segments.push_back(Eigen::Matrix<int,2,1>(id, id-2));
                        rest_length.push_back((x[id]-x[id-2]).norm());
                    }
                }
            }
        }

        driver.sphere_radius = 0.04;
        driver.sphere_center = TV(0.07, -0.045, 0.05);

        driver.helper = [&](T t, T dt) {
            if(t < 4) {
                // driver.sphere_center = TV(0.12, t/4. * 0.15 + (1-t/4.) * (-0.06), 0.05);
                for(size_t i = 0; i < driver.ms.x.size(); ++i) {
                    if (driver.ms.node_is_fixed[i]) {
                        driver.ms.target_x[i](1) -= 0.15 / 4 * dt;
                    }
                }
            }
        };

        driver.test="brush";
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
    driver.ms.target_x = x;
    driver.ms.youngs_modulus = youngs_modulus;
    driver.ms.damping_coeff = damping_coeff;
    driver.ms.node_is_fixed = node_is_fixed;
    driver.ms.rest_length = rest_length;

    driver.run(180);

    return 0;
}
