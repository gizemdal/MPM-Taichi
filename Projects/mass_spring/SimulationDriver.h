#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

#include <sys/stat.h>
#include <iostream>
#include "MassSpringSystem.h"
#include <functional>

template<class T, int dim>
class SimulationDriver{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using SpMat = Eigen::SparseMatrix<T>;
    using Vec = Eigen::Matrix<T,Eigen::Dynamic,1>;

    MassSpringSystem<T,dim> ms;
    T dt;
    T frame_dt;
    TV gravity;

    TV sphere_center;
    T sphere_radius;
    std::string test;

    T collision_stiffness = 0;
    std::vector<bool> has_collision;

    std::function<void(T, T)> helper = [&](T, T) {};

    SimulationDriver()
    : dt((T)0.00001) // 0.0015 for implicit
    {
        gravity.setZero();
        gravity(1) = -9.8;

        sphere_center = TV::Zero();
        sphere_radius = 0.0;

        frame_dt = (T)1/24;
    }

    void run(const int max_frame)
    {
        // TODO: if you are not sure whether your computation is correct. You can run the following two functions to check f and df/dx
        //          when you run the simulation, comment out these functions.
        // ms.checkGradient();
        // ms.checkHessian();
        T accumulate_t = 0;
        mkdir("output/", 0777);
        std::string output_folder = "output/" + test;
        mkdir(output_folder.c_str(), 0777);
        ms.dumpPoly(output_folder + "/" + std::to_string(0) + ".poly");
        for(int frame=1; frame<=max_frame; frame++) {
            std::cout << "=============================== Frame " << frame << " ===============================" << std::endl;
            T remain_dt = frame_dt;
            T current_dt = std::min(dt, frame_dt);
            while (remain_dt > 0) {
                if (remain_dt > dt * 2)
                    current_dt = dt;
                else if (remain_dt > dt)
                    current_dt = remain_dt / 2;
                else
                    current_dt = remain_dt;
                helper(accumulate_t, current_dt);
                advanceOneStepImplicitIntegration(current_dt);
                accumulate_t += current_dt;
                remain_dt -= current_dt;
                std::cout << "Frame " << frame << " finished " << int(100 - 100 * remain_dt/frame_dt) << "%" << std::endl;
            }
            mkdir("output/", 0777);
            std::string output_folder = "output/" + test;
            mkdir(output_folder.c_str(), 0777);
            std::string filename = output_folder + "/" + std::to_string(frame) + ".poly";
            ms.dumpPoly(filename);
            std::cout << std::endl;
        }
    }

    void advanceOneStepImplicitIntegration(T dt)
    {
        int N_points = ms.x.size();
        int N_dof = dim*N_points;

        std::vector<TV> dx(N_points, TV::Zero());
        has_collision = std::vector<bool>(ms.x.size(), false);
        //initial guess and check collision
        int n_collid = 0;
        for(int p=0; p<N_points; p++){
            if (ms.node_is_fixed[p])
                dx[p] = ms.target_x[p] - ms.x[p];
            else
                dx[p] = dt * ms.v[p] + dt * dt * gravity;
            if ((ms.x[p] + dx[p] - sphere_center).norm() < 1.1 * sphere_radius){
                has_collision[p] = true;
                n_collid += 1;
            }
        }
        std::cout << n_collid << " collisions" << std::endl;
        
        int iter = 1;
        T total_length = 0;
        for (size_t i=0; i<ms.segments.size(); ++i) {
            total_length += (ms.x[ms.segments[i][0]] - ms.x[ms.segments[i][1]]).norm();
        }
        T newton_tol = 0.1 * total_length / T(ms.segments.size());
        while (true) {
            SpMat A;
            Vec gradient;
            computeHessian(dx, dt, A);
            computeGradient(dx, dt, gradient);
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<T>> solver;
            solver.compute(A);
            std::vector<TV> ddx(N_points, TV::Zero());
            Eigen::Map<Vec>(ddx[0].data(), N_dof, 1) = solver.solve(gradient);
            T ddx_norm = Eigen::Map<Vec>(ddx[0].data(),N_dof, 1).cwiseAbs().maxCoeff() / dt;
            std::cout << "It.: "<< iter <<  ", ||ddx||_inf / dt: " << ddx_norm << ", target tolerance: " << newton_tol << "                      \r"<<std::flush;
            if (ddx_norm < newton_tol) break;
            T alpha = 1;
            int n_search = 0;
            T E0 = computeEnergy(dx, dt);
            while (true) {
                std::vector<TV> new_dx(N_points, TV::Zero());
                for(int p=0; p<N_points; p++)
                    new_dx[p] = dx[p] - alpha * ddx[p];
                T E = computeEnergy(new_dx, dt);
                if (E < E0) {
                    dx = new_dx;
                    break;
                }
                else {
                    alpha *= 0.5;
                    n_search += 1;
                }
                if (n_search > 100)
                    break;
            }
            if (n_search > 100)
                break;
            iter += 1;
        }

        std::cout << std::endl; // << "# Newton iterations: " << iter << std::endl;
        for(int p=0; p<N_points; p++){
            ms.x[p] += dx[p];
            ms.v[p] = dx[p]/dt;
        }
    }

    T computeEnergy(std::vector<TV>& dx, T dt)
    {
        T total_E = 0;
        int N_points = ms.x.size();
        // inertia term
        for (int i = 0; i < N_points; ++i) {
            total_E += 0.5 * ms.m[i] * (dx[i] - ms.v[i] * dt - gravity * dt * dt).squaredNorm();
        }
        // energy on springs
        for (size_t e=0; e<ms.segments.size(); e++) {
            total_E += dt * dt * ms.evaluateSpringEnergy(e, dx);
            total_E += dt * dt * ms.evaluateDampingEnergy(e, dx, dt);
        }
        // collision
        for (int p = 0; p < N_points; ++p) {
            if (has_collision[p]) {
                T d = (ms.x[p] + dx[p] - sphere_center).norm();
                if (d < sphere_radius)
                    total_E += dt * dt * collision_stiffness * 0.5 * std::pow(1 - d / sphere_radius, 2);
            }
        }
        return total_E; 
    }

    void computeGradient(std::vector<TV>& dx, T dt, Vec& gradient)
    {
        int N_points = ms.x.size();
        int N_dof = dim*N_points;
        gradient.resize(N_dof);
        gradient.setZero();
        std::vector<TV> f_spring;
        ms.evaluateSpringForces(f_spring, dx);
	    std::vector<TV> f_damping;
	    ms.evaluateDampingForces(f_damping, dx, dt);
        for(int p=0; p<N_points; p++) {
            if (!ms.node_is_fixed[p])
                gradient.template segment<dim>(p * dim) = ms.m[p] * (dx[p] - ms.v[p] * dt - gravity * dt * dt) - dt * dt * (f_spring[p] + f_damping[p]);
        }
        // collision
        for (int p = 0; p < N_points; ++p) {
            if (has_collision[p] && !ms.node_is_fixed[p]) {
                T d = (ms.x[p] + dx[p] - sphere_center).norm();
                if (d < sphere_radius) {
                    TV n = (ms.x[p] + dx[p] - sphere_center).normalized();
                    gradient.template segment<dim>(p * dim) += - dt * dt * collision_stiffness * (1 - d/sphere_radius) * n / sphere_radius;
                }
            }
        }
    }

    void computeHessian(std::vector<TV>& dx, T dt, SpMat& A, bool project_spd=true) 
    {
        int N_points = ms.x.size();
        int N_dof = dim*N_points;
        A.resize(N_dof,N_dof);
        A.setZero();
        A.reserve(Eigen::VectorXi::Constant(N_dof,dim*40)); // estimate non-zero entries per column

        // Mass matrix contribution (assembly to global)
        for(int p=0; p<N_points; p++) {
            for (int d = 0; d < dim; d++) {
                int i = p * dim + d; // global dof index
                A.coeffRef(i, i) = ms.m[p];
            }
        }

        for(size_t e=0; e<ms.segments.size(); e++)
        {
            int particle[2]; // global particle index
            particle[0] = ms.segments[e](0);
            particle[1] = ms.segments[e](1);

            // Stiffness matrix contribution
            Eigen::Matrix<T,dim,dim> Ks = ms.evaluateKS(e, dx, project_spd);
            Eigen::Matrix<T,dim*2,dim*2> K_local;
            K_local.template block<dim,dim>(0,0) = Ks;
            K_local.template block<dim,dim>(dim,0) = -Ks;
            K_local.template block<dim,dim>(0,dim) = -Ks;
            K_local.template block<dim,dim>(dim,dim) = Ks;

            // Damping matrix contribution
            Eigen::Matrix<T,dim,dim> Kd = ms.evaluateKD(e, dt);
            Eigen::Matrix<T,dim*2,dim*2> G_local;
            G_local.template block<dim,dim>(0,0) = Kd;
            G_local.template block<dim,dim>(dim,0) = -Kd;
            G_local.template block<dim,dim>(0,dim) = -Kd;
            G_local.template block<dim,dim>(dim,dim) = Kd;

            // Stiffness and damping contribution assembly to global
            for(int p=0; p<2; p++){
                for(int q=0; q<2; q++){
                    if(!ms.node_is_fixed[particle[p]] && !ms.node_is_fixed[particle[q]]) {
                        for (int i = 0; i < dim; i++) {
                            for (int j = 0; j < dim; j++) {
                                A.coeffRef(dim * particle[p] + i, dim * particle[q] + j) -= dt * dt * (K_local(dim * p + i,
                                                                                                    dim * q + j) + G_local(dim*p+i, dim*q+j));
                            }
                        }
                    }
                }
            }
        }

        // collision
        for (int p = 0; p < N_points; ++p) {
            if (has_collision[p] && !ms.node_is_fixed[p]) {
                T d = (ms.x[p] + dx[p] - sphere_center).norm();
                if (d < sphere_radius) {
                    TV n = (ms.x[p] + dx[p] - sphere_center).normalized();
                    Eigen::Matrix<T,dim,dim> C_local = collision_stiffness / sphere_radius / sphere_radius * n * n.transpose() 
                            - collision_stiffness * 2 * (1 - d/sphere_radius) / d / sphere_radius * (Eigen::Matrix<T,dim,dim>::Identity() - n * n.transpose());
                    if (project_spd)
                        makePD(C_local);
                    for (int i = 0; i < dim; i++) {
                        for (int j = 0; j < dim; j++) {
                            A.coeffRef(dim * p + i, dim * p + j) += dt * dt * C_local(i, j);
                        }
                    }
                }
            }
        }

        // process dirichlet-0 nodes
        for (size_t p=0; p<ms.node_is_fixed.size(); p++){
            if(ms.node_is_fixed[p]){
                for(int d=0; d<dim; d++) {
                    A.coeffRef(dim * p + d, dim * p + d) = 1;
                }
            }
        }
        A.makeCompressed();
    }
};
