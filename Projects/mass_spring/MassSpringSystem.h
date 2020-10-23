#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>

template <typename Scalar, int size>
void makePD(Eigen::Matrix<Scalar, size, size>& symMtr)
{
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, size, size>> eigenSolver(symMtr);
    if (eigenSolver.eigenvalues()[0] >= 0.0) {
        return;
    }
    Eigen::DiagonalMatrix<Scalar, size> D(eigenSolver.eigenvalues());
    int rows = ((size == Eigen::Dynamic) ? symMtr.rows() : size);
    for (int i = 0; i < rows; i++) {
        if (D.diagonal()[i] < 0.0) {
            D.diagonal()[i] = 0.0;
        }
        else {
            break;
        }
    }
    symMtr = eigenSolver.eigenvectors() * D * eigenSolver.eigenvectors().transpose();
}

template<class T, int dim>
class MassSpringSystem{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using TM = Eigen::Matrix<T,dim,dim>;
    std::vector<Eigen::Matrix<int,2,1> > segments;
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> target_x;
    std::vector<TV> v;
    T youngs_modulus;
    T damping_coeff;
    std::vector<bool> node_is_fixed;
    std::vector<T> rest_length;

    MassSpringSystem()
    {}

    void evaluateSpringForces(std::vector<TV >& f, std::vector<TV>& dx)
    {
        f.clear();
        f.resize(x.size(), TV::Zero());
        for (size_t i=0; i<segments.size(); ++i)
        {
            TV force = evaluateSpringForce(i, dx);
            int A = segments[i](0);
            int B = segments[i](1);
            f[A] += force;
            f[B] -= force;
        }
    }

    void evaluateDampingForces(std::vector<TV >& f, std::vector<TV>& dx, T dt)
    {
        f.clear();
        f.resize(x.size(), TV::Zero());
        for (size_t i=0; i<segments.size(); ++i)
        {
            TV force = evaluateDampingForce(i, dx, dt);
            int A = segments[i](0);
            int B = segments[i](1);
            f[A] += force;
            f[B] -= force;
        }
    }

    void dumpPoly(std::string filename)
    {
        std::ofstream fs;
        fs.open(filename);
        fs << "POINTS\n";
        int count = 0;
        for (auto X : x) {
            fs << ++count << ":";
            for (int i = 0; i < dim; i++)
                fs << " " << X(i);
            if (dim == 2)
                fs << " 0";
            fs << "\n";
        }
        fs << "POLYS\n";
        count = 0;
        for (const Eigen::Matrix<int, 2, 1>& seg : segments)
            fs << ++count << ": " << seg(0) + 1 << " " << seg(1) + 1 << "\n"; // poly segment mesh is 1-indexed
        fs << "END\n";
        fs.close();
    }

    // for implicit integrator
    T evaluateSpringEnergy(int idx, std::vector<TV>& dx) {
        int A = segments[idx](0);
        int B = segments[idx](1);
        TV xA = x[A] + dx[A];
        TV xB = x[B] + dx[B];
        T l = (xA - xB).norm();
        return 0.5 * rest_length[idx] * youngs_modulus * std::pow(l/rest_length[idx]-(T)1, 2);
    }

    T evaluateDampingEnergy(int idx, std::vector<TV>& dx, T dt) {
        int A = segments[idx](0);
        int B = segments[idx](1);
        TV vA = dx[A]/dt;
        TV vB = dx[B]/dt;
        TV n = (x[A]-x[B]).normalized();
        return 0.5 * dt * damping_coeff * std::pow(n.dot((vA-vB)), 2);
    }

    TV evaluateSpringForce(int idx, std::vector<TV>& dx) {
        // TODO: evaluate fA_spring only, i.e., The spring force on the first node of the segment.
        //       compute xA and xB as follows:
        //       int A = segments[idx](0);
        //       int B = segments[idx](1);
        //       TV xA = x[A] + dx[A];
        //       TV xB = x[B] + dx[B];
        return TV();
    }

    TV evaluateDampingForce(int idx, std::vector<TV>& dx, T dt) {
        // TODO: evaluate fA_damping only, i.e., The damping force on the first node of the segment.
        //       compute vA, vB, and n as follows:
        //       int A = segments[idx](0);
        //       int B = segments[idx](1);
        //       TV vA = dx[A]/dt;
        //       TV vB = dx[B]/dt;
        //       TV n = (x[A]-x[B]).normalized();
        return TV();
    }

    TM evaluateKS(int idx, std::vector<TV>& dx, bool project_pd=true) {
        // TODO: evaluate dfA_spring/dxA and store it in KS, 
        //       i.e. the derivaitive of the spring force on the first node of the segment w.r.t the position of the first node.
        //       compute xA and xB as follows:
        //       int A = segments[idx](0);
        //       int B = segments[idx](1);
        //       TV xA = x[A] + dx[A];
        //       TV xB = x[B] + dx[B];
        // Do not change the project_pd part (-Ks is projected to be PD).

        TM KS; // TODO: fill this matrix

        // write your code above and do not change the following part.
        TM neg_KS = -KS;
        if (project_pd)
            makePD(neg_KS);
        return -neg_KS;
    }

    TM evaluateKD(int idx, T dt) {
        // TODO: evaluate dfA_damping/dxA and store it in KS, 
        //       i.e. the derivaitive of the spring force on the first node of the segment w.r.t the position of the first node.
        //       compute n as follows:
        //       int A = segments[idx](0);
        //       int B = segments[idx](1);
        //       TV n = (x[A]-x[B]).normalized(); 
        //       Be careful that dfA_damping/dxA = dfA_damping/dvA * dvA/dxA = 1/dt * dfA_damping/dvA
        return TM();
    }

    void checkGradient() {
        T eps = 1e-6;
        T dt = 0.01;
        std::vector<TV> dx(x.size());
        T delta = (x[0] - x[1]).norm();
        for (size_t i=0; i<x.size(); ++i) {
            dx[i].setZero();
            for (int d = 0; d < dim; ++d) {
                dx[i](d) += 0.1 * ((T) rand() / RAND_MAX + T(0.5)) * delta;
            }
        }
        for (size_t i=0; i<segments.size(); ++i) {
            int A = segments[i](0);
            TV spring_finite_diff;
            TV spring_analytical = -evaluateSpringForce(i, dx);
            TV damping_finite_diff;
            TV damping_analytical = -evaluateDampingForce(i, dx, dt);
            for(int d = 0; d < dim; ++d) {
                dx[A][d] += eps;
                T se1 = evaluateSpringEnergy(i, dx);
                T de1 = evaluateDampingEnergy(i, dx, dt);
                dx[A][d] -= 2 * eps;
                T se2 = evaluateSpringEnergy(i, dx);
                T de2 = evaluateDampingEnergy(i, dx, dt);
                dx[A][d] += eps;
                spring_finite_diff(d) = (se1 - se2) / (T(2) * eps);
                damping_finite_diff(d) = (de1 - de2) / (T(2) * eps);
            }
            std::cout << "Check gradient on segment " << i << std::endl 
                      << "Analytical spring energy gradient: " << spring_analytical.transpose() << ", finite difference: " << spring_finite_diff.transpose() << std::endl
                      << "Analytical damping energy gradient: " << damping_analytical.transpose() << ", finite difference: " << damping_finite_diff.transpose() << std::endl;
            getchar();
        }
    }

    void checkHessian() {
        T eps = 1e-6;
        T dt = 0.01;
        std::vector<TV> dx(x.size());
        T delta = (x[0] - x[1]).norm();
        for (size_t i=0; i<x.size(); ++i) {
            dx[i].setZero();
            for (int d = 0; d < dim; ++d) {
                dx[i](d) += ((T) rand() / RAND_MAX + T(0.5)) * delta;
            }
        }
        for (size_t i=0; i<segments.size(); ++i) {
            int A = segments[i](0);
            TM spring_finite_diff;
            TM spring_analytical = -evaluateKS(i, dx, false);
            TM damping_finite_diff;
            TM damping_analytical = -evaluateKD(i, dt);
            for(int d = 0; d < dim; ++d) {
                dx[A][d] += eps;
                TV sf1 = -evaluateSpringForce(i, dx);
                TV df1 = -evaluateDampingForce(i, dx, dt);
                dx[A][d] -= 2 * eps;
                TV sf2 = -evaluateSpringForce(i, dx);
                TV df2 = -evaluateDampingForce(i, dx, dt);
                dx[A][d] += eps;
                spring_finite_diff.col(d) = (sf1 - sf2) / (T(2) * eps);
                damping_finite_diff.col(d) = (df1 - df2) / (T(2) * eps);
            }
            std::cout << "Check hessian on segment " << i << std::endl
                      << "analytical spring hessian: \n" << spring_analytical << ", \n finite difference: \n" << spring_finite_diff << std::endl
                      << "analytical damping hessian: \n" << damping_analytical << ", \n finite difference: \n" << damping_finite_diff << std::endl;
            getchar();
        }
    }

};
