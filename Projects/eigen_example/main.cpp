#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <iostream>

using T = double;
constexpr int dim = 3;


int main(int argc, char* argv[])
{
    Eigen::Matrix<T,dim,1> x;
    x << (T).1, (T).2, (T).3;
    std::cout << x.transpose() << std::endl;
    
    return 0;
}
