#include <Partio.h>
#include <vector>
#include <string>

using T = double;
constexpr int dim = 3;

template <class T, int dim>
void writePartio(const std::string& particleFile)
{
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, vH, mH;
    mH = parts->addAttribute("m", Partio::VECTOR, 1);
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    vH = parts->addAttribute("v", Partio::VECTOR, 3);
    for (int i=0; i<3; i++){
        int idx = parts->addParticle();
        float* m = parts->dataWrite<float>(mH, idx);
        float* p = parts->dataWrite<float>(posH, idx);
        float* v = parts->dataWrite<float>(vH, idx);
        m[0] = (T)i;
        for (int k = 0; k < 3; k++)
            p[k] = (T)i;
        for (int k = 0; k < 3; k++)
            v[k] = (T)i;
    }

    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}


int main(int argc, char* argv[])
{
    std::string file="test.bgeo";
    writePartio<T,dim>(file);
    
    return 0;
}
