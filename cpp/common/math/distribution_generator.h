#include "spline_hermite.h"
#include "fastmath.h"
#include <fstream>
#include <iostream>
#include <string>
#include <quaternion.h>
#include <sstream>


class DistributionGenerator{
    public:
    int distributionType = -1;   // 3 - Student's t-distribution with 1 Dof
    double* Eps = 0;
    int n;                  // number of data lines is in data file
    double dx;
    double x0;
    double inv_dx;
    
    DistributionGenerator(){}
    DistributionGenerator( int distributionType_ ){
        init( distributionType_ );
    }
    ~DistributionGenerator(){
        if(Eps) delete Eps;
    }


    void init( int distributionType_ ){
        distributionType = distributionType_;
    }

    inline double generateRandom(){
        double y;
        switch( distributionType ){
            case 0:
                y = randf();
                break;
            case 1:
                y = gauss();
                break;
            case 3:
                y = student(randf());
                break;
        }
        return y;
    }
    inline double student( double x_ ){
        if(!Eps){
            loadEps("data/functions/invStudent_1.dat");
        }
        return splineDistribution( x_ );
    }
    // Box-Muller transform
    inline double gauss(){ 
        return sqrt(-2*log(randf()))*cos(2*M_PI*randf());
    }

    void loadEps( const char* fname ){
        std::ifstream file(fname);
        if (file) {
            std::string line;
            if (std::getline(file, line)) {
                std::istringstream iss(line);
                iss >> n >> dx >> x0;
                Eps = new double[n];
                inv_dx = 1/dx;
            }
            int i = 0;
            while (std::getline(file, line)) {
                std::istringstream iss(line);
                iss >> Eps[i];
                i++;
            }
            file.close();
        }
        else
        {
            std::cout << "Unable to open file";
            exit(1);
        }
    }

    inline double splineDistribution( double x_ ){
        double x;
        double u = (x_-x0)*inv_dx; 
        int i = (int)u;
        double du = u-i;
        if(i<0){ i=0; du=0; }
        if(i>n-3){ i=n-3; du=1; }
        Quat4d c = *(Quat4d*)(Eps+i);
        x = Spline_Hermite::val(x_, c.y, c.z, (c.z-c.x)*0.5, (c.w-c.y)*0.5);
        return x;
    }
    
    void randomPointOnNdimSphere( std::vector<double>& p, int N){
        double sum = 0;
        for(int i=0; i<N; i++){
            p.push_back(randf()-0.5);
            sum += p[i]*p[i];
        }
        sum = sqrt(sum);
        for(int i=0; i<p.size(); i++){
            p[i] /= sum;
        }
    }
};