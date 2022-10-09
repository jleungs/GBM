#include <stdio.h>
#include <cstdlib>
#include <chrono>
#include <random>
#include <cmath>
#include <blitz/array.h>

class BM {
    /* drift coefficent */
    double Mu;
    /* volatility */
    double Sigma;
    /* time in years */
    int T;
    /* number of simulations */
    int M;
    /* number of steps */
    int N = 252;
    /* initial price */
    int S0 = 100;

    void sim_path() {
        int i, j;
        double t0, t1, dt;
        /* generate n*T x m matrix */
        blitz::Array<double, 2> paths(this->N*this->T, this->M);
        /* fill first row with S0 */
        for (i=0; i < this->M; i++)
            paths(0, i) = S0;
        /* compute time steps */
        dt = (double)this->T / (double)this->N;
        /* loop for each column*/
        for (j=0; j < paths.columns(); j++) {
            /* loop for each row */
            for (i=1; i < paths.rows(); i++) {
                /* generate seed from time */
                unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                /* generate a normal random variable N(0, sqrt(dt)) */
                std::default_random_engine generator(seed);
                std::normal_distribution<double> normal(0, std::sqrt(dt));
                /* exponential terms */
                t0 = (this->Mu - (std::pow(this->Sigma, 2)/2))*dt;
                t1 = this->Sigma * normal(generator);
                /* St =  */
                paths(i,j) = paths(i-1,j) * std::exp(t0 + t1);
                //std::cout << paths(j,i) << " " << j << " " << i << std::endl;
            }
        }
        std::cout << paths << std::endl;
    }

    public:
    BM(double mu, double sigma, int t, int m) {
        this->Mu = mu;
        this->Sigma = sigma;
        this->T = t;
        this->M = m;

        sim_path();
    }
};

int
main(int argc, char **argv)
{
    int i;
    if (argc != 5) {
        fprintf(stderr, "usage: %s <DRIFT_COEFFICENT> <VOLATILITY> <YEARS> <N_SIMS>\n", argv[0]);
        return -1;
    }

    for (i=1; i < argc; i++) {
        if (!atof(argv[i])) {
            fprintf(stderr, "argument not valid: %s\n", argv[i]);
            return -1;
        }
    }

    BM(atof(argv[1]), atof(argv[2]), atoi(argv[3]), atoi(argv[4])); 
}
