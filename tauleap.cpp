#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

struct rp {
    double npar;  // number of physical particles represented by the swarm
    double mass;  // mass of physical particles represented by the swarm
};

// Constants
const int ntot = 200;           // Total number of representative particles (swarms)
const double dens = 1.0;        // Dust volume mass density
const double mswrm = 10e20;     // Mass of the swarm
const double m0 = 1.0;          // Mass of the monomer (initial mass of the particles)
const double tend = 20.0;       // Maximum simulation time
const double dtout = 4.0;       // Time step for the output
const int nbins = 100;          // Number of bins for output histograms

int main() {
    std::vector<rp> swarm(ntot);
    double vol = ntot * mswrm / dens;
    double sim_time = 0.0f;
    double tout = dtout;
    std::vector<std::vector<double> > colrates(ntot, std::vector<double>(ntot));
    std::vector<double> colrates_rp(ntot);
    double totrate;
    
    double random_value;                                        // Renamed from `rand`
    double fin;
    std::string fname;
    std::vector<double> colrates_old_col(ntot);                 // For optimization
    std::vector<double> m2fm(nbins);                            // Mass distribution function
    std::vector<double> mgrid(nbins + 1);                       // Mass grid
    int i, j;
    int k = 1;
    double ll, lll, nord;

    // Initialize mass bins for histograms
    mgrid[0] = m0;
    nord = log10(mswrm * ntot / m0);  // Orders of magnitude in mass for histogram
    ll = nord / nbins;
    lll = ll;

    for (i = 1; i < nbins + 1; ++i) {
        mgrid[i] = m0 * pow(10.0, lll);
        lll = ll * (i + 1);
    }
    
    // Initialize swarms
    for (int i = 0; i < ntot; ++i) {
        swarm[i].npar = mswrm / m0;
        swarm[i].mass = m0;
    }

    


    // Calculate initial collision rates
    for (int i = 0; i < ntot; ++i) {
        for (int j = 0; j < ntot; ++j) {
            colrates[i][j] = swarm[j].npar * 0.5 * (swarm[i].mass + swarm[j].mass) / vol;
        }
    }

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    auto start = std::chrono::high_resolution_clock::now();

    while (sim_time < tend) {
        // Calculate total rate
        totrate = 0.0;
        for (double rate : colrates_rp) {
            totrate += rate;
        }
        
        // Determine tau using heuristic (small enough to ensure stability)
        double tau = 0.00320603 *180.0;  // Adjust the denominator to control accuracy
        //std::cout << "tau: " << tau << "\n";

        
        // Estimate number of collisions per pair within tau
        for (int i = 0; i < ntot; ++i) {
            for (int j = 0; j < ntot; ++j) {
                double expected_collisions = colrates[i][j] * tau;
                std::poisson_distribution<int> poisson(expected_collisions);
                int collision_count = poisson(generator);

                for (int c = 0; c < collision_count; ++c) {
                    // Perform collision updates (similar to before)
                    swarm[i].mass += swarm[j].mass;
                    swarm[i].npar = mswrm / swarm[i].mass;

                    for(int p = 0; p< ntot; ++p){
                        colrates_old_col[p] = colrates[p][i];
                    }

                    for (int p = 0; p < ntot; ++p) {
                        colrates[p][i] = swarm[p].npar * 0.5f * (swarm[p].mass + swarm[i].mass) / vol;
                        colrates[i][p] = swarm[i].npar * 0.5f * (swarm[p].mass + swarm[i].mass) / vol;
                    }

                    for (int p = 0; p < ntot; ++p) {
                        colrates_rp[p] += (colrates[p][i] - colrates_old_col[p]);
                    }

                    colrates_rp[i] = 0.0f;
                    for (int p = 0; p < ntot; ++p) {
                        colrates_rp[i] += colrates[i][p];
                    }
                }
            }
        }
        
        sim_time += tau;
         // Output histogram
        if (sim_time > tout) {
            std::cout << "time: " << sim_time << ", producing output " << k << "\n";

            m2fm.assign(nbins, 0.0f);
            for (i = 0; i < ntot; ++i) {
                j = 0;
                while (j < nbins && swarm[i].mass >= mgrid[j]) {
                    j++;
                }
                m2fm[j] += (swarm[i].npar * pow(swarm[i].mass, 2)) / ((mgrid[j + 1] - mgrid[j]) * mswrm * ntot);
            }

            fname = "tau-" + std::to_string(k) + ".dat";
            std::ofstream file(fname);
            for (j = 0; j < nbins; ++j) {
                file << sqrt(mgrid[j + 1] * mgrid[j]) << " " << m2fm[j] << "\n";
            }
            file.close();
            k++;
            tout += dtout;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;


    std::cout << "Tau Simulation Runtime: " << duration.count() << " seconds" << std::endl;



    std::cout << "Simulation complete.\n";
    return 0;
}
