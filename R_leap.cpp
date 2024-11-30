#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <random>
#include <chrono>
#include <numeric>  // For std::accumulate

struct rp {
    double npar;  // number of physical particles represented by the swarm
    double mass;  // mass of physical particles represented by the swarm
};

// Constants
const int ntot = 200;           // Total number of representative particles (swarms)
const double dens = 1.0;        // Dust volume mass density
const double mswrm = 1e21;      // Mass of the swarm
const double m0 = 1.0;          // Mass of the monomer
const double tend = 20.0;       // Maximum simulation time
const double dtout = 4.0;       // Time step for the output
const int nbins = 100;          // Number of bins for output histograms
const int R = 5;                // Number of events per leap

int main() {
    std::vector<rp> swarm(ntot);  // List of the swarms
    std::vector<std::vector<double> > colrates(ntot, std::vector<double>(ntot)); // Collision rate matrix
    std::vector<double> colrates_rp(ntot, 0.0);                                 // Representative particle rates
    std::vector<double> colrates_old_col(ntot);                                 // Old collision rates
    std::vector<double> m2fm(nbins, 0.0);                                       // Mass distribution function
    std::vector<double> mgrid(nbins + 1);                                       // Mass grid
    double vol = ntot * mswrm / dens;                                           // Cell volume
    double sim_time = 0.0;
    double tout = dtout;
    double total_time = 0.0;

    // Initialize mass grid
    double nord = log10(mswrm * ntot / m0);  // Total order of magnitude range
    double ll = nord / nbins;                // Logarithmic bin width
    double lll = 0.0;

    for (int i = 0; i <= nbins; ++i) {
        mgrid[i] = m0 * pow(10.0, lll);
        lll += ll;
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

    // Initialize collision rates for each representative particle
    for (int i = 0; i < ntot; ++i) {
        colrates_rp[i] = std::accumulate(colrates[i].begin(), colrates[i].end(), 0.0);
    }

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    auto start = std::chrono::high_resolution_clock::now();

    // Main simulation loop with R-leaping
    int output_count = 1;
    while (sim_time < tend) {
        double totrate = std::accumulate(colrates_rp.begin(), colrates_rp.end(), 0.0);

        // Determine the time step for R collisions
        double dtime = -log(distribution(generator)) / totrate;
        sim_time += dtime;

        // Event processing loop
        for (int r_event = 0; r_event < R; ++r_event) {
            // Select representative particle
            double random_value = distribution(generator) * totrate;
            double cumulative_rate = 0.0;
            int nri = 0;

            while (cumulative_rate < random_value && nri < ntot) {
                cumulative_rate += colrates_rp[nri++];
            }
            nri--;

            // Select physical particle
            random_value = distribution(generator) * colrates_rp[nri];
            cumulative_rate = 0.0;
            int nrk = 0;

            while (cumulative_rate < random_value && nrk < ntot) {
                cumulative_rate += colrates[nri][nrk++];
            }
            nrk--;

            // Perform the collision
            swarm[nri].mass += swarm[nrk].mass;
            swarm[nri].npar = mswrm / swarm[nri].mass;

            // Update collision rates
            for (int i = 0; i < ntot; ++i) {
                colrates_old_col[i] = colrates[i][nri];
                colrates[i][nri] = swarm[i].npar * 0.5 * (swarm[i].mass + swarm[nri].mass) / vol;
                colrates[nri][i] = colrates[i][nri];
                colrates_rp[i] += colrates[i][nri] - colrates_old_col[i];
            }
            colrates_rp[nri] = std::accumulate(colrates[nri].begin(), colrates[nri].end(), 0.0);
        }

        // Output at intervals
        if (sim_time >= tout) {
            std::cout << "Time: " << sim_time << ", producing output " << output_count << "\n";

            // Calculate histogram
            m2fm.assign(nbins, 0.0);
            for (int i = 0; i < ntot; ++i) {
                int j = 0;
                while (j < nbins && swarm[i].mass >= mgrid[j]) {
                    ++j;
                }
                if (j < nbins) {
                    m2fm[j] += (swarm[i].npar * pow(swarm[i].mass, 2)) / ((mgrid[j + 1] - mgrid[j]) * mswrm * ntot);
                }
            }

            // Write to file
            std::string fname = "R-" + std::to_string(output_count) + ".dat";
            std::ofstream file(fname);
            if (!file) {
                std::cerr << "Error opening file: " << fname << std::endl;
            } else {
                for (int j = 0; j < nbins; ++j) {
                    file << sqrt(mgrid[j + 1] * mgrid[j]) << " " << m2fm[j] << "\n";
                }
                file.close();
            }

            output_count++;
            tout += dtout;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    std::cout << "Simulation runtime: " << duration.count() << " seconds\n";
    return 0;
}
