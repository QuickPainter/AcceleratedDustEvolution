#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <random>

struct rp {
    double npar;  // number of physical particles represented by the swarm
    double mass;  // mass of physical particles represented by the swarm
};

// Constants (equivalent to Fortran `parameter`)
const int ntot = 200;          // total number of representative particles (swarms)
const double dens = 1.0f;       // dust volume mass density
const double mswrm = 10.e20f;   // mass of the swarm
const double m0 = 1.0f;         // mass of the monomer (initial mass of the particles)
const double tend = 20.0f;      // maximum time of the simulation
const double dtout = 4.0f;      // time step for the output
const int nbins = 100;         // number of bins for output histograms

int main() {
    std::vector<rp> swarm(ntot);                               // List of the swarms
    int nri, nrk;
    double vol = ntot * mswrm / dens;                           // Volume of the cell
    double sim_time = 0.0f;                                     // Renamed from `time`
    double dtime;
    double tout = dtout;
    std::vector<std::vector<double> > colrates(ntot, std::vector<double>(ntot)); // Collision rate matrix
    std::vector<double> colrates_rp(ntot);                      // Probability of choosing representative particles
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
    for (i = 0; i < ntot; ++i) {
        swarm[i].npar = mswrm / m0;
        swarm[i].mass = m0;
    }

    // Calculate collision rates
    for (i = 0; i < ntot; ++i) {
        for (j = 0; j < ntot; ++j) {
            colrates[i][j] = swarm[i].npar * 0.5f * (swarm[i].mass + swarm[j].mass) / vol;
        }
    }

    for (i = 0; i < ntot; ++i) {
        colrates_rp[i] = 0.0f;
        for (j = 0; j < ntot; ++j) {
            colrates_rp[i] += colrates[i][j];
        }
    }

    // Main loop
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    while (sim_time < tend) {
        // Calculate total rate for collisions
        totrate = 0.0f;
        for (const auto& rate : colrates_rp) {
            totrate += rate;
        }

        // Determine the time step for the next collision
        random_value = distribution(generator);
        dtime = -1.0f / totrate * log(random_value);

        // Update time
        sim_time += dtime;

        // Select representative particle `nri` for collision
        random_value = distribution(generator) * totrate;
        fin = colrates_rp[0];
        j = 0;
        while (random_value > fin && j < ntot - 1) {
            j++;
            fin += colrates_rp[j];
        }
        nri = j;

        // Select physical particle `nrk` for collision
        random_value = distribution(generator) * colrates_rp[nri];
        fin = colrates[nri][0];
        j = 0;
        while (random_value > fin && j < ntot - 1) {
            j++;
            fin += colrates[nri][j];
        }
        nrk = j;

        // Perform the collision and update the representative particle `nri`
        swarm[nri].mass += swarm[nrk].mass;
        swarm[nri].npar = mswrm / swarm[nri].mass;

        // Update collision rates matrix
        for (i = 0; i < ntot; ++i) {
            colrates_old_col[i] = colrates[i][nri];
        }

        for (i = 0; i < ntot; ++i) {
            colrates[i][nri] = swarm[i].npar * 0.5f * (swarm[i].mass + swarm[nri].mass) / vol;
            colrates_rp[i] += (colrates[i][nri] - colrates_old_col[i]);
        }
        colrates_rp[nri] = 0.0f;
        for (j = 0; j < ntot; ++j) {
            colrates_rp[nri] += colrates[nri][j];
        }

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

            fname = "out-" + std::to_string(k) + ".dat";
            std::ofstream file(fname);
            for (j = 0; j < nbins; ++j) {
                file << sqrt(mgrid[j + 1] * mgrid[j]) << " " << m2fm[j] << "\n";
            }
            file.close();
            k++;
            tout += dtout;
        }
    }

    std::cout << "tend exceeded, finishing simulation\n";
    return 0;
}
