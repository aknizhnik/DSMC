#include <iostream>
#include <chrono>

#include "dsmc.h"

int main_1D() {
    std::array<double,3> Ls{1.0, 1.0, 0.9};
    std::array<int,3> Ns{1, 1, 90};
    int npoint_cell = 1000;
    const double Temp = 300.0;
    const double Pres = 1.0 * 1e6 / 760;  // 1 Torr in CGS
    const double Vel0 = sqrt(3*kB*Temp/(32.0*AUM));
    const double Nden = Pres /(kB * Temp);
    const double grad_P = 0.01 * Pres / 1.0;
    double cell_vol = Ls[0] * Ls[1] * Ls[2] / ( Ns[0] * Ns[1] * Ns[2]);
    double Neff = Nden * cell_vol / npoint_cell;
    double dt = 0.3 * Ls[2] / Ns[2] / Vel0;
    const int nstep = 100000;


    DSMC_Wall zbound(0.5, Temp);

    DSMC_Molecule molec1 = {1e-15, 32.0};
    std::vector<DSMC_Molecule> molecules;
    molecules.push_back(molec1);
    std::vector<double> Neffs;
    Neffs.push_back(Neff);

    DSMC_1D_System syst( Ls,  Ns, npoint_cell, molecules, Neffs, &zbound);

    double rho_avg;
    Vec3 vel_avg;
    syst.init_pos();
    syst.init_vel(Temp);
    syst.calc_average(rho_avg, vel_avg);
    std::cout << "Vel = " << vel_avg.transpose() << "\n";

    auto t1 = std::chrono::high_resolution_clock::now();
    for(int istep=0 ; istep < nstep ; istep++) {
        syst.integrate(dt);
        syst.calc_boundary();
        syst.update_distr();
        syst.collide(dt);
        syst.calc_average(rho_avg, vel_avg);
        syst.set_pres_grad(dt, grad_P);

        if(istep % 1000 == 0) {
            std::cout << "Istep = " << istep << "  Vel = " << vel_avg.transpose() << "\n";
            auto t2 = std::chrono::high_resolution_clock::now();
            auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
            std::cout << "Time = " << dt << " ms\n";
            t1 = std::chrono::high_resolution_clock::now();
        }
    }
    syst.print_distr("distr.dat", rho_avg, vel_avg);

    return 0;
}

int main_1D_tun() {
    std::array<double,3> Ls{1.0, 1.0, 0.1};
    std::array<int,3> Ns{1, 1, 30};
    int npoint_cell = 3000;
    const double Temp = 300.0;
    const double Pres = 1.0 * 1e6 / 760;  // 1 Torr in CGS
    const double Vel0 = sqrt(3*kB*Temp/(32.0*AUM));
    const double Nden = Pres /(kB * Temp);
    const double grad_P = 0.5 * Pres / 20.0;
    double cell_vol = Ls[0] * Ls[1] * Ls[2] / ( Ns[0] * Ns[1] * Ns[2]);
    double Neff = Nden * cell_vol / npoint_cell;
    double dt = 0.3 * Ls[2] / Ns[2] / Vel0;
    const int nstep = 100000;


    DSMC_Wall zbound(0.0, Temp);

    DSMC_Molecule molec1 = {2e-15, 32.0};
    std::vector<DSMC_Molecule> molecules;
    molecules.push_back(molec1);
    std::vector<double> Neffs;
    Neffs.push_back(Neff);

    DSMC_1D_System syst( Ls,  Ns, npoint_cell, molecules, Neffs, &zbound);

    double rho_avg;
    Vec3 vel_avg;
    syst.init_pos();
    syst.init_vel(Temp);
    syst.calc_average(rho_avg, vel_avg);
    std::cout << "Vel = " << vel_avg.transpose() << "\n";

    auto t1 = std::chrono::high_resolution_clock::now();
    for(int istep=0 ; istep < nstep ; istep++) {
        syst.integrate(dt);
        syst.calc_boundary();
        syst.update_distr();
        syst.collide(dt);
        syst.calc_average(rho_avg, vel_avg);
        syst.set_pres_grad(dt, grad_P);

        if(istep % 1000 == 0) {
            std::cout << "Istep = " << istep << "  Vel = " << vel_avg.transpose() << "\n";
            auto t2 = std::chrono::high_resolution_clock::now();
            auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
            std::cout << "Time = " << dt << " ms\n";
            t1 = std::chrono::high_resolution_clock::now();
        }
    }
    syst.print_distr("distr.dat", rho_avg, vel_avg);

    return 0;
}


int main_2D() {
    std::array<double,3> Ls{1.0, 0.1, 0.1};
    std::array<int,3> Ns{1, 10, 10};
    int npoint_cell = 1000;
    const double Temp = 300.0;
    const double Pres = 0.1 * 1e6 / 760;  // 1 Torr in CGS
    const double Vel0 = sqrt(3*kB*Temp/(32.0*AUM));
    const double Nden = Pres /(kB * Temp);
    const double grad_P = 0.5 * Pres / 20.0;
    double cell_vol = Ls[0] * Ls[1] * Ls[2] / ( Ns[0] * Ns[1] * Ns[2]);
    double Neff = Nden * cell_vol / npoint_cell;
    double dt = 0.3 * Ls[2] / Ns[2] / Vel0;
    const int nstep = 50000;
    std::ofstream fres("results.dat");


    DSMC_Wall zbound(0.5, Temp);
    DSMC_Wall ybound(0.5, Temp);

    DSMC_Molecule molec1 = {2.0e-15, 32.0};
    std::vector<DSMC_Molecule> molecules;
    molecules.push_back(molec1);
    std::vector<double> Neffs;
    Neffs.push_back(Neff);

    DSMC_2D_System syst( Ls,  Ns, npoint_cell, molecules, Neffs, &ybound, &zbound);

    double rho_avg;
    Vec3 vel_avg;
    syst.init_pos();
    syst.init_vel(Temp);
    syst.calc_average(rho_avg, vel_avg);
    std::cout << "Vel = " << vel_avg.transpose() << "\n";

    auto t1 = std::chrono::high_resolution_clock::now();
    for(int istep=0 ; istep < nstep ; istep++) {
        syst.integrate(dt);
        syst.calc_boundary();
        syst.update_distr();
        syst.collide(dt);
        syst.calc_average(rho_avg, vel_avg);
        syst.set_pres_grad(dt, grad_P);

        if(istep % 1000 == 0) {
            std::cout << "Istep = " << istep << "  Vel = " << vel_avg.transpose() << "\n";
            fres << "Istep = " << istep << "  Vel = " << vel_avg.transpose() << "\n";
            auto t2 = std::chrono::high_resolution_clock::now();
            auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
            std::cout << "Time = " << dt << " ms\n";
            t1 = std::chrono::high_resolution_clock::now();
        }
    }
    syst.print_distr("distr_2D.dat", rho_avg, vel_avg);

    return 0;
}

int main(){
    return main_2D();
}
