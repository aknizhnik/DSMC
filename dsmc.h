//
// Created by Andrey Knizhnik on 24/04/2021.
//

#ifndef DSMC_2D_DSMC_2D_H
#define DSMC_2D_DSMC_2D_H

#include <vector>
#include <array>
#include <map>
#include <iterator>
#include <fstream>

#include "Eigen/Dense"

inline constexpr double kB = 1.38e-16;
inline constexpr double AUM = 1.6e-24;


using Vec3 = Eigen::Vector3d;

class DSMC_Molecule{
public:
    double sigma;
    double molmass;
};

class DSMC_Point{
public:
    DSMC_Point(){};
    DSMC_Point(const Vec3& pos, const Vec3& vel) : pos(pos), vel(vel)
    {};

    Vec3 pos;
    Vec3 vel;
};

class DSMC_System;

// ===================================================================================

class DSMC_Bounadary{
public:
    Vec3 pos;
    Vec3 normal;
public:
    DSMC_Bounadary() : pos({0.0, 0.0, 0.0}), normal({1.0, 0.0, 0.0}) {};
    DSMC_Bounadary(const Vec3& pos, const Vec3& normal) : pos(pos), normal(normal) {};
    virtual ~DSMC_Bounadary(){};

    bool check(const DSMC_Point& point, double& tau);

    virtual void strike(DSMC_Point& point, double molmass)=0;
};

class DSMC_Wall : public DSMC_Bounadary{
protected:
    Vec3 axes1, axes2;  //!< orthogonal axes in the plane
public:
    DSMC_Wall(double r, double Temp) : DSMC_Bounadary(), reflectance(r), Temp(Temp){};
    DSMC_Wall(const Vec3& pos, const Vec3& normal, double r, double Temp) : DSMC_Bounadary(pos, normal),
        reflectance(r), Temp(Temp){
        set_axes();
    };

    double reflectance = 0.5;
    double Temp;

    void set_axes(){
        normal /= normal.norm();
        Vec3 n = {1.0, 1.0, 1.0};
        axes1 = this->normal.cross(n);
        axes1 /= axes1.norm();
        axes2 = this->normal.cross(axes1);
    }

    virtual void strike(DSMC_Point& point, double molmass);
};

class DSMC_Period : public DSMC_Bounadary{
public:
    DSMC_Period() : DSMC_Bounadary(){};
    DSMC_Period(const Vec3& pos, const Vec3& normal) : DSMC_Bounadary(pos, normal){};


    virtual void strike(DSMC_Point& point, double molmass);
};

// ===================================================================================


class DSMC_Cell{
protected:
    Vec3 r_min;
    Vec3 r_max;
    double rho_sum = 0.0;
    double rho_avg;
    Vec3 vel_sum = {0.0, 0.0, 0.0};
    Vec3 vel_avg;
    int num_avg = 0;
public:
    DSMC_Cell(DSMC_System* psystem);

    DSMC_System* psystem;



    std::vector<std::vector<DSMC_Point>> points;

    std::vector<double> Neffs;    //!< number of real molecules per one in DSMC model

    void init_vel(double Temp);
    void init_pos(const Vec3& r_min, const Vec3& r_max);

    void integrate(double dt);

    void calc_boundary(std::vector<DSMC_Bounadary*>& pboundaries);

    void update_escape(std::vector<std::vector<DSMC_Point>>& escaped_points, std::vector<double>& Neffs_cell);

    void add_point(int molec_ind, DSMC_Point& point, double Neff);

    void collide(double dt);

    void calc_collision(int molec_m, int point_m_i, int molec_n, int point_n_j, double V_ij);

    void set_pres_grad(double dt, double grad_p);

    void calc_average(double& rho, Vec3& vel);

    void get_average(double& rho, Vec3& vel);
};



class DSMC_System{
protected:
    std::array<double, 3> Ls;
    std::array<int, 3> Ns;
    std::vector<DSMC_Cell> cells;
    std::vector<DSMC_Molecule> molecules;
    std::vector<double> Neffs;

    int npoint_cell;
    double cell_volume = 0.0;
    std::vector<DSMC_Bounadary*> boundaries;
    std::map<int, std::vector<DSMC_Bounadary*>> boundary_cells;

    int get_cell_index(int i_x, int i_y, int i_z){
        return i_x * Ns[1]*Ns[2] + i_y * Ns[2] + i_z;
    }
    int get_cell_index(Vec3& pos){
        int i_x = (int)floor(Ns[0]*pos[0]/Ls[0]);
        int i_y = (int)floor(Ns[1]*pos[1]/Ls[1]);
        int i_z = (int)floor(Ns[2]*pos[2]/Ls[2]);
        return get_cell_index(i_x, i_y, i_z);
    }
public:
    DSMC_System(const std::array<double,3>& Ls, const std::array<int,3>& Ns, int npoint_cell,
                std::vector<DSMC_Molecule>& molecules,
                std::vector<double>& Neffs,
                DSMC_Bounadary* yboundary, DSMC_Bounadary* zboundary) ;
    ~DSMC_System(){};

    double V_max = 2e5;

    void init_vel(double Temp);

    void init_pos();

    void integrate(double dt);

    void calc_boundary();

    void update_distr();

    void collide(double dt);

    void set_pres_grad(double dt, double grad_p);

    void calc_average(double& avg_rho, Vec3& avg_vel);

    void print_distr(const std::string& filename, double& avg_rho, Vec3& avg_vel);

    friend class DSMC_Cell;
};


class DSMC_1D_System : public DSMC_System{
public:
    DSMC_1D_System(const std::array<double,3>& Ls, const std::array<int,3>& Ns, int npoint_cell,
                   std::vector<DSMC_Molecule>& molecules,
                   std::vector<double>& Neffs,
                   DSMC_Bounadary* zboundary);
    ~DSMC_1D_System(){};

    int get_index(int i_z){
        return get_cell_index(0, 0, i_z);
    }

};


class DSMC_2D_System : public DSMC_System{
public:
    DSMC_2D_System(const std::array<double,3>& Ls, const std::array<int,3>& Ns, int npoint_cell,
                   std::vector<DSMC_Molecule>& molecules,
                   std::vector<double>& Neffs,
                   DSMC_Bounadary* yboundary, DSMC_Bounadary* zboundary);
    ~DSMC_2D_System(){};

    int get_index(int i_y, int i_z){
        return get_cell_index(0, i_y, i_z);
    }


};

#endif //DSMC_2D_DSMC_2D_H
