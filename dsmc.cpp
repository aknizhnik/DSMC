//
// Created by Andrey Knizhnik on 24/04/2021.
//

#include <random>
#include <iostream>

#include "dsmc.h"

std::random_device rd{};
std::mt19937 gen{rd()};
std::uniform_real_distribution<double> uni(0,1);
std::normal_distribution<> dnorm{0.0,1.0};




bool DSMC_Bounadary::check(const DSMC_Point& point, double& tau){
    Vec3 rd = point.pos - pos;
    double rd_normal = normal.dot(rd);
    if(rd_normal < 0){
        double vel_normal = normal.dot(point.vel);
        tau = rd_normal / vel_normal;
        return true;
    }
    else{
        return false;
    }
}

void DSMC_Period::strike(DSMC_Point &point, double molmass) {
    point.pos += normal;
}

void DSMC_Wall::strike(DSMC_Point &point, double molmass) {
    // find intersection point
    Vec3 pos_old = point.pos;
    Vec3 vel_old = point.vel;

    Vec3 rd = point.pos - pos;
    double rd_normal = normal.dot(rd);
    double vel_normal = normal.dot(point.vel);
    double tau = rd_normal / vel_normal;
    point.pos -= 1.00001*tau * point.vel;

    // select spectral or diffusive scattering
    double r = uni(gen);
    if(r < reflectance){
        point.vel -= 2.0* vel_normal * normal;
        point.pos += tau * point.vel;
    }
    else{// diffusive scattering
        double vel0 = sqrt(kB*Temp/(molmass*AUM));
        point.vel = vel0 * dnorm(gen) * axes1 + vel0 * dnorm(gen) * axes2 +
                vel0 * sqrt(-2.0*log(uni(gen))) * normal;
        point.pos += tau * point.vel;
    }
}

// ===============================================================================

DSMC_Cell::DSMC_Cell(DSMC_System *psystem) : psystem(psystem) {
    int nmolec = psystem->molecules.size();
    points.resize(nmolec);
    for(int i=0 ; i<nmolec ; i++){
        points[i].resize(psystem->npoint_cell);
    }
    Neffs = psystem->Neffs;
};


void DSMC_Cell::init_vel(double Temp) {
    int nmolec = psystem->molecules.size();
    for(int i=0 ; i<nmolec ; i++) {
        double molmass = psystem->molecules[i].molmass;
        double vel0 = sqrt(kB*Temp/(molmass*AUM));
        for(auto& point : points[i]){
            for(int m=0 ; m<3 ; m++){
                point.vel[m] = vel0*dnorm(gen);
            }
        }
    }
}

void DSMC_Cell::init_pos(const Vec3& r_min, const Vec3& r_max){
    this->r_min = r_min;
    this->r_max = r_max;
    for(auto& molec_points : points){
        for(auto& point : molec_points) {
            for(int m=0 ; m<3 ; m++) {
                double r = uni(gen);
                point.pos[m] = r_min[m] + r * (r_max[m] - r_min[m]);
            }
        }
    }
}

void DSMC_Cell::integrate(double dt){
    for(auto& molec_points : points){
        for(auto& point : molec_points) {
            point.pos += point.vel * dt;
        }
    }
}

void DSMC_Cell::calc_boundary(std::vector<DSMC_Bounadary*>& pboundaries) {
    int nmolec = psystem->molecules.size();
    for(int m=0 ; m<nmolec ; m++){
        double molmass = psystem->molecules[m].molmass;
        for (auto &point : points[m]) {
            bool flag;
            int iter = 0;
//            std::vector<Vec3> positions;
//            positions.push_back(point.pos);
//            std::vector<Vec3> velocities;
//            velocities.push_back(point.vel);
//            std::vector<int> ibounds;
//            ibounds.push_back(-1);
            do {
                flag = true;
                double tau_max = 0.0;
                int boundary_max = -1;
                int nboundaries = pboundaries.size();
                for (int ib = 0; ib < nboundaries; ib++) {
                    double tau;
                    if (pboundaries[ib]->check(point, tau)) {
                        if (tau > tau_max) {
                            tau_max = tau;
                            boundary_max = ib;
                            flag = false;
                        }
                    }
                }
                if (!flag) {
                    pboundaries[boundary_max]->strike(point, molmass);
//                    positions.push_back(point.pos);
//                    velocities.push_back(point.vel);
//                    ibounds.push_back(boundary_max);
                }
                iter++;
//                if (iter > 30) {
//                    for(int k=0 ; k<positions.size(); k++){
//                        std::cout << ibounds[k] << "  " << positions[k].transpose() << "  "
//                        << velocities[k].transpose() << std::endl;
//                    }
//                    iter=iter;
//                }
            }while(!flag);
            if((point.pos[0] < 0) || (point.pos[1] < 0) || (point.pos[2] < 0) || (point.pos[1] > 0.1) ||
                    (point.pos[2] > 0.1)){
                point.pos[0] = point.pos[0];
            }
        }
    }
}

void DSMC_Cell::update_escape(std::vector<std::vector<DSMC_Point>>& escaped_points, std::vector<double>& Neffs_cell){
    int nmolec = psystem->molecules.size();
    escaped_points.resize(nmolec);
    Neffs_cell = Neffs;

    for(int m=0 ; m<nmolec ; m++) {
        std::vector<DSMC_Point> new_points;
        new_points.reserve(points[m].size());
        escaped_points[m].clear();
        for (auto &point : points[m]){
            bool flag = true;
            for(int k=0 ; k<3 ; k++){
                if((point.pos[k] > r_max[k])  || (point.pos[k] < r_min[k]) ){
                    flag = false;
                    break;
                }
            }
            if(flag){
                new_points.push_back(point);
            }
            else{
                escaped_points[m].push_back(point);
            }
        }
        std::swap(points[m], new_points);
    }
}

void DSMC_Cell::add_point(int molec_ind, DSMC_Point& point, double Neff){
    int npoint = points[molec_ind].size();

    points[molec_ind].push_back(point);
    //rescale Neff for this molecule type
    //Neffs[molec_ind] = (npoint*Neffs[molec_ind] + Neff) / (npoint + 1);
}

void DSMC_Cell::collide(double dt){
    int nmolec = psystem->molecules.size();
    double cell_volume = psystem->cell_volume;
    double V_max = psystem->V_max;

    for(int m=0 ; m<nmolec ; m++){
        int N_m = points[m].size();
        double sigma_m = psystem->molecules[m].sigma;
        for(int n=0 ; n<nmolec ; n++) {
            int N_n = points[n].size();
            double sigma_n = psystem->molecules[n].sigma;
            double sigma_mn = 0.5 * (sigma_m + sigma_n);
            //number of collisions for the given pair of molecule types
            int M_mn_cand = 0.5 * N_m * N_n * Neffs[n] / cell_volume * psystem->V_max * sigma_mn * dt;

            for (int p = 0; p < M_mn_cand; p++) {
                int point_m_i = N_m * uni(gen);
                int point_n_j = N_n * uni(gen);
                Vec3 V_ij = points[m][point_m_i].vel - points[n][point_n_j].vel;
                double V_ij_norm = V_ij.norm();
                double r = uni(gen);
                if (r < V_ij_norm / V_max) {
                    calc_collision( m, point_m_i, n, point_n_j, V_ij_norm);
                }
            }
        }
    }
}

void DSMC_Cell::calc_collision(int molec_m, int point_m_i, int molec_n, int point_n_j, double V_ij){
    Vec3 V_cm = 0.5*(points[molec_m][point_m_i].vel + points[molec_n][point_n_j].vel);
    double phi = 2.0 * M_PI * uni(gen);
    double cos_theta = 2.0*uni(gen) - 1;
    double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
    Vec3 Vr {sin_theta*cos(phi), sin_theta*sin(phi), cos_theta};
    Vr *= V_ij;
    points[molec_m][point_m_i].vel = V_cm + 0.5* Vr;
    points[molec_n][point_n_j].vel = V_cm - 0.5* Vr;
}

void DSMC_Cell::set_pres_grad(double dt, double grad_p){

    for(auto& molec_points : points){
        for(auto& point : molec_points) {
            point.vel[0] -= grad_p / rho_avg * dt;
        }
    }
}

void  DSMC_Cell::calc_average(double& avg_rho, Vec3& avg_vel){
    int nmolec = psystem->molecules.size();
    double rho = 0.0;
    Vec3 vel = {0.0, 0.0, 0.0};
    for(int m=0 ; m<nmolec ; m++){
        rho += Neffs[m] * points[m].size() * psystem->molecules[m].molmass * AUM;
        for(auto& point : points[m]) {
            vel += point.vel/points[m].size();
        }
    }
    rho /= psystem->cell_volume;
    rho_sum += rho;
    vel_sum += vel;
    num_avg++;
    avg_rho = rho_sum / num_avg;
    rho_avg = avg_rho;
    avg_vel = vel_sum / num_avg;
    vel_avg = avg_vel;
}

void  DSMC_Cell::get_average(double& avg_rho, Vec3& avg_vel){
    avg_rho = rho_avg;
    avg_vel = vel_avg;
}




void DSMC_System::init_vel(double Temp) {
    for(auto& cell : cells){
        cell.init_vel(Temp);
    }
}

void DSMC_System::init_pos() {
    Vec3 r_min, r_max;

    for(int i_x = 0 ; i_x<Ns[0] ; i_x++){
        r_min[0] = Ls[0] / Ns[0] * i_x;
        r_max[0] = Ls[0] / Ns[0] * (i_x + 1);
        for(int i_y = 0 ; i_y<Ns[1] ; i_y++) {
            r_min[1] = Ls[1] / Ns[1] * i_y;
            r_max[1] = Ls[1] / Ns[1] * (i_y + 1);
            for (int i_z = 0; i_z < Ns[2]; i_z++) {
                r_min[2] = Ls[2] / Ns[2] * i_z;
                r_max[2] = Ls[2] / Ns[2] * (i_z + 1);
                int cell_ind = get_cell_index(i_x, i_y, i_z);
                cells[cell_ind].init_pos(r_min, r_max);
            }
        }
    }
}

void DSMC_System::integrate(double dt){
    size_t ncells = cells.size();
#pragma omp parallel for
    for(int i=0; i<ncells  ; i++){
        cells[i].integrate(dt);
    }
}

void DSMC_System::calc_boundary(){
    std::map<int, std::vector<DSMC_Bounadary*>>::iterator itr;
    for (itr = boundary_cells.begin(); itr != boundary_cells.end(); ++itr){
        // boundary cells for the given boundary
        int cell_ind = itr->first;
        std::vector<DSMC_Bounadary*> pboundaries = itr->second;
        cells[cell_ind].calc_boundary(pboundaries);
    }
}

void DSMC_System::update_distr(){
    int ncell = cells.size();
    int nmolec = molecules.size();
    std::vector<std::vector<std::vector<DSMC_Point>>> escaped_points(ncell);
    std::vector<std::vector<double>> cell_Neffs(ncell);

#pragma omp parallel for
    for(int n=0 ; n<ncell ; n++){
        cells[n].update_escape(escaped_points[n], cell_Neffs[n]);
    }

    // redistribute molecules
    for(int n=0 ; n<ncell ; n++){
        for(int m=0 ; m<nmolec ; m++){
            int nescaped = escaped_points[n][m].size();
            for(int l=0 ; l<nescaped ; l++){
                int cell_ind = get_cell_index(escaped_points[n][m][l].pos);
                if(cell_ind < 0 || cell_ind >= ncell){
                    Vec3 pos = escaped_points[n][m][l].pos;
                    int cell_ind = get_cell_index(escaped_points[n][m][l].pos);
                }
                cells[cell_ind].add_point(m, escaped_points[n][m][l], cell_Neffs[n][m]);

            }
        }
    }
}

void DSMC_System::collide(double dt){
    size_t ncells = cells.size();
#pragma omp parallel for
    for(int i=0; i<ncells  ; i++){
        cells[i].collide(dt);
    }
}

void DSMC_System::set_pres_grad(double dt, double grad_p){
    size_t ncells = cells.size();
#pragma omp parallel for
    for(int i=0; i<ncells  ; i++){
        cells[i].set_pres_grad(dt, grad_p);
    }
}

void DSMC_System::calc_average(double &rho, Vec3 &vel) {
    rho = 0.0;
    vel.setZero();
    double rho_cell;
    Vec3 vel_cell;
    for (auto& cell : cells){
        cell.calc_average(rho_cell, vel_cell);
        rho += rho_cell;
        vel += vel_cell;
    }
    rho /= cells.size();
    vel /= cells.size();
}

void DSMC_System::print_distr(const std::string& filename, double& avg_rho, Vec3& avg_vel) {
    std::ofstream fdistr(filename);
    avg_rho = 0.0;
    avg_vel.setZero();
    double rho_cell;
    Vec3 vel_cell;
    for(int i_x = 0 ; i_x<Ns[0] ; i_x++){
        for(int i_y = 0 ; i_y<Ns[1] ; i_y++) {
            for (int i_z = 0; i_z < Ns[2]; i_z++) {
                int cell_ind = get_cell_index(i_x, i_y, i_z);
                cells[cell_ind].get_average(rho_cell, vel_cell);
                avg_rho += rho_cell;
                avg_vel += vel_cell;
                fdistr << i_x <<  "  " << i_y << "  " << i_z << "  " << rho_cell << "  " <<
                vel_cell.transpose() << "\n";
            }
        }
    }
    avg_rho /= cells.size();
    avg_vel /= cells.size();
}


DSMC_System::DSMC_System(const std::array<double,3>& Ls, const std::array<int,3>& Ns, int npoint_cell,
                         std::vector<DSMC_Molecule>& molecules,
                         std::vector<double>& Neffs,
                         DSMC_Bounadary* yboundary, DSMC_Bounadary* zboundary) :
        Ls(Ls), Ns(Ns), npoint_cell(npoint_cell), molecules(molecules), Neffs(Neffs) {
    this->cell_volume = (Ls[0] * Ls[1] * Ls[2]) / (Ns[0] * Ns[1] * Ns[2]);

    //periodic boundary conditions along OX
    boundaries.push_back(new DSMC_Period({0.0, 0.0, 0.0}, {Ls[0], 0.0, 0.0}));
    boundaries.push_back(new DSMC_Period({Ls[0], 0.0, 0.0}, {-Ls[0], 0.0, 0.0}));

    DSMC_Wall *pbound = dynamic_cast<DSMC_Wall *>(yboundary);
    if (pbound) {
        auto pbound1 = new DSMC_Wall(*pbound);
        pbound1->pos = {0.0, 0.0, 0.0};
        pbound1->normal = {0.0, 1.0, 0.0};
        pbound1->set_axes();
        boundaries.push_back(pbound1);
        auto pbound2 = new DSMC_Wall(*pbound);
        pbound2->pos = {0.0, Ls[1], 0.0};
        pbound2->normal = {0.0, -1.0, 0.0};
        pbound2->set_axes();
        boundaries.push_back(pbound2);
    } else {
        boundaries.push_back(new DSMC_Period({0.0, 0.0, 0.0}, {0.0, Ls[1], 0.0}));
        boundaries.push_back(new DSMC_Period({0.0, Ls[1], 0.0}, {0.0, -Ls[1], 0.0}));
    }

    pbound = dynamic_cast<DSMC_Wall *>(zboundary);
    if (pbound) {
        auto pbound1 = new DSMC_Wall(*pbound);
        pbound1->pos = {0.0, 0.0, 0.0};
        pbound1->normal = {0.0, 0.0, 1.0};
        pbound1->set_axes();
        boundaries.push_back(pbound1);
        auto pbound2 = new DSMC_Wall(*pbound);
        pbound2->pos = {0.0, 0.0, Ls[2]};
        pbound2->normal = {0.0, 0.0, -1.0};
        pbound2->set_axes();
        boundaries.push_back(pbound2);
    } else {
        boundaries.push_back(new DSMC_Period({0.0, 0.0, 0.0}, {0.0, 0.0, Ls[2]}));
        boundaries.push_back(new DSMC_Period({0.0, 0.0, Ls[2]}, {0.0, 0.0, -Ls[2]}));
    }
}



DSMC_1D_System::DSMC_1D_System(const std::array<double,3>& Ls, const std::array<int,3>& Ns, int npoint_cell,
                               std::vector<DSMC_Molecule>& molecules,
                               std::vector<double>& Neffs,
                               DSMC_Bounadary* zboundary) : DSMC_System( Ls, Ns, npoint_cell, molecules, Neffs,
                                                                         new DSMC_Period(), zboundary){
    //number of cells
    int ncells = Ns[2];
    for (int i = 0; i < ncells; i++) {
        cells.emplace_back(DSMC_Cell(this));
    }

    // set boundary cells
    for(int i_z=0 ; i_z<Ns[2] ; i_z++) {
        boundary_cells[i_z].push_back(boundaries[0]);
        boundary_cells[i_z].push_back(boundaries[1]);
        //boundary_cells.insert(std::pair<int, DSMC_Bounadary*>( i_z, boundaries[0]));
        //boundary_cells.insert(std::pair<int, DSMC_Bounadary*>( i_z, boundaries[1]));
    }
    for(int i_z=0 ; i_z<Ns[2] ; i_z++) {
        boundary_cells[i_z].push_back(boundaries[2]);
        boundary_cells[i_z].push_back(boundaries[3]);
        //boundary_cells.insert(std::pair<int, DSMC_Bounadary*>( i_z, boundaries[2]));
        //boundary_cells.insert(std::pair<int, DSMC_Bounadary*>( i_z, boundaries[3]));
    }
    boundary_cells[0].push_back(boundaries[4]);
    boundary_cells[Ns[2]-1].push_back(boundaries[5]);
    //boundary_cells.insert(std::pair<int, DSMC_Bounadary*>( 0, boundaries[4]));
    //boundary_cells.insert(std::pair<int, DSMC_Bounadary*>( Ns[2]-1, boundaries[5]));
}

DSMC_2D_System::DSMC_2D_System(const std::array<double,3>& Ls, const std::array<int,3>& Ns, int npoint_cell,
                               std::vector<DSMC_Molecule>& molecules,
                               std::vector<double>& Neffs,
               DSMC_Bounadary* yboundary, DSMC_Bounadary* zboundary) : DSMC_System( Ls, Ns, npoint_cell, molecules,
                                                                                    Neffs, yboundary, zboundary) {
    //gen.seed(1234);

    //number of cells
    int ncells = Ns[1] * Ns[2];
    for (int i = 0; i < ncells; i++) {
        cells.emplace_back(DSMC_Cell(this));
    }

    // set boundary cells
    for(int i_y=0 ; i_y<Ns[1] ; i_y++){
        for(int i_z=0 ; i_z<Ns[2] ; i_z++) {
            int cell_ind = get_index(i_y,i_z);
            boundary_cells[cell_ind].push_back(boundaries[0]);
            boundary_cells[cell_ind].push_back(boundaries[1]);
            boundary_cells[cell_ind].push_back(boundaries[2]);
            boundary_cells[cell_ind].push_back(boundaries[3]);
            boundary_cells[cell_ind].push_back(boundaries[4]);
            boundary_cells[cell_ind].push_back(boundaries[5]);
            //boundary_cells.insert(std::pair<int, DSMC_Bounadary*>( cell_ind, boundaries[0]));
            //boundary_cells.insert(std::pair<int, DSMC_Bounadary*>( cell_ind, boundaries[1]));
        }
    }
    /*
    for(int i_z=0 ; i_z<Ns[2] ; i_z++) {
        int cell_ind = get_index(0,i_z);
        boundary_cells[cell_ind].push_back(boundaries[2]);
        //boundary_cells.insert(std::pair<int, DSMC_Bounadary*>( cell_ind, boundaries[2]));
        cell_ind = get_index(Ns[1]-1,i_z);
        boundary_cells[cell_ind].push_back(boundaries[3]);
        //boundary_cells.insert(std::pair<int, DSMC_Bounadary*>( cell_ind, boundaries[3]));

    }
    for(int i_y=0 ; i_y<Ns[1] ; i_y++){
        int cell_ind = get_index(i_y,0);
        boundary_cells[cell_ind].push_back(boundaries[4]);
        //boundary_cells.insert(std::pair<int, DSMC_Bounadary*>( cell_ind, boundaries[4]));
        cell_ind = get_index(i_y,Ns[2]-1);
        boundary_cells[cell_ind].push_back(boundaries[5]);
        //boundary_cells.insert(std::pair<int, DSMC_Bounadary*>( cell_ind, boundaries[5]));
    }
     */
}
