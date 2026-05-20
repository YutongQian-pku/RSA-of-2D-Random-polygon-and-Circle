#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include "particle.h"
#include "cell.h"
#include "voxel.h"
#include "particle_type.h"
#include "config.h"

struct RsaSummary {
    std::vector<double> avg_counts;
    double avg_packing_fraction = 0.0;
};

inline int choose_particle_type(const std::vector<ParticleType>& particle_types,
                                std::mt19937& gen,
                                std::uniform_real_distribution<>& dis) {
    double rate = dis(gen);
    double cumulative_prob = 0.0;
    for (size_t i = 0; i < particle_types.size(); ++i) {
        cumulative_prob += particle_types[i].prob;
        if (rate < cumulative_prob) return static_cast<int>(i);
    }
    return static_cast<int>(particle_types.size() - 1);
}

inline std::unique_ptr<Particle> make_trial_particle(const ParticleType& p_type,
                                                     const Point& center,
                                                     std::mt19937& gen,
                                                     std::uniform_real_distribution<>& dis) {
    if (p_type.shape == ParticleShape::CIRCLE) {
        return std::make_unique<Circle>(p_type.radius, center);
    }
    double angle = 2.0 * PI * dis(gen);
    return std::make_unique<Polygon>(p_type.vertices, angle, center, p_type.shape_id);
}

inline RsaSummary rsa_cycle(const std::vector<ParticleType>& particle_types,
                            long long area,
                            long long rsa_step,
                            int cycle_num) {
    for (const auto& t : particle_types) {
        if (t.shape == ParticleShape::CIRCLE) {
            std::cout << "Circle: radius=" << t.radius
                      << ", inner_r=" << t.inner_r
                      << ", area=" << t.area
                      << ", prob=" << t.prob << std::endl;
        } else {
            std::cout << "Polygon shape " << t.shape_id
                      << ": vertices=" << t.vertices.size()
                      << ", outer_R=" << t.radius
                      << ", inner_r=" << t.inner_r
                      << ", area=" << t.area
                      << ", prob=" << t.prob << std::endl;
        }
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::vector<std::vector<int>> counts(cycle_num, std::vector<int>(particle_types.size(), 0));
    std::vector<double> packing_fractions(cycle_num, 0.0);

    for (int cycle = 0; cycle < cycle_num; ++cycle) {
        double box_l = std::sqrt(static_cast<double>(area));
        std::vector<std::unique_ptr<Particle>> particles;

        double max_char_L = 0.0;
        for (const auto& t : particle_types) max_char_L = std::max(max_char_L, t.character_L);
        if (max_char_L <= 0.0) {
            std::cerr << "Invalid particle characteristic length." << std::endl;
            std::exit(1);
        }

        double add = 1E100;
        for (const auto& t : particle_types) add = std::min(add, t.inner_r);
        if (add <= 0.0 || add == 1E100) {
            std::cerr << "Invalid inner radius for voxel pruning." << std::endl;
            std::exit(1);
        }

        Cell* grid = nullptr;
        Cell cell(max_char_L, box_l);
        grid = new Cell[cell.cell_N];
        for (int i = 0; i < cell.cell_N; i++) {
            grid[i] = std::move(cell);
            cell = Cell(max_char_L, box_l);
        }
        initial_grid(&grid[0]);

        voxel* vox = nullptr;
        voxel vox1(max_char_L, box_l);
        vox = new voxel[vox1.voxel_N];
        for (int i = 0; i < vox1.voxel_N; i++) vox[i] = vox1;
        initial_vox(&vox[0]);

        long long try_num = 0;
        double dimensionless_t = 0.0;
        int particle_num = 0;
        double occupied_area = 0.0;

        char cha[64];
        std::sprintf(cha, "kinetics-%d.txt", cycle);
        std::ofstream answer(cha, std::ios::app);
        answer << std::fixed << std::setprecision(8);
        answer << "area\t" << area << std::endl;
        answer << "try_num\tt\tparticle_num\tlast_area\tpacking_fraction" << std::endl;

        while (true) {
            try_num++;
            int type_idx = choose_particle_type(particle_types, gen, dis);
            const auto& p_type = particle_types[type_idx];
            std::unique_ptr<Particle> insert_p = make_trial_particle(p_type, get_randloc(box_l), gen, dis);

            if (!check_over(grid, insert_p.get(), vox)) {
                particles.push_back(std::move(insert_p));
                counts[cycle][type_idx]++;
                particle_num++;
                occupied_area += p_type.area;
                packing_fractions[cycle] = occupied_area / static_cast<double>(area);

                dimensionless_t = static_cast<double>(try_num) / static_cast<double>(area);
                answer << try_num << "\t" << dimensionless_t << "\t" << particle_num << "\t"
                       << p_type.area << "\t" << packing_fractions[cycle] << std::endl;

                int id = get_cell_id(grid[0], particles.back()->center.x, particles.back()->center.y);
                auto clonedParticle = particles.back()->clone();
                grid[id].particles.push_back(std::move(clonedParticle));
                grid[id].num++;

                auto p_big = grid[id].particles.back()->enlarge(add);
                if (p_big) p_big->cal_occupy(vox, box_l);
            }

            if (try_num % 1000000 == 0) {
                std::cout << "try num=" << try_num
                          << "\tt=" << dimensionless_t
                          << "\tnum=" << particle_num
                          << "\tphi=" << packing_fractions[cycle] << std::endl;
            }

            if (try_num > rsa_step) {
                std::cout << "Finished cycle " << cycle << std::endl;
                break;
            }

            if (try_num == rsa_step) {
                answer << "final stage add" << std::endl;
                std::cout << "final stage add" << std::endl;
                for (int i = 0; i < vox1.voxel_N; i++) {
                    if (!vox[i].is_occupy) {
                        Point v1;
                        v1.x = vox[i].center.x - vox[i].length * 0.5;
                        v1.y = vox[i].center.y - vox[i].length * 0.5;
                        const long long final_trials_per_voxel = std::max(1LL, static_cast<long long>(FINAL_ADDNUM_MIN));
                        for (long long j = 0; j < final_trials_per_voxel; j++) {
                            int local_type_idx = choose_particle_type(particle_types, gen, dis);
                            const auto& local_type = particle_types[local_type_idx];
                            double randx = v1.x + vox[i].length * double(rand()) / double(RAND_MAX);
                            double randy = v1.y + vox[i].length * double(rand()) / double(RAND_MAX);
                            std::unique_ptr<Particle> local_p = make_trial_particle(local_type, Point(randx, randy), gen, dis);

                            if (!check_over(grid, local_p.get(), vox)) {
                                particles.push_back(std::move(local_p));
                                counts[cycle][local_type_idx]++;
                                particle_num++;
                                occupied_area += local_type.area;
                                packing_fractions[cycle] = occupied_area / static_cast<double>(area);

                                dimensionless_t = static_cast<double>(try_num) / static_cast<double>(area);
                                answer << try_num << "\t" << dimensionless_t << "\t" << particle_num << "\t"
                                       << local_type.area << "\t" << packing_fractions[cycle] << std::endl;

                                int id = get_cell_id(grid[0], particles.back()->center.x, particles.back()->center.y);
                                auto clonedParticle = particles.back()->clone();
                                grid[id].particles.push_back(std::move(clonedParticle));
                                grid[id].num++;

                                auto p_big = grid[id].particles.back()->enlarge(add);
                                if (p_big) p_big->cal_occupy(vox, box_l);
                            }
                        }
                        vox[i].is_occupy = true;
                    }
                }
                int remaining = 0;
                for (int i = 0; i < vox1.voxel_N; i++) if (!vox[i].is_occupy) remaining++;
                std::cout << "all_vox_num = " << vox1.voxel_N << std::endl;
                std::cout << "voxel_unoccupied_num = " << remaining << std::endl;
            }
        }

        answer.close();
        scr(box_l, particles, particle_num, cycle);

        int occupied_voxels = 0;
        for (int i = 0; i < vox1.voxel_N; i++) if (vox[i].is_occupy) occupied_voxels++;
        std::cout << "all_vox_num = " << vox1.voxel_N << std::endl;
        std::cout << "voxel_occupy_num = " << occupied_voxels << std::endl;
        std::cout << "packing_fraction = " << packing_fractions[cycle] << std::endl;

        std::ofstream out("result.txt", std::ios::app);
        out << std::fixed << std::setprecision(8);
        for (const auto& t : particle_types) out << t.label << "\tprob=" << t.prob << "\t";
        out << std::endl;
        out << "counts\t";
        for (size_t i = 0; i < particle_types.size(); ++i) out << counts[cycle][i] << "\t";
        out << "packing_fraction\t" << packing_fractions[cycle] << std::endl;
        out.close();

        delete[] grid;
        delete[] vox;
    }

    RsaSummary summary;
    summary.avg_counts.assign(particle_types.size(), 0.0);
    for (size_t i = 0; i < particle_types.size(); ++i) {
        double sum = 0.0;
        for (int c = 0; c < cycle_num; ++c) sum += counts[c][i];
        summary.avg_counts[i] = sum / static_cast<double>(cycle_num);
    }
    double phi_sum = 0.0;
    for (double phi : packing_fractions) phi_sum += phi;
    summary.avg_packing_fraction = phi_sum / static_cast<double>(cycle_num);
    return summary;
}

#endif
