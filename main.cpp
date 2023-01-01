#include<vector>
#include<string>
#include<fstream>
#include<iostream>
#include <chrono>
#include <omp.h>
#include <array>
#include "data_holders.h"
#include "new_case.h"
#include "simulation.h"
#include "genetic.h"
#include "transformation.h"

using namespace std;

int main(int argc, char** argv)
{
    string new_file = "./ILC./test/dane_171.5_0.1185_80_350_1.37_1_0.073.txt";
    single_file source(new_file);

    vector <pair <bool, double>> pp = {{0, 171.5}, {0,0.1185}, {1,80.0}, {1,350}, {0,1.37}, {0,1.0}, {0, 0.073}};
    vector <double> tol = {100, 100, 100, 100, 100, 100, 100};
    vector <int> free = {0,1,4,5,6};
    vector <unsigned> objectives = {0,2};
    vector <double> v = {340, 341, 342, 343, 344, 345, 346, 347, 348, 349};
    vector <pair<double, double>> v1 = {{340, 0.5}, {341, 0.5}, {342, 0.5}, {343, 3}, {344, 0.5}, {345, 3.0},
                                        {346, 0.5}, {347, 0.5}, {348, 0.5}, {349, 0.5}};

    vector<double> FCC_mw  = {340.67, 341.2, 341.34, 343.32, 343.06, 344.83, 344.66};
    vector<double> ILC_mw  = {341.39, 341.36, 341.39, 343.26, 343.27, 344.65};
    vector<double> CLIC_mw = {341.23, 341.18, 342.9, 343.01, 344.53};
    vector<double> CLIC_my = {340.45, 342.21, 342.87, 342.72, 345.11, 345.29, 350.68, 348.38, 349.08, 350.56};
    vector<double> CLIC_m  = {339.66, 342.71, 342.76, 342.72, 344.86, 344.77};
    auto start = std::chrono::high_resolution_clock::now();
//    ordered_files files("./spektra/spectra_ilc/*");

    map_files *map_files1 = new map_files("./ILC./test/*", "dane", Parameter(pp, tol), false);
//    map_files1->add_background({0.085, 0.082, 0.079, 0.076, 0.073, 0.070, 0.067, 0.064, 0.061, 0}, "./ILC./test/", "dane");
    //  map_files *map_files1 = new map_files("./dane/*", "dane", Parameter(pp, tol), false);
    //map_files1 -> update_dist(new_file, ILC_mw);
    //map_files1 -> normalize_coords(parse_file_name(new_file));
    Case *aCase = new Case(v, source);
    //Simulation simulation(new_file, *map_files1, 30000, 20);
    //Simulation simulation(10, *map_files1, *aCase, new_file, 10000, "where_4D_2.txt");
//    Simulation simulation(v, new_file, *map_files1, 10000, 0);
    Simulation simulation(ILC_mw, new_file, *map_files1, 40.0, 2, 16);
    //Simulation simulation(1, *map_files1, *aCase, new_file, 50000);

    simulation.magic(free);
    //simulation.coolness(0.1, [](const Case & a, const Case & b) -> bool { return a.getSigma()[0] < b.getSigma()[0]; }, "where.txt");
    //simulation.improvment("improvment.txt", 0);
    //simulation.lumos({10.0, 10.0, 10.0, 10.0, 10.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0}, "CLIC_lum_N0.01A0.001Y0_MS20.txt", free);

//    simulation.histograms("rest2.txt", "hist2d.txt");
    //simulation.clean();
    simulation.to_file("masa_400_5D_ILC_mw.txt");
//    simulation.print_avg();
//    simulation.calc_sigma();

//    simulation.map_sigmas({0.001},
//                          {0.05, 0.04, 0.03, 0.02, 0.01, 0.007, 0.006, 0.005, 0.002, 0.001, 0.0007, 0.0005, 0.0002, 0.0001, 0}, "CLIC_mapa1D_n.txt", free);
//      simulation.las_backgroudnas({0.2, 0.15, 0.1, 0.05, 0.03, 0.01, 0.005, 0.002, 0.001, 0.0001, 0.00001},
//                                  {0.1, 0.09, 0.07, 0.06, 0.05, 0.03, 0.02, 0.01, 0.009, 0.008, 0.007, 0.005, 0.003, 0.002, 0.001}, "CLIC_mapa5D_YB_A0.001.txt", free);
//    simulation.los_yukawos({0.002, 0.001, 0.0007, 0.0005, 0.0002, 0.0001, 0.00009, 0.00007, 0.00006, 0.00005, 0.00001},
//                           {0.2, 0.15, 0.1, 0.05, 0.03, 0.01, 0.005, 0.002, 0.001, 0.0001, 0.00001}, "CLIC_mapa5D_SY_N0.01B0.02Y0.1.txt", free);
//    simulation.lumos({10, 15, 20, 25, 30, 35, 40, 45, 50}, "lumos_CLIC.txt", free);
//    vector<double> lumi = {10, 15, 20, 25, 30, 35, 40, 45, 50};
//    fstream out;
//    out.open("lumos_CLIC.txt", ios_base::out);
//
//    for(auto l : lumi)
//    {
//        Simulation simp(v, new_file, *map_files1, 2000, 0);
//        simp.setL(l);
//        simp.magic(free);
//        out << l << " ";
//
//        vector<double> tmp = simp.get_avg_sigma();
//
//        for(auto kt : tmp)
//        {
//            out << kt << " ";
//        }
//
//        tmp = simp.get_avg_results();
//
//        for(auto kt : tmp)
//        {
//            out << kt << " ";
//        }
//        out << endl;
//    }
//    out.close();

//    transformat transform1(*map_files1, files, 10);
//    transform1.save_to_dir("./ILC");

//    NSGA2 alg(v, source, *map_files1,  1000);
//    alg.setReference(*aCase);
//    alg.magic(free);
//    alg.to_file("przed.txt");
//    double mut = 10;
//    for(int i = 0; i < 30; i++)
//    {
//        alg.move_to_next(mut, source, free, objectives);
//        alg.to_file(to_string(i) + ".txt");
//        mut *= 0.9;
//    }
//    alg.to_file("po.txt");
//    alg.points("distCLIC_mw.txt");

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    cout << elapsed.count() << endl;
    return 0;
}