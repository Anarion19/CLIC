//
// Created by kacper on 05.11.18.
//

#ifndef INC_2D_NEW_SIMULATION_H
#define INC_2D_NEW_SIMULATION_H

#include "new_case.h"
#include "transfer.h"
#include "data_holders.h"

#include <iostream>
#include <algorithm>
#include <omp.h>
#include <Eigen/Dense>
#include <Eigen/LU>

#define RES 100.0
using namespace std;

class Simulation
{
public:
    void magic(const vector<int>&);
    void calc_avg();
    void map_sigmas(const vector<double> & sigma_alfa, const vector<double>& sigma_n, const string &, const vector<int>&, double);
    void los_yukawos(const vector<double> & sigma_alfa, const vector<double>& sigma_y, const string &, const vector<int>&, double);
    void las_backgroudnas(const vector<double> & sigma_alfa, const vector<double>& sigma_b, const string &, const vector<int>&, double);
    void lumos(const vector<double> & lumi, const string & file_name, const vector<int> & free);
    void histograms(const string&, const string&);
    void improvment(const string &, int);
    void print() const;
    void print_avg();
    void to_file(const string &, bool) const;
    void coolness(double, function<bool (Case&, Case&)>, string);
    void setL(double l) { for(auto & it : all) { it.setL(l); } }
    void calc_sigma() const;

    vector<double> get_avg_results();
    vector<double> get_avg_sigma();
    vector<double> get_gaussian_results(double, double);
    vector<double> get_gaussian_sigma(double, double);

    Simulation(unsigned,const map_files &, const Case &, const string &, unsigned, const string &);
    Simulation(const Case &, const map_files &, int);
    Simulation(const vector<double> &, const string &, const map_files &, int);
    Simulation(const vector<double> &, const string &, const map_files &, int, double);
    Simulation(const vector<double> &, const string &, const map_files &, double, double, int);
    Simulation(const vector<Case::point>&, const string &, const map_files &, int, double);
    Simulation(const single_file &, const map_files &, int, double);

    virtual ~Simulation() = default;

protected:
    void add(Case & x) { all.emplace_back(x); }

    vector<Case> all;
    map_files files;
    Case reference;
    MatrixXmp basic;

    void calc_base();
};

void Simulation::magic(const vector<int>& free)
{
    reference.score(basic, files, free);
#pragma omp parallel
    {
#pragma omp for schedule(auto)
        for(int i = 0; i < all.size(); i++)
        {
            all[i].score(basic, files, free);
            //all[i].to_file_fit("out5.txt");
        }
    }
}

void Simulation::calc_base()
{
    vector<double> coord;
    int free = files.getSet_up().get_free();
    //base = new CMatrix("base", 2*free + (free - 1)*free*0.5 + 1, 2*free + (free - 1)*free*0.5 + 1);
    basic = MatrixXmp::Zero(2*free + (free - 1)*free*0.5 + 1, 2*free + (free - 1)*free*0.5 + 1);
    double normalization = files.sum_of_dist();

    for(auto & it : files)
    {
        coord.clear();

        double punishment = 1 /( 2 * it.getDistance() + 1);

        auto & tmp = files.getSet_up().getPairs();
        //jeśli wymiar nie jest ustalony jest dodawany do listy coordynatów
        for(int i = 0; i < it.getCoordinats().size(); i++)
        {
            if(! tmp[i].first)
            {
                coord.emplace_back(it.getCoordinats()[i]);
                //cout << it.getCoordinats()[i] << endl;
            }
        }

        int matrix = 0;
        //CMatrix out("out",2*free + (free - 1)*free*0.5 + 1, 1);
        VectorXmp out_e = VectorXmp::Zero(2*free + (free - 1)*free*0.5 + 1);
        //out.m_pData[matrix][0] += 1;
        out_e(matrix) += 1;
        matrix++;

        for(int i = 0; i < coord.size(); i++)
        {
            for(int j = i + 1; j < coord.size(); j++)
            {
                //out.m_pData[matrix][0] += coord[i] * coord[j];
                out_e(matrix) += coord[i] * coord[j];
                //cout << j * i << endl;
                matrix++;
            }
        }

        for(int i = 0; i < coord.size(); i++)
        {
            //out.m_pData[matrix][0] += coord[i];
            out_e(matrix) += coord[i];
            //cout <<  i << endl;
            matrix++;
        }

        for(int i = 0; i < coord.size(); i++)
        {
            //out.m_pData[matrix][0] += coord[i]*coord[i];
            out_e(matrix) += coord[i] * coord[i];
            //cout << i * i << endl;
            matrix++;
        }

        out_e = sqrt(punishment) * out_e;

        //cout << out << endl;
        //cout << matrix << endl;
        basic += (out_e * out_e.transpose());
        //*base = *base + (out * out.Transpose());
    }

    basic = basic.llt().solve(MatrixXmp::Identity(2*free + (free - 1)*free*0.5 + 1,2*free + (free - 1)*free*0.5 + 1));
    //cout << basic << endl;
    //*base = base->Inverse();
    //cout << *base << endl;
}

Simulation::Simulation(const single_file & source, const map_files & files, int how_many, double sigma) : files(files), reference(10, source, 20)
{
    all.reserve(how_many);

    for(int i = 0; i < how_many; i++)
    {
        all.emplace_back(Case(10, source, sigma));
    }

    calc_base();
}

Simulation::Simulation(unsigned case_size, const map_files & files, const Case & reference, const string & source, unsigned size, const string & distr = ""): files(files), reference(reference)
{
    all.reserve(size);

    for(int i = 0; i < size; i++)
    {
        all.emplace_back( distr.empty() ? Case(case_size, source) : Case(case_size, source, distr));
    }

    calc_base();
}

Simulation::Simulation(const Case & aCase, const map_files & files, int how_many = 1): files(files), reference(aCase)
{
    all.reserve(how_many);

    for(int i = 0; i < how_many; i++)
    {
        all.push_back(aCase);
    }

    calc_base();
}

Simulation::Simulation(const vector<double> & config, const string & source_file, const map_files & files, int how_many = 1): files(files), reference(config, source_file)
{
    all.reserve(how_many);
    for(int i = 0; i < how_many; i++)
    {
        all.emplace_back(Case(config, source_file));
    }

    calc_base();
}

void Simulation::map_sigmas(const vector<double> &sigma_alfa, const vector<double> &sigma_n, const string & file, const vector<int>& free, double sigma = 0)
{
    fstream out;
    out.open(file, ios_base::out);
    for (auto it : sigma_alfa)
    {
            for (auto jt : sigma_n)
            {
#pragma omp parallel
                {
#pragma omp for schedule(auto)
                    for (int i = 0; i < all.size(); i++)
                    {
                        all[i].SIGMA_A = it;
                        all[i].SIGMA_N = jt;
                        all[i].min_chi = 0;
                        all[i].score(basic, files, free);
                    }
                }

                out << it << " " << jt << " ";

                vector<double> tmp;

                tmp = (sigma == 0) ? get_avg_sigma() : get_gaussian_sigma(sigma, 340);

                for(auto kt : tmp)
                {
                    out << kt << " ";
                }

                tmp = (sigma == 0) ? get_avg_results() : get_gaussian_results(sigma, 340);

                for(auto kt : tmp)
                {
                    out << kt << " ";
                }
                out << endl;
            }
            out << endl;
    }
    out.close();
}

void Simulation::histograms(const string & file_results, const string & file_map)
{
    fstream file;
    file.open(file_results, ios_base::out);

    map<double, int> dane;
    double size = 0.01;

    for(auto &it : all)
    {
        double key = floor(it.results[0] / size) * size;

        if(dane.find(key) != dane.end())
        {
            dane[key] += 1;
        } else
        {
            dane[key] = 1;
        }
    }

    for(auto &it : dane)
    {
        file << it.first << " " << it.second << endl;
    }

    file.close();
}

void Simulation::improvment(const string & file, int which)
{
    fstream out;
    out.open(file, ios_base::out);
    map<double, int> improvment;
    for(auto it : all)
    {

        int add = 1;
        if(it.sigma.at(which) < reference.sigma.at(which))
        {
            for(auto jt : it.data)
            {
                improvment[jt.first] += add;
            }
        }


    }

    for(auto it : improvment)
    {
        out << it.first << " " << it.second << endl;
    }
    out.close();
}

void Simulation::to_file(const string & file, bool full = false) const
{
    fstream out;
    out.open(file, ios_base::out);
    for(auto it : all)
    {
        it.print(out, full);
        out << endl;
    }
    out.close();
}

void Simulation::print() const
{
    for(auto it : all)
    {
        it.print(cout, false) << endl;
    }
}

void Simulation::calc_avg()
{
    vector<double> avg = all.begin() -> sigma;
    for(int i = 1; i < all.size(); i++)
    {
        for(int j = 0; j < all[i].sigma.size(); j++)
        {
            avg[j] += all[i].sigma[j];
        }
    }

    for(auto & it : avg)
    {
        it = it/all.size();
        cout << it << endl;
    }
}

Simulation::Simulation(const vector<double> & config, const string & source, const map_files & files, int how_many, double sigma): files(files), reference(config, source)
{
    all.reserve(how_many);
    default_random_engine engine(chrono::system_clock::now().time_since_epoch().count());
    normal_distribution<double> gauss(0.0, sigma);

    for(int i = 0; i < how_many; i++)
    {
        vector<double> tmp(config);
        double move = round(gauss(engine)) / RES;
        //cout << move << endl;
        for(auto & it : tmp)
        {
            it += move;
            assert(it >= 335 && it < 357);
        }

        all.emplace_back(Case(tmp, source));
    }

    calc_base();
}

Simulation::Simulation(const vector<double> & config, const string & source, const map_files & files, double sigma, double width, int how_many): files(files), reference(config, source)
{
    all.reserve(how_many * 10000);

    for(int i = 0; i < how_many; i++)
    {
        vector<double> tmp(config);
        double move = round(-sigma*width/2 + ((sigma*width) / how_many) * i ) / RES;

        cout << move << endl;
        for(auto & it : tmp)
        {
            it += move;
            assert(it >= 335 && it < 357);
        }

        for(int j =0; j < 10000; j++)
        {
            all.emplace_back(Case(tmp, source));
        }
    }

    vector<double> tmp(config);
    double move = round( sigma ) / RES;

    cout << move << endl;
    for(auto & it : tmp)
    {
        it += move;
        assert(it >= 335 && it < 357);
    }

    for(int j =0; j < 10000; j++)
    {
        all.emplace_back(Case(tmp, source));
    }

    calc_base();
}

Simulation::Simulation(const vector<Case::point> & config, const string & source, const map_files & files, int how_many, double sigma): files(files), reference(config, source)
{
    all.reserve(how_many);
    default_random_engine engine(chrono::system_clock::now().time_since_epoch().count());
    normal_distribution<double> gauss(0.0, sigma);

    for(int i = 0; i < how_many; i++)
    {
        vector<Case::point > tmp(config);
        double move = round(gauss(engine)) / RES;
        //cout << move << endl;
        for(auto & it : tmp)
        {
            it.first += move;
            assert(it.first >= 335 && it.first < 357);
        }

        all.emplace_back(Case(tmp, source));
    }

    calc_base();
}

void Simulation::los_yukawos(const vector<double> &sigma_alfa, const vector<double> &sigma_y, const string & file, const vector<int> & free, double sigma = 0)
{
    fstream out;
    out.open(file, ios_base::out);
    for (auto it : sigma_alfa)
    {
        for (auto jt : sigma_y)
        {
#pragma omp parallel
            {
#pragma omp for schedule(auto)
                for (int i = 0; i < all.size(); i++)
                {
                    all[i].SIGMA_A = it;
                    all[i].SIGMA_Y = jt;
                    all[i].min_chi = 0;
                    all[i].score(basic, files, free);
                }
            }

            out << it << " " << jt << " ";

            vector<double> tmp;

            tmp = (sigma == 0) ? get_avg_sigma() : get_gaussian_sigma(sigma, 340);

            for(auto kt : tmp)
            {
                out << kt << " ";
            }

            tmp = (sigma == 0) ? get_avg_results() : get_gaussian_results(sigma, 340);

            for(auto kt : tmp)
            {
                out << kt << " ";
            }

            out << endl;
        }
        out << endl;
    }
    out.close();
}

vector<double> Simulation::get_avg_results()
{
    auto tmp = all[0].results;

    for(int j = 1; j < all.size(); j++)
    {
        for(int i = 0; i < tmp.size(); i++)
        {
            tmp[i] += all[j].results[i];
        }
    }

    for(auto & it : tmp)
    {
        it /= all.size();
    }


    return tmp;
}

vector<double> Simulation::get_avg_sigma()
{
    vector<double> tmp(all[0].sigma.size(), 0.0);

    for(auto & it : all)
    {
        for(int i = 0; i < tmp.size(); i++)
        {
            tmp[i] += it.sigma[i];
        }
    }

    for(auto & it : tmp)
    {
        it /= all.size();
    }

    return tmp;
}

vector<double> Simulation::get_gaussian_results(double sigma, double zero)
{
    vector<double> tmp(all[0].results.size(), 0.0);
    double weight = 0;

    for(auto & it : all)
    {
        double x = it.data[0].first - zero;
        weight += gaussian(x, sigma);

        for(int i = 0; i < tmp.size(); i++)
        {
            tmp[i] += it.results[i] * gaussian(x, sigma);
        }
    }

    for(auto & it : tmp)
    {
        it /= weight;
    }


    return tmp;
}

vector<double> Simulation::get_gaussian_sigma(double sigma, double zero)
{
    vector<double> tmp(all[0].sigma.size(), 0.0);
    double weight = 0;

    for(auto & it : all)
    {
        double x = it.data[0].first - zero;
        double w = gaussian(x, sigma);
        weight += w;

        for(int i = 0; i < tmp.size(); i++)
        {
            tmp[i] += it.sigma[i] * w;
        }
    }

    for(auto & it : tmp)
    {
        it /= weight;
    }


    return tmp;
}

void Simulation::coolness(double treshold, function<bool (Case&, Case&)> compare, string file_name)
{
    sort(all.begin(), all.end(), compare);

    map<double, int> mapa;

    for(int i = 0; i < all.size()*treshold; i++)
    {
        for(auto& it : all[i].data)
        {
            double x = it.first - all[i].move/RES;
            if(mapa.count(x) > 0)
            {
                mapa[x]++;
            } else
            {
                mapa[x] = 1;
            }
        }
    }

    fstream file(file_name, ios_base::out);
    for(auto & it : mapa)
    {
        file << it.first << " " << it.second << endl;
    }
    file.close();
}

void Simulation::calc_sigma() const {

}

void Simulation::las_backgroudnas(const vector<double> &sigma_alfa, const vector<double> &sigma_b, const string & file, const vector<int> & free, double sigma = 0)
{
    fstream out;
    out.open(file, ios_base::out);
    for (auto it : sigma_alfa)
    {
        for (auto jt : sigma_b)
        {
#pragma omp parallel
            {
#pragma omp for schedule(auto)
                for (int i = 0; i < all.size(); i++)
                {
                    all[i].SIGMA_A = it;
                    all[i].SIGMA_B = jt;
                    all[i].min_chi = 0;
                    all[i].score(basic, files, free);
                }
            }

            out << it << " " << jt << " ";

            vector<double> tmp;

            tmp = (sigma == 0) ? get_avg_sigma() : get_gaussian_sigma(sigma, 340);

            for(auto kt : tmp)
            {
                out << kt << " ";
            }

            tmp = (sigma == 0) ? get_avg_results() : get_gaussian_results(sigma, 340);

            for(auto kt : tmp)
            {
                out << kt << " ";
            }
            out << endl;
        }
        out << endl;
    }
    out.close();
}

void Simulation::lumos(const vector<double> &lumi, const string &file_name, const vector<int> & free)
{
    fstream out;
    out.open(file_name, ios_base::out);
    for (auto l : lumi)
    {
#pragma omp parallel
        {
#pragma omp for schedule(auto)
            for (int i = 0; i < all.size(); i++)
            {
                all[i].setL(l);
                all[i].min_chi = 0;
                all[i].score(basic, files, free);
            }
        }

        out << l << " ";

        vector<double> tmp = get_avg_sigma();

        for(auto kt : tmp)
        {
            out << kt << " ";
        }

        tmp = get_avg_results();

        for(auto kt : tmp)
        {
            out << kt << " ";
        }
        out << endl;
    }
    out.close();
}

void Simulation::print_avg()
{
    auto tmp = get_avg_results();
    for(auto it : tmp)
    {
        cout << it << " ";
    }
    cout << endl;
    tmp = get_avg_sigma();
    for(auto it : tmp)
    {
        cout << it << " ";
    }
    cout << endl;
}

#endif //INC_2D_NEW_SIMULATION_H
