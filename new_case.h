//
// Created by kacper on 10.10.18.
//

#ifndef INC_2D_NEW_NEW_CASE_H
#define INC_2D_NEW_NEW_CASE_H

#include "data_holders.h"
#include "transfer.h"

#include<cstdlib>
#include<string>
#include<fstream>
#include<iostream>
#include<ctime>
#include<random>
#include<algorithm>
#include <chrono>
#include <omp.h>
#include <array>

#define SIG_BG 0.0
#define S_BG 0.0
using namespace std;

class Case
{
public:
    using point = pair<double, double>;

    Case(const unsigned, const single_file &);
    Case(const unsigned, const single_file &, const string&);
    Case(const vector<double> &, const single_file &);
    Case(const vector<point> &, const single_file &);
    Case(const unsigned, const single_file &, double sigma);
    virtual ~Case() = default;

    std::ostream& print(std::ostream&, bool full) const;

    void score(const MatrixXmp & base, const map_files & files, const vector<int>&);
    void to_file_fit(const string &);
    void to_file(const string &);
    void setL(double l) { L = l; EF = l * 702.0; }
    double fit_diff();

    friend std::ostream& operator<< (std::ostream& stream, const Case& aCase);

    const vector<double> &getSigma() const;
    const vector<double> &getResults() const;
    const vector<point> &getData() const;
    int getMove() const;

    double weight = 1;
protected:
    friend class Simulation;
    using coordinates = vector<double>;

    VectorXmp calc_chi(const map_files &);

    double ALPHA_0 = 0.1185;
    double YUKAWA_0 = 1.0;
    double BACKGROUND_0 = 0.073;
    double SIGMA_A = 0.001;
    double SIGMA_Y = 0.1;
    double SIGMA_N = 0.001;
    double SIGMA_B = 0.01;
    double min_chi = 0;
    double L = 10.0;
    double EF = 702.0 * L;
    int move = 0;

    vector< point > data;
    vector<double> poi;
    vector< double> sigma, results, rand_diff;
    vector<boost::multiprecision::cpp_bin_float_100> coeff;
    map< coordinates , double > chi2;
};

void Case::score(const MatrixXmp & base, const map_files &files, const vector<int>& free = {0})
{

    chi2.clear();
    MatrixXmp coef = base * calc_chi(files);
    coeff.clear();
    sigma.clear();
    results.clear();
    coeff.reserve(coef.rows());

    for(int i = 0; i < coef.rows(); i++)
    {
        coeff.emplace_back(coef(i,0));
        //cout << coef(i,0) << endl;
    }
    //cout << endl;

    unsigned long size =  files.getSet_up().get_free();
    unsigned long b = coeff.size() - size;

    coef = MatrixXmp::Zero(size, size);

    int k = 1;
    for(int i = 0; i < size; i++)
    {
        coef(i,i) = coeff.at(i + b); // Współczynniki kwadratowe
        for(int j = i + 1; j < size; j++)
        {
            coef(j,i) = boost::multiprecision::cpp_bin_float_100(0.5) * coeff.at(k); // Współczynniki mieszane
            coef(i,j) = boost::multiprecision::cpp_bin_float_100(0.5) * coeff.at(k);
            k++;
        }
    }

    coef = coef.llt().solve(MatrixXmp::Identity(coef.rows(), coef.cols()));

    sigma.reserve(size);

    for(int i = 0; i < size; i++)
    {
        sigma.emplace_back(static_cast<double>(sqrt(coef(i,i))));
        //cout << sigma[i] << endl;
    }
    //cout <<endl;

    VectorXmp min = VectorXmp::Zero(size);
    coef = MatrixXmp::Zero(size, size);

    k = 1;
    for(int i = 0; i < size; i++)
    {
        coef(i,i) = boost::multiprecision::cpp_bin_float_100(2.0) * coeff.at(i + b);
        for(int j = i + 1; j < size; j++)
        {
            //cout << "i " << i << " j " << j << " k " << k << endl;
            coef(j,i) = coeff.at(k);
            coef(i,j) = coeff.at(k);
            k++;
        }
    }

    //cout << coef << endl;
    coef = coef.llt().solve(MatrixXmp::Identity(coef.rows(), coef.cols()));

    for(int i = 0; i < size; i++)
    {
        min(i) = -coeff.at(coeff.size() - 2 * size + i);
        //cout << min(i) << endl;
    }

    auto out = coef * min;
    results.reserve(size);
//    cout << coeff.at(0) << endl;
//    cout << out << endl;
    min_chi += static_cast<double>(coeff.at(0));

    for(int i = 0; i < size; i++)
    {
//        cout << coeff.at(coeff.size() - size + i) << " * " <<  out(i) << "^2" << endl;
//        cout << coeff.at(coeff.size() - 2 * size + i) << " * " <<  out(i) << endl;
        min_chi += static_cast<double>(coeff.at(coeff.size() - size + i) * out(i) * out(i));
        min_chi += static_cast<double>(coeff.at(coeff.size() - 2 * size + i) * out(i));
    }

    int f = 1;
    for(int j = 0; j < size; j++)
    {
        for(int  k = j + 1; k < size; k++)
        {
//            cout << coeff.at(f) << " * " <<  out(j) << " * " <<  out(k) << endl;
            min_chi += static_cast<double>(coeff.at(f) * out(j) * out(k));
            f++;
        }
    }

//    cout << endl;
//    cout << min_chi << endl;
//    cout << out << endl;
    for(int i = 0; i < size; i++)
    {
        results.emplace_back(static_cast<double>(out(i)) + files.getMin().at(free.at(i)));
    }
    
}

VectorXmp Case::calc_chi(const map_files & files)
{
    vector < double > expected, coord;
//    chi2.clear();
    int free = files.getSet_up().get_free();
    
    if(ALPHA_0 - files.get_min_apha() > 0)
    {
        //cout << "AL_0 " << ALPHA_0 << endl;
        ALPHA_0 -= files.get_min_apha();
        //cout << "AL_0 " << ALPHA_0 << endl;
    }

    if(YUKAWA_0 - files.getMin().at(YUKAWA_POS) > 0)
    {
        YUKAWA_0 -= files.getMin().at(YUKAWA_POS);
    }
    
    VectorXmp out_e = VectorXmp::Zero(2*free + (free - 1)*free*0.5 + 1);

    //double normalization = files.sum_of_dist();

    for(auto & it : files)
    {
        double chi = 0, tmp1 = 0, tmp2 = 0, alp = 0;
        double punishment = 1 /( 2 * it.getDistance() + 1) ;

        //cout << it.getDistance() << endl;
        expected.clear();
        coord.clear();
        //Wyliczanie alfa
        for(int i = 0; i < data.size(); i++)
        {
            expected.push_back(it.get_value(data[i].first));
            //suma wyników
            tmp1 += data[i].second * EF * weight;
            //suma oczekiwanych rezultatów
            tmp2 += expected[i] * EF * weight;
        }
        tmp1 = tmp1 * SIGMA_N * SIGMA_N + 1;
        tmp2 = tmp2 * SIGMA_N * SIGMA_N + 1;
        alp = tmp1 / tmp2;

        //Wyliczanie chi^2

        for(int i = 0; i < data.size(); i++)
        {
            //różnica między wynikiem a oczekiwanym * alfa
            tmp1 = (data[i].second - alp * expected[i])*(data[i].second - alp * expected[i]) * EF * weight * EF * weight;
            //dzielenie przez alfa * oczekiwane
            tmp1 = tmp1 / poi[i];
            //cout << tmp1 << endl;
            chi += tmp1;
        }

        if(SIGMA_N != 0)
        {
            chi += ((alp - 1)*(alp - 1)) / (SIGMA_N * SIGMA_N);
            //cout << "Norm: " << ((alp - 1)*(alp - 1)) / (SIGMA_N * SIGMA_N) << endl;
        }

        if(SIGMA_A != 0)
        {
            //cout << "Alpha: " << (it.getCoordinats().at(1) - ALPHA_0) * (it.getCoordinats().at(1) - ALPHA_0) / (SIGMA_A * SIGMA_A) << endl;
            chi += (it.getCoordinats().at(ALPHA_POS) - ALPHA_0) * (it.getCoordinats().at(ALPHA_POS) - ALPHA_0) / (SIGMA_A * SIGMA_A);
        }

        if(SIGMA_Y != 0)
        {
            //cout << "Yukawa: " << (it.getCoordinats().at(YUKAWA_POS) - YUKAWA_0) * (it.getCoordinats().at(YUKAWA_POS) - YUKAWA_0) / (SIGMA_Y * SIGMA_Y) << endl;
            //cout << "Yukawa_pos: " << YUKAWA_POS << "; Yukawa_0: " << YUKAWA_0 << "; coordinates: " << it.getCoordinats().at(YUKAWA_POS) << "; sigma: " << SIGMA_Y << endl;
            chi += (it.getCoordinats().at(YUKAWA_POS) - YUKAWA_0) * (it.getCoordinats().at(YUKAWA_POS) - YUKAWA_0) / (SIGMA_Y * SIGMA_Y);
        }

        if(SIGMA_B != 0)
        {
            //cout << "Background: " << (it.getCoordinats().at(BCKG_POS) - BACKGROUND_0) * (it.getCoordinats().at(BCKG_POS) - BACKGROUND_0) / (SIGMA_B * SIGMA_B) << endl;
            chi += (it.getCoordinats().at(BCKG_POS) - BACKGROUND_0) * (it.getCoordinats().at(BCKG_POS) - BACKGROUND_0) / (BACKGROUND_0 * BACKGROUND_0 * SIGMA_B * SIGMA_B);
        }

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
        out_e(matrix) += chi * punishment;
        matrix++;
        for(int i = 0; i < coord.size(); i++)
        {
            for(int j = i + 1; j < coord.size(); j++)
            {
                out_e(matrix) += coord[i] * coord[j] * chi * punishment;
                matrix++;
            }
        }
        for(int i = 0; i < coord.size(); i++)
        {
            out_e(matrix) += coord[i] * chi * punishment;
            matrix++;
        }
        for(int i = 0; i < coord.size(); i++)
        {
            out_e(matrix) += coord[i] * coord[i] * chi * punishment;
            matrix++;
        }

        //chi2[coord] = chi;

        //cout << chi << " " << it.getName() << endl;
    }

    return out_e;
}

Case::Case(const unsigned number, const single_file & file, double sigma)
{
    cout << "ERROR: using outdated constructor class Case" << endl;
    long seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    weight = 10.0 / number;
    data.reserve(number);
    uniform_int_distribution<unsigned long> distribution(0, file.get_size() - 1);
    normal_distribution<double> gauss(0, sigma);
    move = gauss(generator);

    for(int  i = 0; i < number; i++)
    {
        auto & tmp = file.get_pair(static_cast<unsigned int>(distribution(generator) + move) % file.get_size());
        poisson_distribution<int> p(tmp.second * EF + (1 + S_BG) * SIG_BG * L);
        data.emplace_back(point(tmp.first, (p(generator) - SIG_BG * L)/ EF));
    }
}

Case::Case(const vector<point> & where, const single_file & file)
{
    long seed =  chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);

    weight = 10.0 / where.size();
    data.reserve(where.size());
    rand_diff.reserve(where.size());
    poi.reserve(where.size());
    for(auto & it : where)
    {
        auto & tmp1 = file.get_pair(it.first);
        poisson_distribution<int> p(tmp1.second * EF * it.second + (1 + S_BG) * SIG_BG * L * it.second);
        auto tp = p(generator);
        data.emplace_back(point(tmp1.first, (tp - SIG_BG * L * it.second)/ (EF * it.second) ));
        poi.emplace_back(tp);
        //rand_diff.emplace_back(tmp1.second - data[i].second);
        //data.emplace_back(point(tmp1.first, tmp1.second));
    }
}

Case::Case(const unsigned number, const single_file & file)
{
    cout << "ERROR: using outdated constructor class Case" << endl;
    long seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    data.reserve(number);
    uniform_int_distribution<unsigned long> distribution(0, file.get_size() - 1);

    for(int i = 0; i < number; i++)
    {
        auto & tmp1 = file.get_pair(distribution(generator));
        poisson_distribution<int> p(tmp1.second * EF + (1 + S_BG) * SIG_BG * L);
        data.emplace_back(point(tmp1.first, (p(generator) - SIG_BG * L)/ EF));
        //data.emplace_back(point(tmp1.first, tmp1.second));
        //cout << tmp1.first << " " << tmp1.second << " " << data[i].second << endl;
    }
}

Case::Case(const unsigned number, const single_file & file, const string & file_with_custom_distribution)
{
    cout << "ERROR: using outdated constructor class Case" << endl;
    fstream infile;
    vector <double> rozklad;
    rozklad.reserve(file.get_size());
    long seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    infile.open(file_with_custom_distribution, ios_base::in);

    while ( !infile.eof() )
    {
        double tmp;
        infile >> tmp;
        infile >> tmp;
        rozklad.push_back(tmp);
    }
    infile.close();

    discrete_distribution <unsigned long> distribution(rozklad.begin(), rozklad.end());

    for(int i = 0; i < number; i++)
    {
        auto & tmp1 = file.get_pair(distribution(generator));
        poisson_distribution<int> p(tmp1.second * EF + (1 + S_BG) * SIG_BG * L);
        data.emplace_back(point(tmp1.first, (p(generator) - SIG_BG * L)/ EF));
        //cout << tmp1.first << " " << tmp1.second << " " << data[i].second << endl;
    }
}

Case::Case(const vector<double> & where, const single_file & file)
{
    weight = 10.0 / where.size();
    long seed =  chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    data.reserve(where.size());
    rand_diff.reserve(where.size());
    poi.reserve(where.size());
    for(int i = 0; i < where.size(); i++)
    {
        auto & tmp1 = file.get_pair(where[i]);
        poisson_distribution<int> p(tmp1.second * EF * weight + (1 + S_BG) * SIG_BG * L * weight);
        auto tp = p(generator);
        data.emplace_back(point(tmp1.first, (tp - SIG_BG * L * weight) / (EF * weight)));
        poi.emplace_back(tp);
        //rand_diff.emplace_back(tmp1.second - data[i].second);
        //data.emplace_back(point(tmp1.first, tmp1.second));
    }
}

const vector<double> &Case::getSigma() const { return sigma; }

const vector<double> &Case::getResults() const { return results; }

const vector<Case::point> &Case::getData() const { return data; }

void Case::to_file_fit(const string & name)
{
    fstream file;
    file.open(name, ios_base::app);
    for(auto & it: chi2)
    {

        for(auto & i : it.first)
        {

            file << i << " ";
        }

        file << it.second << endl;
    }
    file << endl;
    file.close();
}

void Case::to_file(const string & name)
{
    fstream file;
    file.open(name, ios_base::out);
    for(auto & it: coeff)
    {
        //file << it << " ";
    }
    file << endl;
    for(auto & it: data)
    {
        file << it.first << " " << it.second << endl;
    }
    for(auto & it: chi2)
    {
        for(auto & i : it.first)
        {
            file << i << " ";
        }

        file << it.second << endl;
    }
    file.close();
}

std::ostream& Case::print(std::ostream& stream, bool full) const {
    if(full)
    {
        stream << SIGMA_A << " " << SIGMA_N << endl;
        stream << "Chi^2:" << endl;
        for(auto & it: chi2)
        {
            for(auto & i : it.first)
            {
                stream << i << " ";
            }

            stream << it.second << endl;
        }
        stream << "Value: " << endl;
    }

    for(auto &it : results){
        stream << it << " ";
    }
    stream << "| ";

    if(full) stream << endl << "Sigma: " << endl;
    for(auto &it : sigma){
        stream << it << " ";
    }
    stream << "| ";

    if(full) {
        stream << endl << "Coefficients:" << endl;
        for (auto &it : coeff) {
            stream << it << " ";
        }
    }

    if(full) stream << endl << "Data:" << endl;
    stream << data.size() << " ";
    for(auto &it : data)
    {
        stream << it.first << " ";
        if(full) stream << it.second << endl;
    }

    if(full)
    {
        stream << endl << "chi^2:" << endl;

        for(auto &it : chi2)
        {
            for(auto & jt: it.first)
            {
                stream << jt << " ";
            }
             stream << it.second << endl;
        }
        stream << "min_chi^2: " << endl;
        stream << min_chi << endl;
    }

    return stream;
}

std::ostream &operator<<(std::ostream &stream, const Case &aCase)
{
    return aCase.print(stream, true);
}

double Case::fit_diff()
{
    int size = chi2.begin()->first.size();
    double result = 0;
    //vector<char> var = {'x', 'y', 'z', 't', 's'};
    for(auto & it : chi2)
    {
        auto out = it.first;
        boost::multiprecision::cpp_bin_float_100 chi = coeff.at(0);
        for (int i = 0; i < size; i++)
        {
            chi += coeff.at(coeff.size() - size + i) * out.at(i) * out.at(i);
            chi += coeff.at(coeff.size() - 2 * size + i) * out.at(i);
        }

        int f = 1;
        for (int j = 0; j < size; j++)
        {
            for (int k = j + 1; k < size; k++)
            {
                //cout << coeff.at(f) << " * " <<  var.at(j) << " * " <<  var.at(k) << " + ";
                chi += coeff.at(f) * out.at(j) * out.at(k);
                f++;
            }
        }

        //cout << coeff.at(0) << endl;

        result += static_cast<double>(abs(it.second - chi));
    }

    return result;
}

int Case::getMove() const { return move; }

#endif //INC_2D_NEW_NEW_CASE_H
