//
// Created by kacper on 01.11.18.
//

#ifndef INC_2D_NEW_TRANSFER_H
#define INC_2D_NEW_TRANSFER_H

#include <iostream>
#include <vector>
#include <cmath>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <Eigen/Dense>
using namespace std;

typedef Eigen::Matrix<boost::multiprecision::cpp_bin_float_100 ,Eigen::Dynamic,Eigen::Dynamic>  MatrixXmp;
typedef Eigen::Matrix<boost::multiprecision::cpp_bin_float_100 ,Eigen::Dynamic,1>               VectorXmp;

template<typename T>
std::string tostring(const T &n) {
    std::ostringstream oss;
    oss << n;
    string s =  oss.str();
    int dotpos = s.find_first_of('.');
    if(dotpos!=std::string::npos){
        int ipos = s.size()-1;
        while(s[ipos]=='0' && ipos>dotpos){
            --ipos;
        }
        s.erase ( ipos + 1, std::string::npos );
    }
    return s;
}

bool out_range(double in)
{
    return in < 335 || in >=  357;
}

double gaussian(double x, double sigma, double mean = 0)
{
    return (1 / sqrt(2.0 * M_PI * sigma * sigma)) * exp(- ( (x - mean) * (x-mean) ) / (2.0 *sigma*sigma) );
}

vector<double> parse_file_name(const string & file_name)
{
    vector<string> parameter;
    boost::split(parameter, file_name, boost::is_any_of("/"));

    string list = parameter[parameter.size() -1];
    parameter.clear();

    boost::split(parameter, list, boost::is_any_of("_"));

    vector<double> coordinats;
    coordinats.reserve(parameter.size());

    parameter.erase(parameter.begin());

    for(auto & p : parameter)
    {
        coordinats.emplace_back(stod(p));
    }

    return coordinats;
}

class Parameter
{
public:
    using const_iterator = vector< pair<bool, double> >::const_iterator;

    Parameter(const vector<pair<bool, double>> &pairs, const vector<double>&);
    Parameter() = default;

    const bool compare (const vector <double> &) const;
    const vector<pair<bool, double>> &getPairs() const { return pairs; }
    int get_fixed() const { return fixed; }
    const unsigned long getSize() const { return pairs.size(); }
    const unsigned long get_free() const { return pairs.size() - fixed; }

    const_iterator begin() const { return pairs.begin(); }
    const_iterator end() const { return pairs.end(); }

private:
    vector< pair<bool, double> > pairs;
    vector< double > tolerance;
    int fixed = 0;
};

Parameter::Parameter(const vector<pair<bool, double>> &pairs, const vector<double>& t) : pairs(pairs), tolerance(t)
{
    for(auto & i : pairs)
    {
        if(i.first)
        {
            fixed++;
        }
    }
}

const bool Parameter::compare(const vector<double> & in) const
{
    if(in.size() != pairs.size()) { return false; }

    for(int i = 0; i < pairs.size(); i++)
    {
        if(pairs[i].first)
        {
            if( pairs[i].second != in[i])
            {
                return false;
            }
        } else if(abs(pairs[i].second - in[i]) >= tolerance[i])
        {
            return false;
        }
    }

    return true;
}

#endif //INC_2D_NEW_TRANSFER_H
