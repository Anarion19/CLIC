//
// Created by kacper on 14.09.18.
//

#include "transfer.h"

#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iostream>
#include <map>
#include <experimental/filesystem>
#include <glob.h>
#include <utility>
#include <set>
#include <boost/algorithm/string.hpp>

#ifndef INC_2D_NEW_DATA_HOLDERS_H
#define INC_2D_NEW_DATA_HOLDERS_H

#define DEPTH 6
#define ALPHA_POS 1
#define YUKAWA_POS 5
#define BCKG_POS 6
#define EFF 7020.0
#define EPSILON 0.000001

using namespace std;

class single_file
{
protected:
    vector < pair < double, double > > data;
    vector < double > coordinats;
    string name = "";
    double delta;
    double distance = 0;
    void calc_delta();

public:
    void setCoordinats(const vector<double> &coordinats);
    const string &getName() const;
    const pair<double, double> & get_pair(unsigned long) const;
    const pair<double, double> get_pair(double) const;
    const vector<pair<double, double>> &getData() const;
    const vector<double> &getCoordinats() const;

    double multiply(const single_file &) const;
    void calc_dist(const single_file &, const vector<double> & );
    void show_yourself() const;
    void save_to_file(const string &) const;
    void add_background(const string &, double value) const;
    single_file(const string &, const vector <double> &, bool);
    single_file(const string &);
    single_file() = default;

    unsigned long get_size() const;
    const double get_value(double) const;
    const double get_value(unsigned long) const;
    double getDistance() const;
    double get_data(unsigned);

};

class map_files
{
public:
    using iterator = vector<single_file>::iterator;
    using const_iterator = vector<single_file>::const_iterator;

    map_files(const string &, const string &, const Parameter &, bool);
    map_files(const string &);
    map_files(const vector<single_file> & data) : data(data) {}

    const_iterator begin() const { return data.begin(); }
    const_iterator end() const { return data.end(); }
    const_iterator cbegin() const { return data.cbegin(); }
    const_iterator cend() const { return data.cend(); }

    const Parameter & getSet_up() const;
    void save_to_dir(const string & ) const;
    const vector<double> &getMin() const;
    unsigned long get_size() const { return data.size(); }
    double get_min_apha() const { return min.at(ALPHA_POS); }
    void normalize_coords();
    void normalize_coords(const vector<double> & zero);
    void update_dist(const single_file &, const vector<double> &);
    void add_background(const vector<double> &, const string &, const string &) const;
    double sum_of_dist() const;
private:
    vector<single_file> data;
    vector<double> min = {0,0,0,0,0,0,0};
    Parameter set_up;
};

void map_files::normalize_coords()
{
    min = data.begin() -> getCoordinats();
    // Znalezienie mnimalnych
    for(auto & it : data)
    {
        auto tmp = it.getCoordinats();
        for(int i = 0; i < tmp.size(); i++)
        {
            if(min.at(i) > tmp.at(i))
            {
                min[i] = tmp[i];
            }
        }
    }

    // Odejmowanie
    for(auto & it : data)
    {
        auto tmp = it.getCoordinats();
        for(int i = 0; i < tmp.size(); i++)
        {
            tmp[i] -= min[i];

            assert(tmp[i] >= 0);
        }
        it.setCoordinats(tmp);

    }
}

void map_files::normalize_coords(const vector<double> &zero)
{
    min = zero;
    for(auto & it : data)
    {
        assert(zero.size() == it.getCoordinats().size());
        auto tmp = it.getCoordinats();
        for(int i = 0; i < zero.size(); i++)
        {
            tmp.at(i) -= zero[i];
        }
        it.setCoordinats(tmp);
    }
}

void map_files::save_to_dir(const string & dir) const
{
    for(auto it : data)
    {
        it.save_to_file(dir + it.getName());
    }
}

map_files::map_files(const string &folder, const string & key, const Parameter & set_up, bool show_files = false): set_up(set_up)
{
    glob_t glob_result;
    glob(folder.data(), GLOB_TILDE, NULL, & glob_result);
    string tmp;

    for(int i = 0; i < glob_result.gl_pathc; i++)
    {
        vector<string> parameter;
        boost::split(parameter, glob_result.gl_pathv[i], boost::is_any_of("/"));

        tmp = parameter[parameter.size() -1];
        parameter.clear();

        boost::split(parameter, tmp, boost::is_any_of("_"));
        vector<double> stuff;
        stuff.reserve(parameter.size());

        if(parameter[0] == key && parameter.size() > DEPTH)
        {
            parameter.erase(parameter.begin());

            for(auto & p : parameter)
            {
                stuff.push_back(stod(p));
            }

            if(set_up.compare(stuff))
            {
                data.emplace_back(single_file(glob_result.gl_pathv[i], stuff, show_files));
                cout << glob_result.gl_pathv[i] << endl;
            }
        }
        else cerr << "file name incompatible: " <<glob_result.gl_pathv[i] << endl;
    }
}

map_files::map_files(const string & folder)
{
    glob_t glob_result;
    glob(folder.data(), GLOB_TILDE, NULL, & glob_result);
    string tmp;

    for(int i = 0; i < glob_result.gl_pathc; i++)
    {
        vector<string> parameter;
        boost::split(parameter, glob_result.gl_pathv[i], boost::is_any_of("/"));

        tmp = parameter[parameter.size() -1];
        parameter.clear();

        boost::split(parameter, tmp, boost::is_any_of("_"));
        vector<double> stuff;
        stuff.reserve(parameter.size());

        parameter.erase(parameter.begin());

        for(auto & p : parameter)
        {
            stuff.push_back(stod(p));
        }

        data.emplace_back(single_file(glob_result.gl_pathv[i], stuff, false));
    }
}

const Parameter &map_files::getSet_up() const { return set_up; }

const vector<double> &map_files::getMin() const { return min; }

void map_files::update_dist(const single_file & source, const vector<double> & where)
{
    for(auto & it : data)
    {
        it.calc_dist(source, where);
    }
}

double map_files::sum_of_dist() const
{
    double out = 0;

    for(auto & it : data)
    {
        out += 1 /( 2 * it.getDistance() + 1);
    }

    return out;
}

void map_files::add_background(const vector<double> & values, const string & folder, const string & name) const
{
    for(auto & val : values)
    {
        for(auto & file : data)
        {
            string fname = folder + name;
            for(auto & par : file.getCoordinats())
            {
                fname = fname + "_" + tostring(par);
            }
            fname = fname + "_" + tostring(val) + ".txt";
            file.add_background(fname, val);
        }
    }
}

single_file::single_file(const string &file_name, const vector<double> & coordinates, bool show_after_reading = false): name(file_name), coordinats(coordinates)
{
    for(auto it : coordinates)
    {
        cout << it << " ";
    }
    cout << endl;
    fstream input;
    input.open(file_name, ios_base::in);

    if(!input.good() && input.is_open())
    {
        cerr << "Cannot open single_file" << endl;
    }
    else
    {
        double key, value;
        pair<double , double > tmp;
        data.reserve(305);

        while (!input.eof())
        {
            input >> key;
            input >> value;

            tmp.first = key;
            tmp.second = value;
            //cout << key << " " << value << endl;
            data.push_back(tmp);
        }
    }

    data.shrink_to_fit();
    calc_delta();

    if(show_after_reading) show_yourself();
}

single_file::single_file(const string & file_name) : name(file_name)
{
    fstream input;
    input.open(file_name, ios_base::in);

    coordinats = parse_file_name(file_name);

    if(!input.good() || !input.is_open())
    {
        cerr << "Cannot open single_file: " << file_name << endl;
    }
    else
    {
        pair<double , double > tmp;
        data.reserve(305);

        while (!input.eof())
        {
            input >> tmp.first;
            input >> tmp.second;
            //cout << key << " " << value << endl;

            data.push_back(tmp);
        }
    }

    data.shrink_to_fit();
}

double single_file::multiply(const single_file & file) const
{
    if(data.size() != file.data.size())
    {
        cout << this->name << " " << file.name << endl;
        cout << data.size() << " " << file.data.size() << endl;
        throw invalid_argument("different size files cannot be multiplied");
    }

    //cout << data.size() << " " << file.data.size() << endl;

    double out = 0;

    for(int i = 0; i < data.size(); i++)
    {
        if(data[i].first != file.data[i].first)
        {
            //cout <<i << " " << data[i].first << " " << file.data[i].first << endl;
        } else {
            out += data[i].second * file.data[i].second;
        }
    }

    return out;
}

void single_file::show_yourself() const
{
    if ( name.empty() )
    {
        cout << "Nothing inside" << endl;
    } else
    {
        if(data.size() == 0)
        {
            cout << "file name: " << name << "is empty" << endl;
        } else
        {
            cout << "file name: " << name << endl;
            cout << "delta: " << delta <<endl;

            for(auto item : data)
            {
                cout << item.first << " " << item.second << endl;
            }
        }
    }
}

void single_file::save_to_file(const string & file_name) const
{
    //cout << file_name << endl;

    fstream file;
    file.open(file_name, ios_base::out);

    for(auto it : data)
    {
        file << it.first;
        file << " ";
        file << it.second;
        file << endl;
    }

    file.close();
}

void single_file::calc_delta()
{
    auto tmp = data.cbegin();
    delta = tmp -> first;
    tmp++;
    delta = tmp -> first - delta;
}

const pair<double, double> & single_file::get_pair(unsigned long key) const
{
    if(key < data.size())
    {
        return data.at(key);
    }
    else {
        cerr << "single_file::get_pair(unsigned) invalid input: " << key << endl;
        assert(false);
    }
}

const pair<double, double> single_file::get_pair(double key) const
{
//    map < double, double > tmp(data.begin(), data.end());
//    auto it = tmp.find(key);
    auto it = data.begin();

    for(; it != data.end(); it++)
    {
        if(it->first == key) break;
        else if(abs(it-> first - key) < EPSILON)
        {
            //cout << "Epsilon comparison used " << key << " " << it-> first << endl;
            break;
        }
    }

    if(it != data.end())
    {
        return pair<double, double> (it->first, it->second);
    } else
    {
        cerr<< "single_file::get_pair(double) Invalid key: " << key << endl;
        assert(false);
    }
}

double single_file::get_data(unsigned key)
{
    if(key < data.size())
    {
        return data[key].second;
    } else
        {
        cerr << "invalid input" << endl;
    }
}

unsigned long single_file::get_size() const { return data.size(); }

const string &single_file::getName() const { return name; }

const vector<pair<double, double>> &single_file::getData() const { return data; }

const double single_file::get_value(double in) const { return get_pair(in).second; }

const double single_file::get_value(unsigned long in) const { return get_pair(in).second; }

const vector<double> &single_file::getCoordinats() const { return coordinats; }

void single_file::setCoordinats(const vector<double> &coordinats) { single_file::coordinats = coordinats; }

void single_file::calc_dist(const single_file & file, const vector<double> & where)
{
    double out = 0;

    try
    {
        if(data.size() != file.data.size())
        {
            cerr << "Files size doesn't mach" << endl;
            throw "Files size doesn't mach";
        }
        for(auto & it : where)
        {
            auto & pair1 = this->get_pair(it);
            auto & pair2 = file.get_pair(it);

            double tmp = (pair1.second - pair2.second);
            tmp *= tmp;
            tmp /= pair2.second;

            out += tmp * EFF;
        }

        out += where.size();
    }
    catch (exception & e)
    {
        out = 0;
        cerr << "single_file::calc_dist vector where contains illegal position" << endl;
    }

    distance = out;

}

double single_file::getDistance() const { return distance; }

void single_file::add_background(const string & file_name, double value) const
{
    fstream file;
    file.open(file_name, ios_base::out);
    for(auto it : data)
    {
        file << it.first << " " << it.second + value << endl;
    }
    file.close();
}


#endif //INC_2D_NEW_DATA_HOLDERS_H
