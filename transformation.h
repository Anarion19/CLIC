//
// Created by kacper on 22.11.18.
//

#ifndef INC_2D_NEW_TRANSFORMATION_H
#define INC_2D_NEW_TRANSFORMATION_H
#include <string>
#include <map>
#include <glob.h>
#include "data_holders.h"

using namespace std;

class ordered_files
{
private:
    map< double, single_file > data;

public:
    const single_file & get_value(double key) const;
    ordered_files(const string &);
    unsigned long get_size() const { return data.size(); }

    using const_iterator = map< double , single_file >::const_iterator;

    const_iterator begin() const { return data.begin(); }
    const_iterator end() const { return data.end(); }
};

class transformed_file: public single_file
{
public:
    transformed_file(const ordered_files &, const single_file &, int interpolation = 1);
};

class transformat
{
public:
    transformat(const map_files &, const ordered_files &, int interpolation = 1);
    void save_to_dir(const string &);
private:
    vector<single_file> data;
};

transformat::transformat(const map_files & map, const ordered_files & ordered_files1, int interpolation)
{
    data.reserve(map.get_size());
    for(auto it: map)
    {
        data.emplace_back(transformed_file(ordered_files1, it, interpolation));
    }
}

void transformat::save_to_dir(const string & folder)
{
    map_files map(data);
    map.save_to_dir(folder);
}

transformed_file::transformed_file(const ordered_files & files, const single_file & source, int interpolation)
{
    this -> name = source.getName();
    data.reserve(files.get_size());
    auto i = files.begin();

    for(auto it : source.getData())
    {
        if(it.first == i -> first)
        {
            for(int k = 0; k < interpolation; k++)
            {
                pair<double, double> tmp;
                //cout << it.first << " " << i -> first;
                tmp.first = i->first;
                tmp.second = source.multiply(i->second);
                //cout << tmp.first<< " " << tmp.second << endl;
                data.emplace_back(tmp);
                i++;
            }
        }
    }

    //show_yourself();
}

const single_file & ordered_files::get_value(double key) const
{
    auto it = data.find(key);
    if(it != data.end())
    {
        return it -> second;
    } else
    {
        throw invalid_argument("no such key " + to_string(key));
    }
}

ordered_files::ordered_files(const string & folder)
{
    glob_t glob_result;
    glob(folder.data(), GLOB_TILDE, NULL, & glob_result);
    string tmp;

    for(int i = 0; i < glob_result.gl_pathc - 1; i++)
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

        //cout << stuff[stuff.size() - 1] << " " << glob_result.gl_pathv[i] << endl;
        data[stuff[stuff.size() - 1]] = single_file(glob_result.gl_pathv[i]);
    }
}
#endif //INC_2D_NEW_TRANSFORMATION_H
