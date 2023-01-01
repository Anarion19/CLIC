//
// Created by kacper on 07.11.2019.
//

#ifndef CLIC_GENETIC_H
#define CLIC_GENETIC_H

#include "simulation.h"

using namespace std;

class Individual : public Case
{
public:
    Individual(const vector<point> & genotype, const single_file & file) : Case(genotype, file) { weight = 10.0 / genotype.size(); }

    [[nodiscard]] auto is_healthy() const -> bool;
    [[nodiscard]] auto is_dominated(const shared_ptr<Individual> & a, const vector<unsigned> & obj) const -> bool;
    [[nodiscard]] auto breed(const shared_ptr<Individual> & spouse, long mutation, const single_file & file) const -> shared_ptr<Individual>;
    void setSigma(const vector<double> & s);

    double rank = 0;
    int front = -1;
    vector<shared_ptr<Individual>> dominators;
};

class NSGA2 : public Simulation
{
public:
    NSGA2(const single_file & source, const map_files & files, unsigned points, unsigned population);
    NSGA2(const single_file & source, const map_files & files, unsigned population);
    NSGA2(const vector<double> v, const single_file & source, const map_files & files, unsigned population);
    //NSGA2(const string &, const map_files &);

    void magic(const vector<int> & free);
    void magic(const vector<shared_ptr<Individual>> & list, const vector<int> &);
    void to_file(const string & file) const;
    void points(const string & file) const;
    void setReference(const Case & ref);
    void move_to_next(int mutation, const single_file &, const vector<int> &, const vector<unsigned> &);

protected:
    vector<shared_ptr<Individual>> generation;

    unsigned number = 0;
};

NSGA2::NSGA2(const single_file & source, const map_files & files, unsigned points, unsigned population)
: Simulation(source, files, 0, 0)
{
    generation.reserve(population);

    long seed =  chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    uniform_int_distribution<long> distr(0, source.get_size() - 1);

    vector<Individual::point> genotype(points);

    for(unsigned i = 0; i < population; i++)
    {
        for(auto  & it : genotype)
        {
            it = Individual::point(source.get_pair(static_cast<unsigned long>(distr(generator))).first, 10.0 / points);
        }

        generation.emplace_back(make_shared<Individual>(genotype, source));
    }

    calc_base();
}

NSGA2::NSGA2(const vector<double> v, const single_file &source, const map_files &files, unsigned population)
: Simulation(source, files, 0, 0)
{
    generation.reserve(population);

    for(unsigned i = 0; i < population; i++)
    {
        vector<Case::point> genotype(v.size());

        for(int j = 0; j < v.size(); j++)
        {
            genotype[j].first = v[j];
            genotype[j].second = 10.0 / v.size();
        }


        generation.emplace_back(make_shared<Individual>(genotype, source));
    }
}

NSGA2::NSGA2(const single_file &source, const map_files &files, unsigned population)
: Simulation(source, files, 0, 0)
{
    generation.reserve(population);

    long seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    uniform_int_distribution<long> distribution(4, 30);
    uniform_int_distribution<long> distr2(0, source.get_size() - 1);

    for(unsigned i = 0; i < population; i++)
    {
        auto size = distribution(generator);
        vector<Case::point> genotype(size);

        for(auto & it : genotype)
        {
            it.first = source.get_pair(static_cast<unsigned long>(distr2(generator))).first;
            it.second = 10.0 / size;
        }

        generation.emplace_back(make_shared<Individual>(genotype, source));
    }

    calc_base();
}

void NSGA2::magic(const vector<int> & free)
{
    reference.score(basic, files, free);

#pragma omp parallel for
    for(int i = 0; i < generation.size(); i++)
    {
        generation[i]->score(basic, files, free);
    }
}

void NSGA2::to_file(const string & file) const
{
    fstream out(file,ios_base::out);

    reference.print(out, false);
    out << endl;
    out << endl;

    for(auto & it : generation)
    {
        it->print(out, false);
        out << endl;
    }

    out.close();
}

void NSGA2::magic(const vector<shared_ptr<Individual>> &list, const vector<int> & free)
{

#pragma omp parallel for
    for(int i = 0; i < list.size(); i++)
    {
        list[i]->score(basic, files, free);
    }
}

void NSGA2::move_to_next(int mutation, const single_file & file, const vector<int> & free, const vector<unsigned> & obj)
{
    vector<shared_ptr<Individual>> children;
    children.reserve(generation.size() * 4);

    long seed =  chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    uniform_int_distribution<long> distribution(0, generation.size() - 1);

//    vector< vector<shared_ptr<Individual>> > lists(obj.size());
//    for(int i = 0; i < obj.size(); i++)
//    {
//        lists[i] = generation;
//        sort(lists[i].begin(), lists[i].end(), [obj, i](const shared_ptr<Individual> & a,
//                const shared_ptr<Individual> & b) -> bool { return a->getSigma()[obj[i]] < b->getSigma()[obj[i]]; });
//    }
//
//    for(int i = 0; i < generation.size(); i++)
//    {
//        for(int j = 0; j < obj.size(); j++)
//        {
//            for(int k = j + 1; k < obj.size(); k++)
//            {
//                children.emplace_back(lists[j][i]->breed(lists[k][i], mutation, file));
//                children.emplace_back(lists[j][i]->breed(lists[k][i], mutation, file));
//            }
//        }
//    }
    for(auto & it : generation)
    {
        long partner = distribution(generator);
        children.emplace_back(it->breed(generation[partner], mutation, file));

        partner = distribution(generator);
        children.emplace_back(it->breed(generation[partner], mutation, file));

        partner = distribution(generator);
        children.emplace_back(it->breed(generation[partner], mutation, file));
    }

    vector<vector<double>> sup(children.size());
    for(auto & it : sup)
    {
        it = vector<double>({0,0,0,0,0});
    }

    for(int i = 0; i < 3; i++)
    {
        magic(children, free);

        children.erase(remove_if(children.begin(), children.end(), [](const shared_ptr<Individual> & a)->bool
                                { return !a->is_healthy(); }), children.end());
#pragma omp parallel for
        for(int j = 0; j < children.size(); j++)
        {
            auto & sig = children[j]->getSigma();

            for(int k = 0; k < sig.size(); k++)
            {
                sup[j][k] = sig[k] > sup[j][k] ? sig[k] : sup[j][k];
            }
        }
    }

    for(int j = 0; j < children.size(); j++)
    {
        children[j]->setSigma(sup[j]);
    }

    children.insert(children.end(), generation.begin(), generation.end());

    sort(children.begin(), children.end(), [obj](const shared_ptr<Individual> & a,
            const shared_ptr<Individual> & b) -> bool { return a->getSigma()[obj[0]] < b->getSigma()[obj[0]]; });

    if(obj.size() > 1)
    {
        for(int i = 0; i < children.size(); i++)
        {
            for(int j = i + 1; j < children.size(); j++)
            {
                if(children[j]->is_dominated(children[i], obj))
                {
                    //cout << i << " dominates " << j << endl;
                    children[i]->dominators.emplace_back(children[j]);
                    children[j]->rank++;
                } else if(children[i]->is_dominated(children[j], obj))
                {
                    //cout << j << " dominates " << i << endl;
                    children[j]->dominators.emplace_back(children[i]);
                    children[i]->rank++;
                }
            }
        }

        vector<shared_ptr<Individual>> fronts;
        fronts.reserve(100);
        int front = 1;
        int how_many = children.size();

        while (how_many > 0 || front > children.size())
        {
            for(auto & p : children)
            {
                if(p->front == -1 && p->rank == 0)
                {

                    fronts.emplace_back(p);
                    how_many--;
                    p->front = front;
                }
            }

            for(auto & p : fronts)
            {
                for(auto & q : p->dominators)
                {
                    q->rank--;
                }
            }

            front++;
//        cout << fronts.size() << " ";
            fronts.clear();
        }

        sort(children.begin(), children.end(), [](const shared_ptr<Individual> & a,
                const shared_ptr<Individual> & b) -> bool { return a->front > b->front; });
    }

    //generation.erase(generation.begin() + generation.size() * 0.9, generation.end());

#pragma omp parallel for
    for(int i = 0; i < generation.size(); i++)
    {
        generation[i] = children[i];
        generation[i]->rank = 0;
        generation[i]->dominators.clear();
        generation[i]->front = -1;
    }

    number++;
}

void NSGA2::points(const string &file) const
{
    map<double, int> mapa;

    for(int i = 0; i < generation.size(); i++)
    {
        for(auto& it : generation[i]->getData())
        {
            double x = it.first - generation[i]->getMove() / RES;
            if(mapa.count(x) > 0)
            {
                mapa[x]++;
            } else
            {
                mapa[x] = 1;
            }
        }
    }

    fstream out(file, ios_base::out);
    for(auto & it : mapa)
    {
        out << it.first << " " << it.second << endl;
    }
    out.close();
}

void NSGA2::setReference(const Case &ref) { reference = ref; }

shared_ptr<Individual>
Individual::breed(const shared_ptr<Individual> &spouse, long mutation, const single_file &file) const
{
    auto & father = this->getData();
    auto & mother = spouse->getData();

    vector<point> genotype;

    long seed =  chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    uniform_int_distribution<long> distribution(-mutation, mutation);
    uniform_int_distribution<unsigned> coin(0, 1);

    auto max = father.size() > mother.size() ? father.size() : mother.size();

    genotype.reserve(max);
    //cout << mother.size() << " " << father.size() << " " << max << endl;
    for(unsigned i = 0; i < max; i++)
    {
        try
        {
            if(distribution(generator) < mutation * 0.9 || genotype.size() <= 2 )
            {
                double dna = coin(generator) ? father.at(i).first : mother.at(i).first;
                double move = static_cast<double>(distribution(generator)) / RES;
                if (out_range(dna + move)) { move *= -1; }

                genotype.emplace_back(dna + move, 0);
            }
        }
        catch (const exception & e) {}
    }

    if(distribution(generator) > mutation * 0.8)
    {
        double dna = genotype[0].first;
        double move = static_cast<double>(distribution(generator)) / RES;
        if (out_range(dna + move)) { move *= -1; }

        genotype.emplace_back(dna + move, 0);
    }


    for(auto & it : genotype)
    {
        it.second = 10.0 / genotype.size();
    }


    return make_shared<Individual>(genotype, file);
}

bool Individual::is_healthy() const
{
    for(auto & it : results)
    {
        if(it < 0.01) return false;
    }

    for(auto & it : sigma)
    {
        if(it < 0) return false;
    }

    if(results.size() - 1 >= YUKAWA_POS)
    {
        if(results[YUKAWA_POS] < 0.5 * YUKAWA_0) return false;
        if(sigma[YUKAWA_POS] < 0.1) return false;
    }

    return true;
}

bool Individual::is_dominated(const shared_ptr<Individual> &a, const vector<unsigned> & obj) const
{
    for(auto & it : obj)
    {
        if(sigma[it] > a->sigma[it]) return false;
    }
    for(auto & it : obj)
    {
        if(sigma[it] < a->sigma[it]) return true;
    }

    return false;
}

void Individual::setSigma(const vector<double> &s) { sigma = s; }


#endif //CLIC_GENETIC_H
