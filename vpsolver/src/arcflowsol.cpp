/**
This code is part of the Arc-flow Vector Packing Solver (VPSolver).

Copyright (C) 2013-2016, Filipe Brandao
Faculdade de Ciencias, Universidade do Porto
Porto, Portugal. All rights reserved. E-mail: <fdabrandao@dcc.fc.up.pt>.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <climits>
#include <cstring>
#include <ctime>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include "graph.hpp"
#include "common.hpp"
#include "instance.hpp"
#include "arcflowsol.hpp"
using namespace std;

/* Class ArcflowSol */

ArcflowSol::ArcflowSol(const Instance &_inst, const map<Arc, int> &_flow,
                       int _S, const vector<int> &_Ts, int _LOSS):
        inst(_inst), flow(_flow), S(_S), Ts(_Ts), LOSS(_LOSS) {
    vector<int> dem(inst.m);
    for (int i = 0; i < inst.m; i++) {
        dem[i] = inst.demands[i];
    }

    objvalue = 0;
    sols.resize(inst.nbtypes);
    nbins.resize(inst.nbtypes);
    for (int t = 0; t < inst.nbtypes; t++) {
        sols[t] = extract_solution(&dem, Ts[t]);
        if (!is_valid(sols[t], t)) {
            throw_error("Invalid solution! (capacity)");
        }

        for (const pattern_pair &pat : sols[t]) {
            objvalue += pat.first * inst.Cs[t];
            nbins[t] += pat.first;
        }
    }

    for (int i = 0; i < inst.m; i++) {
        if (dem[i] > 0) {
            throw_error("Invalid solution! (demand)"); // causes error
        }
    }

    int fs = 0;
    for (const auto &kvpair : flow) {
        fs += kvpair.second;
    }
    if (fs != 0) {
        throw_error("Invalid solution! (flow)"); // causes error
    }
}

vector<pattern_pair> ArcflowSol::remove_excess(const vector<pattern_int> &sol,
                                               vector<int> *_dem) const {
    vector<int> &dem = *_dem;
    vector<pattern_pair> tmp;
    for (const pattern_int &pat : sol) {
        map<int, int> count;
        for (const int &it : pat.second) {
            count[it] += 1;
        }
        vector<int> rm;
        int rep = pat.first;
        while (rep > 0) {
            rm.clear();
            for (auto &kvpair : count) {
                int type = inst.items[kvpair.first].type;
                kvpair.second = min(kvpair.second, dem[type]);
                if (kvpair.second == 0) {
                    rm.push_back(kvpair.first);
                }
            }
            for (const int &ind : rm) {
                count.erase(ind);
            }

            int f = rep;
            for (const auto &kvpair : count) {
                int type = inst.items[kvpair.first].type;
                f = min(f, dem[type]/kvpair.second);
            }
            rep -= f;

            tmp.push_back(MP(f, vector<int_pair>(all(count))));
            for (const auto &kvpair : count) {
                int type = inst.items[kvpair.first].type;
                dem[type] -= f * kvpair.second;
            }
        }
    }

    map<vector<int_pair>, int> mp;
    for (pattern_pair &pp : tmp) {
        sort(all(pp.second));
        mp[pp.second] += pp.first;
    }

    vector<pattern_pair> finalsol;
    for (const auto &kvpair : mp) {
        finalsol.push_back(MP(kvpair.second, kvpair.first));
    }
    return finalsol;
}

vector<pattern_pair> ArcflowSol::extract_solution(vector<int> *_dem, int T) {
    vector<int> &dem = *_dem;
    set<int> nodes;
    map<int, vector<Arc>> adj;
    for (const auto &kvpair : flow) {
        int u = kvpair.first.u;
        int v = kvpair.first.v;
        nodes.insert(u);
        nodes.insert(v);
        if (v != S) {
            adj[v].push_back(kvpair.first);
        }
    }

    int &zflow = flow[Arc(T, S, LOSS)];

    vector<int> lst(all(nodes));

    vector<pattern_int> sol;
    //int lars = 0;
    while (true) {
        map<int, Arc> pred;
        map<int, int> dp;
        dp[S] = zflow;
        for (const int &v : lst) {
            int &val = dp[v];
            Arc &p = pred[v];
            for (const Arc &a : adj[v]) {
                
                //if(dp.count(a.u) == 0) {
                    //printf(" %d   dp.count: %ld\n", lars++, dp.count(a.u));
                    //continue;
                //}
                
                throw_assert(dp.count(a.u) != 0); // causes error
                int mf = min(dp[a.u], flow[a]);
                if (mf > val) {
                    p = a;
                    val = mf;
                }
            }
        }
        int f = dp[T];
        zflow -= f;
        if (f == 0) {
            break;
        }
        int v = T;
        sol.push_back(pattern_int());
        pattern_int &patt = sol.back();
        patt.first = f;
        while (v != S) {
            Arc a = pred[v];
            int u = a.u;
            int lbl = a.label;
            if (lbl != LOSS) {
                patt.second.push_back(lbl);
            }
            flow[a] -= f;
            v = u;
        }
    }
    return remove_excess(sol, &dem);
}

bool ArcflowSol::is_valid(const vector<pattern_pair> &sol, int btype) const {
    for (const pattern_pair &pat : sol) {
        vector<int> w(inst.ndims);
        for (const auto &itpair : pat.second) {
            if (inst.binary && itpair.second > 1) {
                return false;
            }
            const Item &it = inst.items[itpair.first];
            for (int i = 0; i < inst.ndims; i++) {
                w[i] += it[i] * itpair.second;
            }
        }
        for (int i = 0; i < inst.ndims; i++) {
            if (w[i] > inst.Ws[btype][i]) {
                return false;
            }
        }
    }
    return true;
}

void ArcflowSol::print_solution(bool print_inst = true, bool pyout = false) {
    std::ofstream destFile("result.txt");
    destFile<<objvalue<<"\n";
   
    for (int t = 0; t < inst.nbtypes; t++) {
        if ((inst.nbtypes > 1) && (nbins[t]!=0)) {
            destFile<<t<<" ";
         //   printf("Machine %d\n", t);
        }
        vector<pattern_pair> &sol = sols[t];
        for (const pattern_pair &pat : sol) {
            vector<int_pair> tmp;
            for (const int_pair &itpair : pat.second) {
                int t = inst.items[itpair.first].type;
                int opt = inst.items[itpair.first].opt;
                for (int i = 0; i < itpair.second; i++) {
                    tmp.push_back(MP(t, opt));
                }
            }
            sort(all(tmp));
          
         //   printf("%d x [", pat.first);
            bool first = true;
            for (const int_pair &p : tmp) {
                if (first) {
                    first = false;
                } else {
                    destFile<<" ";
             //       printf(", ");
                }
                if (p.second == -1) {
                    destFile<<p.first+1;
               //     printf("i=%d", p.first+1);
                } else {
                    printf("i=%d opt=%d", p.first+1, p.second+1);
                }
            }
             destFile<<"\n";
        }
    }
     destFile.close();

    if (print_inst) {
        inst.print();
    }

    if (pyout) {
        printf("PYSOL=(%d,[", objvalue);
        for (int t = 0; t < inst.nbtypes; t++) {
            printf("[");
            vector<pattern_pair> &sol = sols[t];
            for (const pattern_pair &pat : sol) {
                vector<int_pair> tmp;
                for (const int_pair &itpair : pat.second) {
                    int t = inst.items[itpair.first].type;
                    int opt = inst.items[itpair.first].opt;
                    for (int i = 0; i < itpair.second; i++) {
                        tmp.push_back(MP(t, opt));
                    }
                }
                sort(all(tmp));

                printf("(%d,[", pat.first);
                for (const int_pair &p : tmp) {
                    if (p.second == -1) {
                        printf("(%d, 0),", p.first);
                    } else {
                        printf("(%d, %d),", p.first, p.second);
                    }
                }
                printf("]),");
            }
            printf("],");
        }
        printf("])\n");
    }
}
