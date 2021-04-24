#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <cmath>
#include <map>
#include <functional>
#include <cassert>
#include <Windows.h>
#include <limits>
#include <algorithm>
using namespace std;

typedef long double ld;

const ld eps = 1e-10;
const ld M_VALUE = 100000;

struct constraint {
    vector<ld> coeff;
    ld c;

    constraint(vector<ld> coeff, ld c)
    : coeff(move(coeff)), c(c) {
    }

    constraint(const constraint& cnst)
    : coeff(cnst.coeff), c(cnst.c) {
    }

    constraint& operator *=(double m) {
        for (int i = 0; i < (int) coeff.size(); i++) {
            coeff[i] *= m;
        }
        c *= m;
        return *this;
    }

    constraint& operator -= (const constraint &cnst) {
        for (int i = 0; i < (int) coeff.size(); i++) {
            coeff[i] -= cnst.coeff[i];
        }
        c -= cnst.c;
        return *this;
    }

    friend constraint operator * (const constraint &cnst, double m) {
        constraint loc(cnst);
        loc *= m;
        return loc;
    }

    friend ostream& operator << (ostream& out, const constraint& cnst) {
        for (ld val : cnst.coeff) {
            out << val << " ";
        }
        out << "  = " << cnst.c << "\n";
        return out;
    }
};

void modify_constraints(vector<constraint>& constraints, const vector<int>& basis_vars) {
    assert(constraints.size() == basis_vars.size());

    for (int i = 0; i < (int) basis_vars.size(); i++) {
//        cout << "constraints !!! :\n";
//        for (const constraint& c : constraints) {
//            cout << c;
//        }
        int basis = basis_vars[i];
        int max_abs_val = i;
        for (int j = i; j < (int) constraints.size(); j++) {
            if (abs(constraints[j].coeff[basis]) > abs(constraints[max_abs_val].coeff[basis])) {
                max_abs_val = j;
            }
        }
        if (abs(constraints[max_abs_val].coeff[basis]) < eps) {
            throw new runtime_error("Incorrect basis - infinite number of solutions.");
        }
        ld c = constraints[max_abs_val].coeff[basis];
        swap(constraints[i], constraints[max_abs_val]);

        constraints[i] *= 1 / c;

//        cout << "constraints !!! :\n";
//        for (const constraint& c : constraints) {
//            cout << c;
//        }

        for (int j = 0; j < (int) constraints.size(); j++) {
            if (j != i) {
                constraints[j] -= constraints[i] * constraints[j].coeff[basis];
            }
        }
//        cout << "constraints !!! :\n";
//        for (const constraint& c : constraints) {
//            cout << c;
//        }
//        cout << "------\n";
    }

    for (constraint& cnst : constraints) {
        for (ld& val : cnst.coeff) {
            if (abs(val) < eps) {
                val = 0;
            }
        }
        if (abs(cnst.c) < eps) {
            cnst.c = 0;
        }
    }

    for (int i = 0; i < (int) constraints.size(); i++) {
        if (constraints[i].c < 0) {
            cout << "constraints !!! :\n";
            for (const constraint& c : constraints) {
                cout << c;
            }

            throw new runtime_error("Incorrect basis - no solutions.");
        }
    }
}

pair<vector<ld>, ld> simplex_method(const vector<ld>& func_coeff,
                          vector<constraint> constraints,
                          vector<int> basis_vars) {
    if (constraints.size() == 0) {
        for (ld val : func_coeff) {
            if (val < 0) {
                return {{}, -INFINITE};
            }
        }
        return {vector<ld>(func_coeff.size(), 0), 0};
    }
    vector<int> not_basis_vars;
    {
        set<int> set_basis_vars;
        for (int v : basis_vars) {
            set_basis_vars.insert(v);
        }
        for (int i = 0; i < (int) func_coeff.size(); i++) {
            if (set_basis_vars.count(i) == 0) {
                not_basis_vars.push_back(i);
            }
        }
    }

    while (true) {
        modify_constraints(constraints, basis_vars);
        {
            cout << "\n";
            cout << "constraints :\n";
            for (const constraint& c : constraints) {
                cout << "  " << c;
            }

            cout << "basis_vars :\n  ";
            for (int val : basis_vars) {
                cout << val + 1 << " ";
            }
            cout << "\n";

            vector<ld> res(func_coeff.size(), 0);
            ld func_res = 0;
            for (int i = 0; i < basis_vars.size(); i++) {
                res[basis_vars[i]] = constraints[i].c;
                func_res += func_coeff[basis_vars[i]] * constraints[i].c;
            }

            cout << "current result :\n  ";
            for (ld val : res) {
                cout << val << " ";
            }
            cout << "\n";

            cout << "func_res : " << func_res << "\n";
        }

        if (basis_vars.size() == func_coeff.size()) {
            vector<ld> res(basis_vars.size(), 0);
            ld func_val = 0;
            for (int i = 0; i < (int) basis_vars.size(); i++) {
                res[basis_vars[i]] = constraints[i].c;
                func_val = constraints[i].c * func_coeff[basis_vars[i]];
            }
            return {move(res), func_val};
        }

        vector<ld> func_in_not_basis_coeff(func_coeff.size(), 0);
        for (int v : not_basis_vars) {
            func_in_not_basis_coeff[v] += func_coeff[v];
        }
        for (int i = 0; i < (int) basis_vars.size(); i++) {
            for (int v : not_basis_vars) {
                func_in_not_basis_coeff[v] -= func_coeff[basis_vars[i]] * constraints[i].coeff[v];
            }
        }

        cout << "func_in_not_basis_coeff :\n  ";
        for (ld val : func_in_not_basis_coeff) {
            cout << val << " ";
        }
        cout << "\n";

        vector<pair<ld, int>> candidates;
        for (int i = 0; i < not_basis_vars.size(); i++) {
            candidates.emplace_back(func_in_not_basis_coeff[not_basis_vars[i]], i);
        }
        sort(candidates.begin(), candidates.end());

//        int best_v = not_basis_vars.front();
//        int best_v_pos = 0;
//        for (int i = 0; i < (int) not_basis_vars.size(); i++) {
//            if (func_in_not_basis_coeff[not_basis_vars[i]] < func_in_not_basis_coeff[best_v]) {
//                best_v = not_basis_vars[i];
//                best_v_pos = i;
//            }
//        }
        int next_candidate = 0;
        while (next_candidate < candidates.size()) {
            if (candidates[next_candidate].first >= 0) {
                break;
            }
            int not_basis_candidate = not_basis_vars[candidates[next_candidate].second];

            int best_basis_var_pos = -1;
            ld best_time = 0;

            bool is_zero_positive_case = false;

            for (int i = 0; i < (int) basis_vars.size(); i++) {
                if ((constraints[i].coeff[not_basis_candidate] > 0) && (constraints[i].c > 0)) {
                    if ((best_basis_var_pos == -1) || (best_time > constraints[i].c / constraints[i].coeff[not_basis_candidate])) {
                        best_basis_var_pos = i;
                        best_time = constraints[i].c / constraints[i].coeff[not_basis_candidate];
                    }
                }
                if ((constraints[i].coeff[not_basis_candidate] > 0) && (constraints[i].c == 0)) {
                    is_zero_positive_case = true;
                }
            }

            cout << "not_basis_candidate : " << not_basis_candidate + 1 << "    basis_candidate : " << basis_vars[best_basis_var_pos] + 1 << "\n";
//            cout << "best_time : " << best_time << "\n";

            if (is_zero_positive_case) {
                next_candidate++;
            } else {
                if (best_basis_var_pos == -1) {
                    return {{}, -INFINITE};
                }
                swap(basis_vars[best_basis_var_pos], not_basis_vars[candidates[next_candidate].second]);
                break;
            }
        }

        if ((next_candidate == candidates.size()) || (candidates[next_candidate].first >= 0)) {
            vector<ld> res(func_coeff.size(), 0);
            ld func_val = 0;
            for (int i = 0; i < (int) basis_vars.size(); i++) {
                res[basis_vars[i]] = constraints[i].c;
                func_val += constraints[i].c * func_coeff[basis_vars[i]];
            }
            return {move(res), func_val};
        }
    }
    return {{}, numeric_limits<ld>::quiet_NaN()};
}

pair<vector<ld>, ld> m_simplex_method(const vector<ld>& func_coeff,
                          vector<constraint> constraints) {
    vector<ld> loc_func_coeff = func_coeff;
    vector<int> basis_vars(constraints.size(), 0);
    for (int i = 0; i < constraints.size(); i++) {
        basis_vars[i] = func_coeff.size() + i;
        loc_func_coeff.push_back(M_VALUE);
        if (constraints[i].c < 0) {
            constraints[i] *= -1;
        }

        for (int j = 0; j < constraints.size(); j++) {
            if (j == i) {
                constraints[i].coeff.push_back(1);
            } else {
                constraints[i].coeff.push_back(0);
            }
        }
    }
    return simplex_method(loc_func_coeff, move(constraints), basis_vars);
}

void solve(string input, string output) {
    ifstream in(input);
    int n, m, with_basis;
    in >> n >> m >> with_basis;

    vector<ld> func_coeff(n, 0);
    for (int i = 0; i < n; i++) {
        in >> func_coeff[i];
    }

    vector<constraint> constraints;

    for (int i = 0; i < m; i++) {
        vector<ld> coeff(n, 0);
        ld c;
        for (int j = 0; j < n; j++) {
            in >> coeff[j];
        }
        in >> c;
        constraints.emplace_back(move(coeff), c);
    }

    pair<vector<ld>, ld> res;

    if (with_basis != 0) {
        vector<int> basis_vars(m, 0);
        for (int i = 0; i < m; i++) {
            in >> basis_vars[i];
            basis_vars[i]--;
        }
        res = simplex_method(func_coeff, constraints, basis_vars);
    } else {
        res = m_simplex_method(func_coeff, constraints);
    }

    ofstream out(output);
    for (ld val : res.first) {
        out << val << " ";
    }
    out << "\n" << res.second << "\n";
}

int main() {
    #ifdef Vlad_kv
        freopen("input.txt", "r", stdin);
//        freopen("log.txt", "w", stdout);
    #endif // Vlad_kv

    try {
        solve("tests/test_7.txt", "output.txt");
    } catch (runtime_error *err) {
        cout << "runtime_error : " << err->what() << "\n";
    }

    return 0;
}
