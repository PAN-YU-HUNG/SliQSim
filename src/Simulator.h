#ifndef _SIMULATOR_H_
#define _SIMULATOR_H_

#include <iostream>
#include <stdio.h> // FILE
#include <unordered_map>
#include <sys/time.h> //estimate time
#include <fstream> //fstream
#include <sstream> // int to string
#include <cstdlib> //atoi
#include <string> //string
#include <sstream>
#include <random>
#include <cmath>
#include <vector>
#include "../cudd/cudd/cudd.h"
#include "../cudd/cudd/cuddInt.h"
#include "../cudd/util/util.h"


#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899

class Simulator
{
public:
    // constructor and destructor
    Simulator(int type, int nshots, int seed, int bitSize, bool reorder, bool alloc) :
    n(0), r(bitSize), w(4), k(0), inc(3), shift(0), error(0),
    normalize_factor(1), rus_normalize_factor(1), gatecount(0), NodeCount(0), isMeasure(0), shots(nshots), isReorder(reorder), isAlloc(alloc)
    , sim_type(type), statevector("null"), gen(std::default_random_engine(seed)){
    }
    Simulator(int nshots, int seed, int bitSize, bool reorder, bool alloc) :
    n(0), r(bitSize), w(4), k(0), inc(3), shift(0), error(0),
    normalize_factor(1), rus_normalize_factor(1), gatecount(0), NodeCount(0), isMeasure(0), shots(nshots), isReorder(reorder), isAlloc(alloc)
    , sim_type(0), statevector("null"), gen(std::default_random_engine(seed)){
    }
    ~Simulator()  {
        clear();
    }

    /* gates */
    void Toffoli(int targ, std::vector<int> cont, std::vector<int> ncont);
    void Fredkin(int swapA , int swapB, std::vector<int> cont);
    void Peres(int a, int b, int c);
    void Peres_i(int a, int b, int c);
    void Hadamard(int iqubit);
    void rx_pi_2(int iqubit);
    void ry_pi_2(int iqubit);
    void Phase_shift(int phase, int iqubit); // phase can only be 2 to the power of an integer
    void Phase_shift_dagger(int phase, int iqubit);
    void PauliX(int iqubit);
    void PauliY(int iqubit);
    void PauliZ(std::vector<int> iqubit); // Z or CZ
    void measure(int qreg, int creg);
    void RUS(std::vector<int> mqubits, std::vector<int> cond);
    void rz_RUS(int iqubit, double angle);

    /* measurement */
    void getExpectVal();
    void measurement();
    void getStatevector();

    /* simulation */
    void init_simulator(int n);
    void sim_qasm_file(std::string qasm);
    void sim_qasm(std::string qasm);
    void print_results();

    /* misc */
    void reorder();
    void decode_entries();
    void print_info(double runtime, size_t memPeak);
    void setUpRUS(int _rus);    // using RUS
    void check_and_build_rus(const std::string& epsilon, const std::string& delta);

private:
    DdManager *manager;
    DdNode ***All_Bdd;
    int n; // # of qubits
    int r; // resolution of integers
    int w; // # of integers
    int k; // k in algebraic representation
    int inc; // add inc BDDs when overflow occurs, used in alloc_BDD
    int shift; // # of right shifts
    int shots;
    int sim_type; // 0: statevector, 1: measure
    bool isMeasure;
    bool isReorder;
    bool isAlloc;
    int nClbits;
    std::vector<int> expval_qubits;
    double expval;
    std::vector<std::vector<int>> measured_qubits_to_clbits; // empty if not measured
    std::string measure_outcome;
    double normalize_factor; // normalization factor used in measurement
    double rus_normalize_factor; // normalization factor used in RUS
    DdNode *bigBDD; // big BDD used if measurement
    std::default_random_engine gen; // random generator
    std::unordered_map<DdNode *, double> Node_Table; // key: node, value: summed prob
    std::unordered_map<std::string, int> state_count;
    std::string statevector;
    std::string run_output; // output string for Qiskit

    unsigned long gatecount;
    unsigned long NodeCount;
    double error;

    /* measurement */
    double get_total_prob(DdNode *node, int kd2, int nVar, int nAnci_fourInt);
    double measure_probability(DdNode *node, int kd2, int nVar, int nAnci_fourInt, int edge);
    void measure_one(int position, int kd2, double H_factor, int nVar, int nAnci_fourInt, std::string *outcome);
    void measure_and_collapse(std::unordered_map<int,int>& measured_qubits_to_clbits);
    void build_bigBDD(int nAnci_oneInt, int nAnci_fourInt);
    void collapse_to(std::unordered_map<int, int>& qubit_to_state, bool reset_to_zero=true);

    /* misc */
    void init_state(int *constants);
    void init_state_by_matrix(int state_k, std::vector<std::vector<int>>& state);
    void alloc_BDD(DdNode ***Bdd, bool extend);
    void dropLSB(DdNode ***Bdd);
    int overflow3(DdNode *g, DdNode *h, DdNode *crin);
    int overflow2(DdNode *g, DdNode *crin);
    void nodecount();
    std::string find_rz_RUS_gate(const double& angle);

    /* Using RUS */
    bool usingRUS;
    int rus_anci;// rus available when rus_anci >= 0
    const std::string rus_gates_path="./rus/rus_precompile";
    const std::string rus_syn_executable="./rus/rusSyn";

    // Clean up Simulator
    void clear() {
        for (int i = 0; i < w; i++)
            for (int j = 0; j < r; j++)
                Cudd_RecursiveDeref(manager, All_Bdd[i][j]);
        for (int i = 0; i < w; i++)
            delete[] All_Bdd[i];
        delete [] All_Bdd;
        measured_qubits_to_clbits.clear();
        measure_outcome.clear();
        Node_Table.clear();
        state_count.clear();
        Cudd_Quit(manager);
    };
};

#endif
