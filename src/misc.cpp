#include "Simulator.h"
#include "util_sim.h"
#if __cplusplus >= 201703L
    #include <filesystem>
    namespace filesystem = std::filesystem;
#else
    #include <experimental/filesystem>
    namespace filesystem = std::experimental::filesystem;
#endif


/**Function*************************************************************

  Synopsis    [initialize state vector by a basis state]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::init_state(int *constants)
{
    DdNode *var, *tmp;
    All_Bdd = new DdNode **[w];
    for (int i = 0; i < w; i++)
        All_Bdd[i] = new DdNode *[r];

    for (int i = 0; i < r; i++)
    {
        if (i == 0)
        {
            for (int j = 0; j < w - 1; j++)
            {
                All_Bdd[j][i] = Cudd_Not(Cudd_ReadOne(manager));
                Cudd_Ref(All_Bdd[j][i]);
            }
            All_Bdd[w - 1][i] = Cudd_ReadOne(manager);
            Cudd_Ref(All_Bdd[w - 1][i]);
            for (int j = n - 1; j >= 0; j--)
            {
                var = Cudd_bddIthVar(manager, j);
                if (constants[j] == 0)
                    tmp = Cudd_bddAnd(manager, Cudd_Not(var), All_Bdd[w - 1][i]);
                else
                    tmp = Cudd_bddAnd(manager, var, All_Bdd[w - 1][i]);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(manager, All_Bdd[w - 1][i]);
                All_Bdd[w - 1][i] = tmp;
            }
        }
        else
        {
            for (int j = 0; j < w; j++)
            {
                All_Bdd[j][i] = Cudd_Not(Cudd_ReadOne(manager));
                Cudd_Ref(All_Bdd[j][i]);
            }
        }
    }
}

/**Function*************************************************************
  Synopsis    [initialize state vector by a initial state matrix]
  Description [state matrix : w x 2^n integer matrix.
               E.g. state[3][0] means the d component of coefficient of the first basis state]
  SideEffects []
  SeeAlso     []
***********************************************************************/
void Simulator::init_state_by_matrix(int state_k, std::vector<std::vector<int>>& state)
{
    this->k = state_k;

    // DdNode *var, *tmp;
    All_Bdd = new DdNode **[w];
    for (int i = 0; i < w; i++)
        All_Bdd[i] = new DdNode *[r];

    // For each a,b,c,d
    for (int i = 0; i < w; i++)
    {
        // Initialize r BDDs to constant 0
        for (int b=0; b<r; b++)
        {
            All_Bdd[i][b] = Cudd_Not(Cudd_ReadOne(manager));
            Cudd_Ref(All_Bdd[i][b]);
        }
        // For each basis state
        for (int j = 0;j < pow(2,n);j++)
        {
            int state_int = state[i][j];
            // For each bit in the integer
            for (int b=0; b<r; b++)
            {
                // Add minterm if needed
                if ((state_int >> b)%2 == 0)
                    continue;
                else
                {
                    // Convert j to its bit string
                    std::string state_str(n,'0');
                    for (int x = n-1; x >= 0; x--)
                    {
                        if ((j >> x)%2 == 0)
                            continue;
                        else
                        {
                            state_str[x] = '1';
                        }
                    }
                    // Create the minterm and OR it to the BDD All_Bdd[i][b]
                    DdNode *var, *tmp, *tmp_term;
                    tmp_term = Cudd_ReadOne(manager);
                    Cudd_Ref(tmp_term);
                    for (int x = n-1; x >= 0; x--)
                    {
                        var = Cudd_bddIthVar(manager, x);
                        if (state_str[x] == '0')
                            tmp = Cudd_bddAnd(manager, Cudd_Not(var), tmp_term);
                        else
                            tmp = Cudd_bddAnd(manager, var, tmp_term);
                        Cudd_Ref(tmp);
                        Cudd_RecursiveDeref(manager, tmp_term);
                        tmp_term = tmp;
                    }
                    tmp = Cudd_bddOr(manager,All_Bdd[i][b], tmp_term);
                    Cudd_Ref(tmp);
                    Cudd_RecursiveDeref(manager, All_Bdd[i][b]);
                    All_Bdd[i][b] = tmp;
                }
            }
        }
    }
}

/**Function*************************************************************

  Synopsis    [allocate new BDDs for each integer vector]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::alloc_BDD(DdNode ***Bdd, bool extend)
{
    DdNode *tmp;

    DdNode ***W = new DdNode **[w];
    for (int i = 0; i < w; i++)
        W[i] = new DdNode *[r];

    for (int i = 0; i < r - inc; i++)
        for (int j = 0; j < w; j++)
            W[j][i] = Bdd[j][i];

    for (int i = 0; i < w; i++)
        delete[] Bdd[i];

    for (int i = 0; i < w; i++)
        Bdd[i] = W[i];

    if (extend)
    {
        for (int i = r - inc; i < r; i++)
        {
            for (int j = 0; j < w; j++)
            {
                Bdd[j][i] = Cudd_ReadOne(manager);
                Cudd_Ref(Bdd[j][i]);
                tmp = Cudd_bddAnd(manager, Bdd[j][r - inc - 1], Bdd[j][i]);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(manager, Bdd[j][i]);
                Bdd[j][i] = tmp;
            }
        }
    }
}

/**Function*************************************************************

  Synopsis    [Drop LSB and shift right by 1 bit]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::dropLSB(DdNode ***Bdd)
{
    DdNode *tmp;

    for (int i = 0; i < w; i++)
    {
        Cudd_RecursiveDeref(manager, Bdd[i][0]); // drop LSB
        // right shift
        for (int j = 1; j < r; j++)
        {
            Bdd[i][j - 1] = Bdd[i][j];
        }
        // sign extension
        Bdd[i][r - 1] = Cudd_ReadOne(manager);
        Cudd_Ref(Bdd[i][r - 1]);
        tmp = Cudd_bddAnd(manager, Bdd[i][r - 2], Bdd[i][r - 1]);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, Bdd[i][r - 1]);
        Bdd[i][r - 1] = tmp;
    }
}

/**Function*************************************************************

  Synopsis    [detect overflow in integer vectors]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int Simulator::overflow3(DdNode *g, DdNode *h, DdNode *crin)
{
    DdNode *tmp, *dd1, *dd2;
    int overflow;

    dd1 = Cudd_bddXor(manager, g, crin);
    Cudd_Ref(dd1);

    dd2 = Cudd_bddXnor(manager, g, h);
    Cudd_Ref(dd2);

    tmp = Cudd_bddAnd(manager, dd1, dd2);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(manager, dd1);
    Cudd_RecursiveDeref(manager, dd2);

    if (Cudd_CountPathsToNonZero(tmp))
        overflow = 1;
    else
        overflow = 0;
    Cudd_RecursiveDeref(manager, tmp);

    return overflow;
}

/**Function*************************************************************

  Synopsis    [detect overflow in integer vectors -- for the case that h is 0]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int Simulator::overflow2(DdNode *g, DdNode *crin){
    DdNode *tmp;
    int overflow;

    tmp = Cudd_bddAnd(manager, Cudd_Not(g), crin);
    Cudd_Ref(tmp);

    if (Cudd_CountPathsToNonZero(tmp))
        overflow = 1;
    else
        overflow = 0;
    Cudd_RecursiveDeref(manager, tmp);

    return overflow;
}

/**Function*************************************************************

  Synopsis    [decode and print each entry of the state vector]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::decode_entries()
{
    double oneroot2 = 1 / sqrt(2);
    double H_factor = pow(oneroot2, k);
    double re = 0, im = 0;
    int *assign = new int[n];
    unsigned long long nEntries = pow(2, n);
    int oneEntry;
    long long int_value = 0;
    DdNode *tmp;

    for (int i = 0; i < n; i++) //initialize assignment
        assign[i] = 0;

    std::cout << "Amplitudes of the Computational Basis States:" << std::endl;

    for (unsigned long long i = 0; i < nEntries; i++) // compute every entry
    {
        re = 0;
        im = 0;
        for (int j = 0; j < w; j++) // compute every complex value
        {
            int_value = 0;
            for (int h = 0; h < r; h++) // compute every integer
            {
                tmp = Cudd_Eval(manager, All_Bdd[j][h], assign);
                Cudd_Ref(tmp);
                oneEntry = !(Cudd_IsComplement(tmp));
                Cudd_RecursiveDeref(manager, tmp);
                if (h == r - 1)
                    int_value -= oneEntry * pow(2, h);
                else
                    int_value += oneEntry * pow(2, h);
            }
            /* translate to re and im */
            re += int_value * cos((double) (w - j - 1)/w * PI);
            im += int_value * sin((double) (w - j - 1)/w * PI);
        }
        re *= H_factor;
        im *= H_factor;
        std::cout << i << ": " << re << " + " << im << "i" << std::endl;
        full_adder_plus_1(n, assign);
    }

    delete[] assign;
}

/**Function*************************************************************

  Synopsis    [reorder BDDs]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::reorder()
{
    int reorder_signal = Cudd_ReduceHeap(manager, CUDD_REORDER_SYMM_SIFT, 0);
    if (!reorder_signal)
        std::cout << "reorder fails" << std::endl;
}

/**Function*************************************************************

  Synopsis    [update max #nodes]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::nodecount()
{
    unsigned long NodeCount_new = Cudd_ReadNodeCount(manager);
    if (NodeCount_new > NodeCount)
         NodeCount = NodeCount_new;
}

/**Function*************************************************************

  Synopsis    [print statistics]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::print_info(double runtime, size_t memPeak)
{
    std::cout << "  Runtime: " << runtime << " seconds" << std::endl;
    std::cout << "  Peak memory usage: " << memPeak << " bytes" << std::endl; //unit in bytes
    std::cout << "  #Applied gates: " << gatecount << std::endl;
    std::cout << "  Max #nodes: " << NodeCount << std::endl;
    std::cout << "  Integer bit size: " << r << std::endl;
    std::cout << "  Accuracy loss: " << error << std::endl;
    // std::cout << "  #Integers: " << w << std::endl;

    // std::unordered_map<std::string, int>::iterator it;
    // std::cout << "  Measurement: " << std::endl;
    // for(it = state_count.begin(); it != state_count.end(); it++)
    //     std::cout << "      " << it->first << ": " << it->second << std::endl;
}

// using RUS
void Simulator::setUpRUS(int _rus)
{
    usingRUS = _rus >= 0;
    rus_anci = _rus;
}

/**Function*************************************************************
  Synopsis    [find the path of Rz gate file]
  Description []
  SideEffects []
  SeeAlso     []
***********************************************************************/
std::string Simulator::find_rz_RUS_gate(const double& angle){

    std::ifstream config_file(rus_gates_path + "/config.yaml");
    std::string angle_unit;
    if (config_file.is_open()){
        angle_unit = find_in_yaml(config_file, "angle_unit");
    }
    else{
        throw std::runtime_error("Error: cannot open " + rus_gates_path + "/config.yaml");
    }
    if(angle_unit.empty()){
        throw std::runtime_error("Error: angle_unit not found in " + rus_gates_path + "/config.yaml");
    }
    double angle_unit_d = parse_theta(angle_unit);
    if(angle_unit_d <= 0){
        throw std::runtime_error("Error: wrong format of angle_unit in " + rus_gates_path + "/config.yaml");
    }
    double normalized_angle = normalize_theta(angle);
    int angle_int = static_cast<int>(round(normalized_angle / angle_unit_d));

    config_file.clear();
    config_file.seekg(0);
    std::string file_path = find_in_yaml(config_file, std::to_string(angle_int));
    if(file_path.empty()){
        throw std::runtime_error("Error: angle " + std::to_string(angle) + " not found in " + rus_gates_path + "/config.yaml");
    }
    return rus_gates_path + "/" + file_path;
}


/**Function*************************************************************
  Synopsis    [check and build RUS]
  Description []
  SideEffects []
  SeeAlso     []
***********************************************************************/

void Simulator::check_and_build_rus(const std::string& epsilon, const std::string& delta){
    if (filesystem::exists(rus_gates_path))
    {
        // rus_gates_path exists
        bool fail = false;
        std::ifstream config_file(rus_gates_path + "/config.yaml");
        if (config_file.is_open()){
            std::stringstream buffer;
            buffer << config_file.rdbuf();
            std::string epsilon_str = find_in_yaml(buffer, "epsilon");
            buffer.clear();
            buffer.seekg(0);
            std::string delta_str = find_in_yaml(buffer, "angle_unit");

            if(!(epsilon_str.empty() || delta_str.empty())){

                double epsilon_config = parse_theta(epsilon_str);
                double delta_config = parse_theta(delta_str);
                double epsilon_d = parse_theta(epsilon);
                double delta_d = parse_theta(delta);
                if(abs(epsilon_config - epsilon_d) > 1e-10 || abs(delta_config - delta_d) > 1e-10) {
                    std::cout
                        << "Warning: epsilon or delta in " + rus_gates_path + "/config.yaml"
                        << " is different from the input epsilon or delta. \n"
                        << "Rebuilding RUS gates..."
                        << std::endl;
                    fail = true;
                }
            } else {
                std::cerr
                    << "Error: epsilon or delta not found in " + rus_gates_path + "/config.yaml\n"
                    << "Rebuilding RUS gates..."
                    << std::endl;
                fail = true;
            }
        }
        else{
            std::cerr
                << "Error: cannot open " + rus_gates_path + "/config.yaml\n"
                << "Rebuilding RUS gates..."
                << std::endl
            ;
            fail = true;
        }

        if(!fail)return;
        else filesystem::remove_all(rus_gates_path);
    } else {
        // rus_gates_path does not exist
        std::cout << "Built RUS gate not found. \nBuilding RUS gates..."<<std::endl;
        filesystem::create_directories(rus_gates_path);
    }
    std::string command = rus_syn_executable
        + " -O " + rus_gates_path
        + " -E " + epsilon
        + " -P " + delta
        + " -F 4"
        + " -D -1"
        + " -C g-count"
        + " --qubit-name q[__compuation_bit__]"
        + " --ancil-name q[__ancilla_bit__]"
        + " --database"
        ;
    int result = system(command.c_str());
    if (result != 0) {
        throw std::runtime_error("Error: failed to execute rus_syn_executable");
    }
}
