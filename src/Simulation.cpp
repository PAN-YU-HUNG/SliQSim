#include "Simulator.h"
#include "util_sim.h"


/**Function*************************************************************

  Synopsis    [Initailize simulator]

  Description [This function will set #qubits n, construct initial state, and enable dynamic reordering]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::init_simulator(int nQubits)
{
    n = nQubits; // set the number n here
    manager = Cudd_Init(n, n, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);

    int *constants = new int[n];
    for (int i = 0; i < n; i++)
        constants[i] = 0; // TODO: custom initial state
    measured_qubits_to_clbits = std::vector<std::vector<int>>(n, std::vector<int>(0));

    // Initialize the state
    init_state(constants);

    delete[] constants;
    if (isReorder) Cudd_AutodynEnable(manager, CUDD_REORDER_SYMM_SIFT);
}


/**Function*************************************************************

  Synopsis    [parse and simulate the qasm file]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::sim_qasm_file(std::string qasm)
{
    std::string inStr;
    std::stringstream inFile_ss(qasm);
    while (getline(inFile_ss, inStr))
    {
        inStr = inStr.substr(0, inStr.find("//"));
        if (inStr.find_first_not_of("\t\n ") != std::string::npos)
        {
            std::stringstream inStr_ss(inStr);
            getline(inStr_ss, inStr, ' ');
            if (inStr == "qreg")
            {
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                init_simulator(stoi(inStr));
            }
            else if (inStr == "creg")
            {
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                nClbits = stoi(inStr);
            }
            else if (inStr == "OPENQASM"){;}
            else if (inStr == "include"){;}
            else if (inStr == "measure")
            {
                isMeasure = 1;
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                int qIndex = stoi(inStr);
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                int cIndex = stoi(inStr);
                measure(qIndex, cIndex);
            }
            else if (inStr == "initial_state") // Add initial state
            {
                // Construct the initial state matrix
                const int INITIAL_STATE_PRECISION = 20;
                int state_k = INITIAL_STATE_PRECISION*2;
                std::vector<std::vector<int>> state_matrix(w, std::vector<int>(pow(2,n),0));
                for (int i = 0 ; i < pow(2,n) ; i++)
                {
                    // Read the coefficient of the i-th basis state
                    getline(inFile_ss, inStr);

                    // Get the real and imaginary part
                    std::string doubleStr;
                    std::stringstream inDouble_ss(inStr);
                    double real_part, imag_part;
                    bool use_minus_delimiter = false;
                    bool neg_real = false, neg_imag = false;

                    // Determine whether the real or imaginary part is negative
                    if (inStr[0] == '-')
                    {
                        neg_real = true;
                        inDouble_ss.str(inStr.substr(1));
                    }
                    else if (inStr[0] == '+')
                    {
                        inDouble_ss.str(inStr.substr(1));
                    }
                    if (inDouble_ss.str().find_first_of('-') != std::string::npos)
                    {
                        use_minus_delimiter = true;
                        neg_imag = true;
                    }

                    // Split the string to obtain the real and imaginary part
                    if (use_minus_delimiter)
                        getline(inDouble_ss, doubleStr, '-');
                    else
                        getline(inDouble_ss, doubleStr, '+');
                    real_part = (neg_real)? -stod(doubleStr) : stod(doubleStr);
                    getline(inDouble_ss, doubleStr, 'i');
                    imag_part = (neg_imag)? -stod(doubleStr) : stod(doubleStr);

                    // Add the read coefficient to the initial state matrix
                    // real part to d, imaginary part to b
                    state_matrix[3][i] = int(std::round(real_part*pow(2,INITIAL_STATE_PRECISION)));
                    state_matrix[1][i] = int(std::round(imag_part*pow(2,INITIAL_STATE_PRECISION)));
                }

                // Add initial state
                init_state_by_matrix(state_k,state_matrix);

                // Delete the initial state matrix
                state_matrix.clear();
            }
            else if (inStr == "exp_val") // Calculate expectation value
            {
                getline(inStr_ss, inStr, '[');
                while(getline(inStr_ss, inStr, ']'))
                {
                    expval_qubits.push_back(stoi(inStr));
                    getline(inStr_ss, inStr, '[');
                }
            }
            else if (inStr == "rus")
            {
                if (!usingRUS)
                {
                    std::cerr << std::endl
                              << "[warning]: Syntax \'" << inStr << "\' can only be supported under --rus option. The line is ignored ..." << std::endl;
                    continue;
                }
                std::vector<int> buffer;
                getline(inStr_ss, inStr, '[');
                while (getline(inStr_ss, inStr, ']'))
                {
                    // std::cout<<inStr<<std::endl;
                    buffer.push_back(stoi(inStr));
                    getline(inStr_ss, inStr, '[');
                }
                if(buffer.size() % 2 == 0){
                    std::vector<int> mqubits(buffer.begin(), buffer.begin() + buffer.size() / 2);
                    std::vector<int> cond(buffer.begin() + buffer.size() / 2, buffer.end());
                    assert(mqubits.size() == cond.size());
                    RUS(mqubits, cond);
                } else {
                    assert(buffer.size() == 1);
                    auto pos = inStr.find("==");
                    inStr = inStr.substr(pos+2);
                    std::vector<int> cond(1, stoi(inStr));
                    RUS(buffer, cond);
                }
            }
            else
            {
                if (inStr.find_first_of("rz(") == 0 && inStr.back() == ')')
                {
                    if (!usingRUS)
                    {
                        std::cerr << std::endl
                                  << "[warning]: Syntax \'" << inStr << "\' can only be supported under --rus option. The line is ignored ..." << std::endl;
                        continue;
                    }
                    std::string angle = inStr.substr(3, inStr.size() - 4);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    // std::cout<<"inStr: "<<inStr<<std::endl;
                    // std::cout<<"angle: "<<angle<<std::endl;
                    rz_RUS(stoi(inStr), parse_theta(angle));
                }
                else if (inStr == "x")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    PauliX(stoi(inStr));
                }
                else if (inStr == "y")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    PauliY(stoi(inStr));
                }
                else if (inStr == "z")
                {
                    std::vector<int> iqubit(1);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    iqubit[0] = stoi(inStr);
                    PauliZ(iqubit);
                    iqubit.clear();
                }
                else if (inStr == "h")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    Hadamard(stoi(inStr));
                }
                else if (inStr == "s")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift(2, iqubit);
                }
                else if (inStr == "sdg")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift_dagger(-2, iqubit);
                }
                else if (inStr == "t")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift(4, iqubit);
                }
                else if (inStr == "tdg")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift_dagger(-4, iqubit);
                }
                else if (inStr == "rx(pi/2)")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    rx_pi_2(stoi(inStr));
                }
                else if (inStr == "ry(pi/2)")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    ry_pi_2(stoi(inStr));
                }
                else if (inStr == "cx")
                {
                    std::vector<int> cont(1);
                    std::vector<int> ncont(0);
                    int targ;
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    cont[0] = stoi(inStr);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    targ = stoi(inStr);
                    Toffoli(targ, cont, ncont);
                    cont.clear();
                    ncont.clear();
                }
                else if (inStr == "cz")
                {
                    std::vector<int> iqubit(2);
                    for (int i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        iqubit[i] = stoi(inStr);
                    }
                    PauliZ(iqubit);
                    iqubit.clear();
                }
                else if (inStr == "swap")
                {
                    int swapA, swapB;
                    std::vector<int> cont(0);
                    for (int i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        if (i == 0)
                            swapA = stoi(inStr);
                        else
                            swapB = stoi(inStr);
                    }
                    Fredkin(swapA, swapB, cont);
                    cont.clear();
                }
                else if (inStr == "cswap")
                {
                    int swapA, swapB;
                    std::vector<int> cont(1);
                    for (int i = 0; i < 3; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        if (i == 0)
                            cont[i] = stoi(inStr);
                        else if (i == 1)
                            swapA = stoi(inStr);
                        else
                            swapB = stoi(inStr);
                    }
                    Fredkin(swapA, swapB, cont);
                    cont.clear();
                }
                else if (inStr == "ccx" || inStr == "mcx")
                {
                    std::vector<int> cont(0);
                    std::vector<int> ncont(0);
                    int targ;
                    getline(inStr_ss, inStr, '[');
                    while(getline(inStr_ss, inStr, ']'))
                    {
                        cont.push_back(stoi(inStr));
                        getline(inStr_ss, inStr, '[');
                    }
                    targ = cont.back();
                    cont.pop_back();
                    Toffoli(targ, cont, ncont);
                    cont.clear();
                    ncont.clear();
                }
                else
                {
                    std::cerr << std::endl
                            // << "[warning]: Gate \'" << inStr << "\' is not supported in this simulator. The gate is ignored ..." << std::endl;
                            << "[warning]: Syntax \'" << inStr << "\' is not supported in this simulator. The line is ignored ..." << std::endl;
                }
            }
        }
    }
    if (isReorder) Cudd_AutodynDisable(manager);
}

/**Function*************************************************************

  Synopsis    [simulate the circuit described by a qasm file]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::sim_qasm(std::string qasm)
{
    sim_qasm_file(qasm); // simulate

    if (sim_type == 0 && isMeasure == 0 && expval_qubits.empty())
    {
        std::cout << "Error: no measurement detected. Cannot do sampling.\n" << std::flush;
        assert(sim_type != 0 || isMeasure != 0 || !expval_qubits.empty());
    }
    if (sim_type == 1)
    {
        if (isMeasure == 1)
        {
            std::cout << "Warning: measurement detected. The final statevector will collapse based on the measurement outcome.\n" << std::flush;
            if (shots != 1)
            {
                shots = 1;
                std::cout << "Warning: shot number is limited to 1 in all_amplitude mode.\n" << std::flush;
            }
        }
        else
        {
            if (shots != 1)
            {
                std::cout << "Warning: no measurement detected. The --shots argument is ignored.\n" << std::flush;
            }
        }
    }

    // measure based on simulator type
    if (sim_type == 0) // sampling mode
    {
        measurement();
    }
    else if (sim_type == 1) // all_amplitude mode
    {
        if (isMeasure == 1)
        {
            measurement();
        }
        getStatevector();
    }

    if (!expval_qubits.empty())
    {
        getExpectVal();
    }
    print_results();

    /* Write the DDManager Info
    FILE *pFile;
    pFile = fopen("statistic.txt", "w");
    Cudd_PrintInfo(manager, pFile);
    fclose(pFile); */
}



/**Function*************************************************************

  Synopsis    [print state vector and distribution of sampled outcomes]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::print_results()
{
    // write output string based on state_count and statevector
    std::unordered_map<std::string, int>::iterator it;

    run_output = "{";
    if (state_count.begin() != state_count.end())
    {
        run_output += "\"counts\": { ";
        for (it = state_count.begin(); it != state_count.end(); it++)
        {
            if (std::next(it) == state_count.end())
                run_output = run_output + "\"" + it->first + "\": " + std::to_string(it->second);
            else
                run_output = run_output + "\"" + it->first + "\": " + std::to_string(it->second) + ", ";
        }
        run_output += " }";
        run_output += (statevector != "null") ? ", " : "";
    }

    run_output += (statevector != "null") ? "\"statevector\": " + statevector + " }" : " }";
    std::cout << run_output << std::endl;

    if (!expval_qubits.empty())
    {
        std::cout << "The expectation value is " << expval << std::endl;
    }
}
