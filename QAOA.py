from pennylane import qaoa
from pennylane import numpy as np
from matplotlib import pyplot as plt
import networkx as nx
import pennylane as qml
from pennylane.tape import QuantumScript
from qiskit import QuantumCircuit
from mqt.ddsim.primitives.estimator import Estimator
from mqt.ddsim import CircuitSimulator
import os
import sys
import copy
import ast
import multiprocessing

if (len(sys.argv) != 4):
    print("Wrong number of inputs !!!")
    sys.exit()

debug = False
# Debug mode

case_index = int(sys.argv[3])
# Index of the case

if case_index == 0:
    n_vertex = 6
    edges = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 4), (3, 4), (3, 5)]
elif case_index == 1:
    n_vertex = 6
    edges = [(0, 2), (1, 2), (2, 3), (2, 4), (2, 5)]
elif case_index == 2:
    n_vertex = 10
    edges = [(0, 9), (1, 9), (2, 9), (3, 9), (4, 9), (5, 9), (6, 9), (7, 9), (8, 9)]
else:
    print("The case number doesn't exist !!!")
    sys.exit()
# select case and define graph

depth = 1
# repeat number of QAOA layers

n_step = 20
# optimization iterations

op_type = [qaoa.min_vertex_cover,
           qaoa.max_independent_set,
           qaoa.maxcut,
           qaoa.min_vertex_cover,
           ][0]
# select one optimization type
# https://docs.pennylane.ai/en/stable/code/qml_qaoa.html

params = np.array([[0.1]*depth, [0.5]*depth], requires_grad=True)
# initial parameters

if sys.argv[1] == '1':
    use_outer_simulator = True
    outer_simulator_index = int(sys.argv[2])
else:
    use_outer_simulator = False
    dev_index = int(sys.argv[2])
# select simulator with program inputs

outer_simualtor_string = ["SliQSim", "DDSIM" , "SliQSim_RUS"]
# outer_simulator_index = 0
if use_outer_simulator:
    outer_simulator = outer_simualtor_string[outer_simulator_index]
else :
    outer_simulator = ""
# outer simulator used

rz_precision = np.pi/512
rz_data_point = int(2 * np.pi / rz_precision + 0.5)
if outer_simulator == "SliQSim" or debug:
    data_folder = "gridsynth_precompile"
    data_name = "gridsynth_%d.data" % rz_data_point
    assert(data_name in os.listdir(data_folder))
    with open(data_folder + "/" + data_name) as file:
        precompile = file.read().split('\n')[:rz_data_point]
# set rz gate precision (temporary approach) and read needed file for SliQSim

if outer_simulator == "SliQSim_RUS" or debug:
    qasm_path_sliqsim = "temp.qasm"
    os.popen("echo 'qreg q[0];' > " + qasm_path_sliqsim)
    os.popen("./SliQSim --sim_qasm "
    + qasm_path_sliqsim
    + (" --type 1 --rus %d --rus_epsilon 1e-6 --rus_delta %.20f" %(0 , rz_precision))).read().split('\n')

dev_string = ["default.qubit", "qiskit.aer", "lightning.gpu"]
# dev_index = 0
# list and select the device

if not use_outer_simulator:
    dev = qml.device(dev_string[dev_index], wires=n_vertex)
else:
    dev = qml.device(dev_string[0], wires=n_vertex)
# device used for pennylane, can be "qiskit.aer", "qulacs.simulator", etc.

if use_outer_simulator:
    output_fig_name = "./fig/QAOA_" + outer_simulator + "_" + str(case_index) + ".png"
else:
    output_fig_name = "./fig/QAOA_" + dev_string[dev_index].replace('.','_') + "_" + str(case_index) + ".png"
# Output figure name

# ==========================================================================

def loop_task(i, basic_file_content, pstr_ops, pstr_coeffs, q_count, rz_precision):
    file_content = basic_file_content

    # expval
    for op in pstr_ops[i]:
        qubit_wire_z = int(op[1:].strip('()'))
        file_content += ("exp_val q[%d]; \n" % (qubit_wire_z))

    qasm_path_sliqsim = f"sliqsim{i}.qasm"
    with open(qasm_path_sliqsim, 'w') as file:
        file.write(file_content)
    # print("command: ")
    # print("./SliQSim --sim_qasm "
    # + qasm_path_sliqsim
    # + (" --rus %d --rus_epsilon 1e-6 --rus_delta %.20f" %(q_count - 1 , rz_precision)))
    # execute
    exe_result = os.popen("./SliQSim --sim_qasm "
    + qasm_path_sliqsim
    + (" --rus %d --rus_epsilon 1e-6 --rus_delta %.20f" %(q_count - 1 , rz_precision))).read().split('\n')
    for line in exe_result:
        if line.startswith("The expectation value is"):
            exp_v = float(line.split()[-1])
    return pstr_coeffs[i] * exp_v

def qaoa_layer(gamma, alpha):
    qaoa.cost_layer(gamma, cost_h)
    qaoa.mixer_layer(alpha, mixer_h)

@qml.qnode(dev)
def circuit(params, **kwargs):
    for w in range(n_vertex):
        qml.Hadamard(wires=w)
    qml.layer(qaoa_layer, depth, params[0], params[1])
    return qml.expval(cost_h)

def circuit_outer(params):
    # export qasm file
    circuit.construct([params], {})

    # read qasm file and execute the circuit
    # Pauli strings are stored in "pauli_strings"
    # process pauli string
    pstr_coeffs = pauli_strings.coeffs
    pstr_ops = [[s.strip() for s in str(ops).split('@')] for ops in pauli_strings.ops]

    if outer_simulator == "SliQSim":
        qasm_string_original = circuit.qtape.to_openqasm()
        basic_file_content = ""

        # QAOA circuit
        for qasm_str in qasm_string_original.split('\n'):
            # print(qasm_str)
            if 'measure' in qasm_str:
                continue
            if 'rz' not in qasm_str:
                basic_file_content += qasm_str + '\n'
            else :
                qasm_str = qasm_str.split(' ')
                rz_theta = float(qasm_str[0].strip('rz()'))
                qubit_wire = int(qasm_str[1].strip('q[];'))
                multiple = round(float( rz_theta / rz_precision)) % rz_data_point
                for gate in precompile[multiple]:
                    if gate == "Q":
                        gate = "sdg"
                    gate = gate.lower()
                    basic_file_content += ("%s q[%d]; \n" % (gate, qubit_wire))

        result = 0
        for i in range(len(pstr_ops)):
            file_content = basic_file_content

            # expval
            for op in pstr_ops[i]:
                qubit_wire_z = int(op[1:].strip('()'))
                file_content += ("exp_val q[%d]; \n" % (qubit_wire_z))

            qasm_path_sliqsim = "sliqsim.qasm"
            with open(qasm_path_sliqsim, 'w') as file:
                file.write(file_content)

            # execute
            exe_result = os.popen("./SliQSim --sim_qasm " + qasm_path_sliqsim).read().split('\n')
            for line in exe_result:
                if line.startswith("The expectation value is"):
                    exp_v = float(line.split()[-1])
            # print(exp_v)
            result += pstr_coeffs[i] * exp_v

        return np.array(result)
    elif outer_simulator == "DDSIM":
        qasm_string_original = circuit.qtape.to_openqasm()

        # QAOA circuit
        qc = QuantumCircuit(n_vertex)
        for qasm_str in qasm_string_original.split('\n'):
            # print(qasm_str)
            if 'measure' in qasm_str:
                continue
            if 'rz' in qasm_str:
                rz_theta, qubit_wire = qasm_str.split(' ')
                rz_theta = float(rz_theta.strip('rz()'))
                qubit_wire = int(qubit_wire.strip('q[];'))
                qc.rz(rz_theta, qubit_wire)
            else :
                qasm_str = qasm_str.split(' ')
                if (qasm_str[0] == 'h'):
                    qubit_wire = int(qasm_str[1].strip('q[];'))
                    qc.h(qubit_wire)
                elif (qasm_str[0] == 'cx'):
                    qubit_wire = [int(w.strip('q[];')) for w in qasm_str[1].split(',')]
                    qc.cx(qubit_wire[0],qubit_wire[1])


        result = 0
        estimator = Estimator()
        for i in range(len(pstr_ops)):
            # execute
            p_str = ["I"] * n_vertex
            for op in pstr_ops[i]:
                qubit_wire_z = int(op[1:].strip('()'))
                p_str[n_vertex - 1 - qubit_wire_z] = "Z"
            p_str = ''.join(p_str)
            job = estimator.run([qc],[p_str])
            exp_v = job.result().values[0]
            # print(exp_v)
            result += pstr_coeffs[i] * exp_v

        return np.array(result)
    elif outer_simulator == "SliQSim_RUS":
        # TODO : Combine SliQSim and RUS to calculate the expectation value of a Pauli string
        qasm_string_original = circuit.qtape.to_openqasm()
        basic_file_content = ""

        q_count = None
        # QAOA circuit
        for qasm_str in qasm_string_original.split('\n'):
            # print(qasm_str)
            # add ancilla qubit
            if 'qreg' in qasm_str:
                # print(qasm_str)
                before_q = qasm_str.split('q[')[0]
                q_count = int(qasm_str.split('q[')[1].split(']')[0])
                q_count = q_count + 1
                basic_file_content += ("%s q[%d]; \n" % (before_q, q_count))
            elif 'measure' in qasm_str:
                continue
            else :
                basic_file_content += qasm_str + '\n'
        assert(q_count is not None)
        assert(q_count == n_vertex + 1)
        result = 0

        args = [(i, basic_file_content, pstr_ops, pstr_coeffs, q_count, rz_precision) for i in range(len(pstr_ops))]
        with multiprocessing.Pool(8) as pool:
            result = sum(pool.starmap(loop_task, args))
        return np.array(result)

        # for i in range(len(pstr_ops)):
        #     file_content = basic_file_content

        #     # expval
        #     for op in pstr_ops[i]:
        #         qubit_wire_z = int(op[1:].strip('()'))
        #         file_content += ("exp_val q[%d]; \n" % (qubit_wire_z))

        #     qasm_path_sliqsim = "sliqsim.qasm"
        #     with open(qasm_path_sliqsim, 'w') as file:
        #         file.write(file_content)
        #     # print("command: ")
        #     # print("./SliQSim --sim_qasm "
        #     # + qasm_path_sliqsim
        #     # + (" --rus %d --rus_epsilon 1e-6 --rus_delta %.20f" %(q_count - 1 , rz_precision)))
        #     # execute
        #     exe_result = os.popen("./SliQSim --sim_qasm "
        #     + qasm_path_sliqsim
        #     + (" --rus %d --rus_epsilon 1e-6 --rus_delta %.20f" %(q_count - 1 , rz_precision))).read().split('\n') # todo: user adjustable epsilon and delta
        #     for line in exe_result:
        #         if line.startswith("The expectation value is"):
        #             exp_v = float(line.split()[-1])
        #             # Weird bug when pstr_ops[i] has 2 elements (temporary fix)
        #             if len(pstr_ops[i]) == 2:
        #                 exp_v = -exp_v
        #     # print(exp_v)
        #     result += pstr_coeffs[i] * exp_v

        # return np.array(result)
        # return np.array(0) # Output the expectation value here
    else:
        return np.array(0)

@qml.qnode(dev)
def prob_circuit(params):
    for w in range(n_vertex):
        qml.Hadamard(wires=w)
    qml.layer(qaoa_layer, depth, params[0], params[1])
    return qml.probs(wires=range(n_vertex))

def prob_circuit_outer(params, shots=1000):
    # export qasm file
    circuit.construct([params], {})

    # read qasm file and execute the circuit
    # Pauli strings are stored in "pauli_strings"
    if outer_simulator == "SliQSim":
        qasm_string_original = circuit.qtape.to_openqasm()
        basic_file_content = ""

        # QAOA circuit
        for qasm_str in qasm_string_original.split('\n'):
            # print(qasm_str)
            if 'measure' in qasm_str:
                basic_file_content += qasm_str + '\n'
            elif 'rz' not in qasm_str:
                basic_file_content += qasm_str + '\n'
            else :
                qasm_str = qasm_str.split(' ')
                rz_theta = float(qasm_str[0].strip('rz()'))
                qubit_wire = int(qasm_str[1].strip('q[];'))
                multiple = round(float( rz_theta / rz_precision)) % rz_data_point
                for gate in precompile[multiple]:
                    if gate == "Q":
                        gate = "sdg"
                    gate = gate.lower()
                    basic_file_content += ("%s q[%d]; \n" % (gate, qubit_wire))

        qasm_path_sliqsim = "sliqsim.qasm"
        with open(qasm_path_sliqsim, 'w') as file:
            file.write(basic_file_content)

        # execute
        exe_result = os.popen("./SliQSim --sim_qasm " + qasm_path_sliqsim + (" --shots %d" %(shots))).read().split('\n')
        count = ast.literal_eval(exe_result[0])['counts']
        return [count.setdefault(bin(i).lstrip('0b')[::-1].ljust(n_vertex,'0'),0)/shots for i in range(2**n_vertex)]
    elif outer_simulator == "DDSIM":
        qasm_sim_path = 'ddsim.qasm'
        with open(qasm_sim_path, 'w') as file:
            file.write(circuit.qtape.to_openqasm())

        sim = CircuitSimulator(qasm_sim_path)
        count = sim.simulate(shots = shots)
        return [count.setdefault(bin(i).lstrip('0b').zfill(n_vertex)[::-1],0)/shots for i in range(2**n_vertex)]
    elif outer_simulator == "SliQSim_RUS":
        # TODO: Combine SliQSim and RUS to simulate the QAOA circuit and get the probability distribution
        qasm_string_original = circuit.qtape.to_openqasm()
        basic_file_content = ""

        # QAOA circuit
        q_count = None
        c_count = None
        add_measure = False
        for qasm_str in qasm_string_original.split('\n'):
            # print(qasm_str)
            # add ancilla qubit
            if 'qreg' in qasm_str:
                before_q = qasm_str.split('q[')[0]
                q_count = int(qasm_str.split('q[')[1].split(']')[0])
                q_count = q_count + 1
                basic_file_content += ("%s q[%d]; \n" % (before_q, q_count))
            elif 'creg' in qasm_str:
                before_q = qasm_str.split('c[')[0]
                c_count = int(qasm_str.split('c[')[1].split(']')[0])
                c_count = c_count + 1
                basic_file_content += ("%s c[%d]; \n" % (before_q, c_count))
            else :
                basic_file_content += qasm_str + '\n'
            # elif 'measure' in qasm_str and not add_measure:
            #     add_measure = True
            #     basic_file_content += qasm_str + '\n'
            #     basic_file_content += "measure q[%d] -> c[%d]; \n" % (n_vertex, n_vertex)

        basic_file_content += "measure q[%d] -> c[%d]; \n" % (n_vertex, n_vertex)
        qasm_path_sliqsim = "sliqsim.qasm"
        with open(qasm_path_sliqsim, 'w') as file:
            file.write(basic_file_content)

        assert(q_count is not None)
        assert(q_count == c_count)
        assert(q_count == n_vertex + 1)
        # print("command: ")
        # print("./SliQSim --sim_qasm "
        # + qasm_path_sliqsim
        # + (" --shots %d" %(shots))
        # + (" --rus %d --rus_epsilon 1e-6 --rus_delta %.20f" %(q_count - 1, rz_precision)))
        # execute
        exe_result = os.popen("./SliQSim --sim_qasm "
        + qasm_path_sliqsim
        + (" --shots %d" %(shots))
        + (" --rus %d --rus_epsilon 1e-6 --rus_delta %.20f" %(q_count - 1, rz_precision))).read().split('\n') # todo: user adjustable epsilon and delta
        count = ast.literal_eval(exe_result[0])['counts']
        print("histogram : " + str(count))
        return [count.setdefault('0'+(bin(i).lstrip('0b').zfill(n_vertex)[::-1]),0)/shots for i in range(2**n_vertex)]
        # return [0.1]*(2**n_vertex) # return the distribution here
    else:
        return [0.1]*(2**n_vertex)

# Updating parameter when using outer simulator
def step_outer(circuit_outer, params, delta, stepsize):
    updated_params = copy.deepcopy(params)
    for i in range(len(params)):
        for j in range(len(params[i])):
            params[i][j] += delta[i]
            grad = circuit_outer(params)

            params[i][j] -= 2 * delta[i]
            grad -= circuit_outer(params)

            params[i][j] += delta[i]
            grad /= 2 * delta[i]
            print("Computed gradient : " + str(grad))
            updated_params[i][j] = params[i][j] - grad * stepsize

    # print(updated_params)
    return updated_params

# ==========================================================================

# Define the graph, its corresponding QAOA circuit, and used Pauli string
graph = nx.Graph(edges)
cost_h, mixer_h = op_type(graph, constrained = False)
pauli_strings = qml.expval(cost_h).obs

# Debug mode
if debug:
    params = [[-0.2],[0.4]]
    print("Cost(default qubit) : " + str(circuit(params)))

    outer_simulator = "SliQSim_RUS"
    print("Cost(SliQSim_RUS) : " + str(circuit_outer(params)))
    # print(prob_circuit(params))
    # outer_simulator = "SliQSim"
    # print("Cost(SliQSim) : " + str(circuit_outer(params)))
    # print(prob_circuit_outer(params, 100000))
    outer_simulator = "DDSIM"
    print("Cost(DDSIM) : " + str(circuit_outer(params)))
    # print(prob_circuit_outer(params, 100000))

    sys.exit()

# Train the parameter
optimizer = qml.GradientDescentOptimizer()
for i in range(n_step):
    print("%03d"%i, end = ' ', flush=True)
    if use_outer_simulator:
        # TODO(probably solved): Find a good way to update the parameter
        # Let the result of SliQSim converge to a good solution
        stepsize = 0.01
        if outer_simulator == "SliQSim":
            delta = [stepsize, stepsize] # Can be modified
        elif outer_simulator == "SliQSim_RUS":
            delta = [stepsize, stepsize] # Can be modified
        else :
            delta = [stepsize, stepsize]
        params = step_outer(circuit_outer, params, delta, stepsize)

        # Discretization for SliQSim
        if outer_simulator == "SliQSim":
            for i in range(len(params)):
                for j in range(len(params[i])):
                    params[i][j] = round(float( params[i][j] / rz_precision)) * rz_precision

        print("Cost : " + str(circuit_outer(params)))
    else:
        params = optimizer.step(circuit, params)
        print("Cost : " + str(circuit(params)))
print()

print("Optimal Parameters:")
print(params)
print()

# Compute the probability distribution
if use_outer_simulator:
    shots = 100000
    probs = prob_circuit_outer(params, shots)
else:
    probs = prob_circuit(params)

# Draw the figure
plt.bar(range(2 ** n_vertex), probs)
plt.savefig(output_fig_name)
# plt.show()

sys.exit()