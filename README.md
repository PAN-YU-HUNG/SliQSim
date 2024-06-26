# SliQSim - A BDD-based Quantum Circuit Simulator

## Introduction
`SliQSim` is a BDD-based quantum circuit simulator implemented in C/C++ on top of [CUDD](http://web.mit.edu/sage/export/tmp/y/usr/share/doc/polybori/cudd/cuddIntro.html) package. In `SliQSim`, a bit-slicing technique based on BDDs is used to represent quantum state vectors. For more details of the simulator, please refer to the [paper](https://arxiv.org/abs/2007.09304).

## Build
To build the simulator, one needs to first configure `CUDD`:
```commandline
cd cudd
./configure --enable-dddmp --enable-obj --enable-shared --enable-static
cd ..
```
Next, build the binary file `SliQSim`:
```commandline
make
```

## Execution
The circuit format being simulated is `OpenQASM` used by IBM's [Qiskit](https://github.com/Qiskit/qiskit), and the gate set supported in this simulator now contains Pauli-X (x), Pauli-Y (y), Pauli-Z (z), Hadamard (h), Phase and its inverse (s and sdg), π/8 and its inverse (t and tdg), Rotation-X with phase π/2 (rx(pi/2)), Rotation-Y with phase π/2 (ry(pi/2)), Controlled-NOT (cx), Controlled-Z (cz), Toffoli (ccx and mcx), SWAP (swap), and Fredkin (cswap). One can find some example benchmarks in [examples](https://github.com/NTU-ALComLab/SliQSim/tree/master/examples) folder.

For simulation types, we provide both "sampling" and "all_amplitude" simulation options. The help message states the details:

```commandline
$ ./SliQSim --help
Options:
--help                produce help message
--sim_qasm arg        simulate qasm file string
--seed [=arg(=1)]     seed for random number generator
--print_info          print simulation statistics such as runtime, memory, etc.
--type arg (=0)       the simulation type being executed.
                      0: sampling mode (default option), where the sampled outcomes will be provided.
                      1: all_amplitude mode, where the final state vector will be shown.
--shots arg (=1)      the number of outcomes being sampled in "sampling mode" .
--r arg (=32)         integer bit size.
--reorder arg (=1)    allow variable reordering or not.
                      0: disable reordering.
                      1: enable reordering (default option).
--alloc arg (=1)      allocate new BDDs when overflow is detected.
                      0: do not allocate new BDDs. This may lead to numerical errors.
                      1: allocate new BDDs (default option).

```
To use the sampling mode (default), it is required to have measurement operations included in the qasm file. Conversely, in all_amplitude mode, measurement operations are generally omitted, but if they are present in the qasm file, the final state vector will collapse based on the measurement result. It is important to note that all_amplitude mode is not recommended for simulations involving a large number of qubits, as it could result in a significantly long runtime.

For example, simulating [example/bell_state_measure.qasm](https://github.com/NTU-ALComLab/SliQSim/blob/master/examples/bell_state_measure.qasm), which is a 2-qubit bell state circuit with measurement gates at the end, with the sampling mode simulation option can be executed by
```commandline
./SliQSim --sim_qasm examples/bell_state_measure.qasm --type 0 --shots 1024
```

Then the sampled results will be shown:
```commandline
{ "counts": { "11": 542, "00": 482 } }
```

If option `--print_info` is used, simulation statistics such as runtime and memory usage will also be provided:
```commandline
  Runtime: 0.014433 seconds
  Peak memory usage: 12611584 bytes
  #Applied gates: 2
  Max #nodes: 13
  Precision of integers: 32
  Accuracy loss: 2.22045e-16
```

To demostrate the all_amplitude mode simulation, we use [example/bell_state.qasm](https://github.com/NTU-ALComLab/SliQSim/blob/master/examples/bell_state.qasm), which is the same circuit as in the sampling mode simulation example except that the measurement gates are removed:
```commandline
./SliQSim --sim_qasm examples/bell_state.qasm --type 1
```

This will show the resulting state vector:
```commandline
{"statevector": ["0.707107", "0", "0", "0.707107"] }
```

One may also execute our simulator as a backend option of Qiskit through [SliQSim Qiskit Interface](https://github.com/NTU-ALComLab/SliQSim-Qiskit-Interface).


## Citation
Please cite the following paper if you use our simulator for your research:

<summary>
  <a href="https://ieeexplore.ieee.org/document/9586191">Y.-H. Tsai, J.-H. R. Jiang, and C.-S. Jhang, “Bit-slicing the hilbert space:  Scaling up accurate quantum circuit simulation,” in <em>Design Automation Conference (DAC)</em>, 2021, pp. 439–444.</a>
</summary>

```bibtex
@INPROCEEDINGS{9586191,
  author={Tsai, Yuan-Hung and Jiang, Jie-Hong R. and Jhang, Chiao-Shan},
  booktitle={Design Automation Conference (DAC)},
  title={Bit-Slicing the Hilbert Space: Scaling Up Accurate Quantum Circuit Simulation},
  year={2021},
  pages={439-444},
  doi={10.1109/DAC18074.2021.9586191}
}
```

## Contact
If you have any questions or suggestions, feel free to [create an issue](https://github.com/NTU-ALComLab/SliQSim/issues), or contact us through matthewyhtsai@gmail.com.

## Project : SliQSim on QAOA

Thanks to [homuch](https://github.com/homuch) for implementing RUS synthesis, merging RUS gates into SliQSim, and adding multiprocessing in the QAOA.py

For the detail of RUS usage, see [https://github.com/homuch/rus](https://github.com/homuch/rus)

First, make and compile the SliQSim executable.

Then run QAOA.py (need pennylane package) (to test the code you can set use_outer_simulator to False)

To execute python file, type
```commandline
python QAOA.py <use_outer_simulator> <simulator index> <case number>
```
For example, **python QAOA.py 1 0 0** will use the outer simulator 0 (SliQSim) to run Case 1

Currently support simulator

* qml.device 0 : defaulit.qubit
* qml.device 1 : qiskit.aer (need <= Qiskit 0.45.0)
* qml.device 2 : lightning.gpu (need to find a way to download plugin)
* Outer simulator 0 : SliQSim
* Outer simulator 1 : DDSIM
* Outer simulator 2 : SliQSim_RUS

Currently working :

1. (Implemented!!!) Combine SliQSim and RUS : We wish to deal with the following line in the QAOA qasm file : **rz(0.1) q[0];**. We will precompile the RUS circuits in the folder */rus/precompile* first. Then there are at least two ways to combine them. The first way is
directly read in the command rz(0.1) q[0]; and the SliQSim code can find the nearest angle among the precompile file and perform RUS for it (You may add the option --rus to indicate the use of RUS).
The second way is find the nearest angle among the precompile file in the QAOA.py python code and write the corresponding circuit into the qasm file that will be executed by SliQSim.

2. (Probably solved) Find the way to update parameters with SliQSim : Currently the performance of SliQSim is bad (not getting good solution). Need to find a good way to update parameters (such as the way to compute gradient, the choice of the stepsize etc.)

3. Maybe run some circuits with arbitrary angle RZ gates ?