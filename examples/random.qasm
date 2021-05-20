OPENQASM 2.0;
include "qelib1.inc";
qreg q[25];
creg c[25];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];
h q[9];
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
h q[15];
h q[16];
h q[17];
h q[18];
h q[19];
h q[20];
h q[21];
h q[22];
h q[23];
h q[24];
cz q[3],q[10];
x q[7];
rx(pi/2) q[16];
rx(pi/2) q[21];
cz q[10],q[1];
x q[15];
cx q[2],q[24];
rx(pi/2) q[4];
t q[8];
y q[18];
cz q[9],q[8];
cz q[10],q[21];
cz q[23],q[5];
t q[10];
h q[16];
rx(pi/2) q[22];
ry(pi/2) q[10];
y q[12];
ccx q[23],q[24],q[19];
z q[14];
ry(pi/2) q[15];
t q[8];
t q[9];
rx(pi/2) q[14];
x q[13];
cz q[3],q[23];
s q[17];
cswap q[13],q[6],q[7];
ccx q[3],q[22],q[1];
cx q[5],q[13];
z q[4];
z q[16];
ccx q[10],q[17],q[0];
ry(pi/2) q[23];
ry(pi/2) q[18];
cx q[8],q[16];
ccx q[0],q[17],q[13];
x q[17];
y q[22];
s q[1];
z q[18];
cz q[22],q[12];
cswap q[1],q[13],q[12];
z q[5];
ccx q[10],q[12],q[3];
ry(pi/2) q[5];
cx q[6],q[2];
cswap q[1],q[20],q[22];
t q[7];
h q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
measure q[8] -> c[8];
measure q[9] -> c[9];
measure q[10] -> c[10];
measure q[11] -> c[11];
measure q[12] -> c[12];
measure q[13] -> c[13];
measure q[14] -> c[14];
measure q[15] -> c[15];
measure q[16] -> c[16];
measure q[17] -> c[17];
measure q[18] -> c[18];
measure q[19] -> c[19];
measure q[20] -> c[20];
measure q[21] -> c[21];
measure q[22] -> c[22];
measure q[23] -> c[23];
measure q[24] -> c[24];
 