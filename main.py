import scrambling_and_chaos as scch
import random

####################################################################
######################## Lists of gate sets ########################
####################################################################
Paulis = ['X', 'Y', 'Z']
Incoherent = ['X', 'Y', 'Z', 'T', 'S', 'CNOT']

Clifford_separable = ['H', 'S']
Clifford_separable_plus_T = ['H', 'S', 'T']

Clifford = ['H', 'S', 'CNOT']
Clifford_plus_Pauli = ['H', 'S', 'CNOT','X', 'Y', 'Z']
Clifford_plus_T = ['H', 'S', 'CNOT', 'T']

####################################################################
###### Variable defintions (change only variables in here) #########
####################################################################
nq = 10  # number of qubits
nlayers = 10  # number of layes between each otoc evaluation
ncircs = 10  # number of simulated circuits
gset = Clifford_separable  #choose the gateset from list above
evals = 100  # unitary evolution steps (nr. of otoc evaluations)
# t_inject = []  # times to inject cat gadget
t_inject = random.sample(range(5,61), 36)
t_inject.sort()
save = f"/Users/fcrperes/Documents/3._Academic/PhD/Thesis/My_work/Cat_state_injection"
####################################################################

scch.execute_code(nq, nlayers, ncircs, gset, evals, t_inject, save)