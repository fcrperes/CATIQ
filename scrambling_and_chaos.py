import random
import json
import numpy as np
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, execute, Aer


def random_circuits(nr_qubits, nr_layers, gate_set):
    """Function that evolves the state of a randomly created quantum
    circuit with nr_qubits qubits and nr_layers layers, and gates
    drawn from a gate set defined by gate_set.

    Args:
        nr_qubits (int): Number of qubits in the quantum circuit
        nr_layers (int): Number of layers responsible for the state
        evolution
        gate_list (List[str]): List indicating the chosen gate set

    Returns:
        QuantumCircuit: The corresponding quantum circuit
    """
    qc = QuantumCircuit(nr_qubits)

    # Apply the gates in the gate list to the quantum circuit
    for _ in range(nr_layers):
        used_qubits = []
        unused_qubits = list(range(nr_qubits))
        for qubit in range(nr_qubits):
            if qubit in used_qubits:
                continue

            unused_qubits.remove(qubit)
            # Choose a random gate from the gate list
            gate = random.choice(gate_set)
            if len(used_qubits) == nr_qubits - 1:
                while gate == 'CNOT':
                    gate = random.choice(gate_set)

            if gate == 'H':
                qc.h(qubit)
            elif gate == 'T':
                qc.t(qubit)
            elif gate == 'X':
                qc.x(qubit)
            elif gate == 'Y':
                qc.y(qubit)
            elif gate == 'Z':
                qc.z(qubit)
            elif gate == 'S':
                qc.s(qubit)
            elif gate == 'CNOT':
                coin = random.randint(0, 1)
                if coin == 0:
                    ctrl = qubit
                    targ = ctrl + 1
                    used_qubits.append(targ)
                    unused_qubits.remove(targ)
                else:
                    targ = qubit
                    ctrl = targ + 1
                    used_qubits.append(ctrl)
                    unused_qubits.remove(ctrl)
                qc.cx(ctrl, targ)

            used_qubits.append(qubit)

    return qc


def hadamard_layer(nr_qubits, randh):
    """Function that creates a Hadamard layer.

    Args:
        nr_qubits (int): Number of qubits in the quantum circuit
        randh (bool): A flag that determines if the Hadamard layer
        is complete or if there's a random number of gates in random
        locations. If randh = True, the layer is incomplete; other-
        wise, the layer is complete.

    Returns:
        QuantumCircuit: The corresponding quantum circuit
    """
    qc = QuantumCircuit(nr_qubits)

    if randh:
        had_count = random.randint(1, nr_qubits)
        for qb in random.sample(list(range(nr_qubits)), had_count):
            qc.h(qb)

    else:
        for qb in range(nr_qubits):
            qc.h(qb)

    return qc


def phase_layer(nr_qubits, rands):
    """Function that creates a layer with phase gates.

    Args:
        nr_qubits (int): Number of qubits in the quantum circuit
        rands (bool): A flag that determines if the phase layer is
        complete or if there's a random number of gates in random
        locations. If rands = True, the layer is incomplete; other-
        wise, the layer is complete.

    Returns:
        QuantumCircuit: The corresponding quantum circuit
    """
    qc = QuantumCircuit(nr_qubits)

    if rands:
        had_count = random.randint(1, nr_qubits)
        for qb in random.sample(list(range(nr_qubits)), had_count):
            qc.s(qb)

    else:
        for qb in range(nr_qubits):
            qc.s(qb)

    return qc


def kron_product(nr_qubits, Op):
    """Function for making the scrambling operators.

    Args:
        nr_qubits (int): Number of qubits
        Op (np.array): The numpy array defining a given operator

    Returns:
        np.array: A numpy array which describes the scrambled
        version of Op
    """
    Identity = np.eye(Op.shape[0])
    result = Op
    for _ in range(nr_qubits - 1):
        result = np.kron(result, Identity)
    return result


def otoc(U_op, V_op, W_op, dens_mat):
    """Function to calculate the OTOC

    Args:
        U_op (np.array): Circuit unitary
        V_op (np.array): V of the otoc
        W_op (np.array): W of the otoc (butterfly operator)
        dens_mat (np.array): Circuit input state

    Returns:
        complex: The complex-valued F(t)
    """
    # Evolution of the operator W under the unitary dynamics U
    W_t = np.dot(np.dot(np.transpose(np.conjugate(U_op)), W_op), U_op)

    # OTOC
    F_t = np.trace(
        np.transpose(np.conjugate(W_t)) @ np.transpose(np.conjugate(V_op))
        @ W_t @ V_op @ dens_mat)
    return F_t


def execute_code(nr_qubits,
                 nr_layers,
                 nr_circuits,
                 gate_set,
                 t_steps,
                 inj_times,
                 location,
                 rand=True,
                 flagh=True,
                 flags=True):
    """Main function that runs the entire code.

    Args:
        nr_qubits (int): Number of qubits
        nr_layers (int): Number of layers between each otoc evaluation
        nr_circuits (int): Number of circuits
        gate_set (list["str"]): Desired gate set
        t_steps (int): Number of time steps to consider (i.e. number
        of otoc evaluations)
        inj_times (list[int]): Injection times (i.e., times where we
        inject the cat gadget)
        location (str): Path to where we store the data and plots
        rand (bool, opt): If we are to use randomly generated layers
        or only pseudo-random layers
        flagh (bool, opt): If rand is False, then we need to decide
        if each of the non-randomly chosen Hadamrd layers are created
        in a random way (if True) or are complete (if False)
        flags (bool, opt): If rand is False, then we need to decide
        if each of the non-randomly chosen phase layers are created
        in a random way (if True) or are complete (if False)
    """
    dimension_HS = 2**nr_qubits
    with open(f'{location}/Time_steps.txt', 'w') as data_file:
        json.dump(list(range(t_steps)), data_file)

    # Define the initial state rho
    rho = np.zeros((dimension_HS, dimension_HS))
    rho[0, 0] = 1

    # Define the operators V and W
    PauliX = np.array([[0, 1], [1, 0]])
    PauliZ = np.array([[1, 0], [0, -1]])
    V = kron_product(nr_qubits, PauliX)
    W = np.kron(kron_product(nr_qubits - 1, np.eye(2)), PauliZ)
    # For making a single plot
    F_total_re = []
    F_total_im = []

    for circ in range(nr_circuits):
        # Create two quantum circuit variables:
        circuit = QuantumCircuit(nr_qubits)
        block = QuantumCircuit(nr_qubits)
        F_values_re = []
        F_values_im = []
        save_circ = ""
        if rand:
            for t in range(t_steps):
                if t == 0:
                    # Execute the circuit on a simulator
                    backend = Aer.get_backend('unitary_simulator')
                    job = execute(circuit, backend)
                    U = job.result().get_unitary()
                    F_values_re.append(otoc(U, V, W, rho).real)
                    F_values_im.append(otoc(U, V, W, rho).imag)
                    continue

                if t in inj_times:
                    central_qubit = random.randint(0, nr_qubits - 1)
                    if central_qubit >= int(nr_qubits / 2):
                        for j in range(0, central_qubit+1):
                            circuit.t(j)
                        for j in range(0, central_qubit):
                            circuit.cx(j, j + 1)
                        circuit.t(central_qubit)
                        for j in range(central_qubit, 1, -1):
                            circuit.cx(j - 1, j)
                    else:
                        for j in range(central_qubit,nr_qubits):
                            circuit.t(j)
                        for j in range(nr_qubits - 1, central_qubit, -1):
                            circuit.cx(j, j - 1)
                        circuit.t(central_qubit)
                        for j in range(central_qubit, nr_qubits - 1):
                            circuit.cx(j + 1, j)
                    save_circ += circuit.qasm()
                else:
                    block = random_circuits(nr_qubits, nr_layers, gate_set)
                    save_circ += block.qasm()
                    circuit.append(block, [*range(nr_qubits)])

                # Execute the circuit on a simulator
                backend = Aer.get_backend('unitary_simulator')
                job = execute(circuit, backend)
                U = job.result().get_unitary()
                F_values_re.append(otoc(U, V, W, rho).real)
                F_values_im.append(otoc(U, V, W, rho).imag)

        else:
            which = 0
            for t in range(t_steps):
                if t == 0:
                    # Execute the circuit on a simulator
                    backend = Aer.get_backend('unitary_simulator')
                    job = execute(circuit, backend)
                    U = job.result().get_unitary()
                    F_values_re.append(otoc(U, V, W, rho).real)
                    F_values_im.append(otoc(U, V, W, rho).imag)
                    continue

                if t in inj_times:
                    for j in range(nr_qubits):
                        circuit.t(j)
                    central_qubit = random.randint(0, nr_qubits - 1)

                    if central_qubit >= int(nr_qubits / 2):
                        for j in range(0, central_qubit):
                            circuit.cx(j, j + 1)
                        circuit.t(central_qubit)
                        for j in range(central_qubit, 1, -1):
                            circuit.cx(j - 1, j)
                    
                    else:
                        for j in range(nr_qubits - 1, central_qubit, -1):
                            circuit.cx(j, j - 1)
                        circuit.t(central_qubit)
                        for j in range(central_qubit, nr_qubits - 1):
                            circuit.cx(j + 1, j)
                    
                    save_circ += circuit.qasm()

                else:
                    if which % 4 == 0:
                        block = hadamard_layer(nr_qubits, flagh)

                    else:
                        block = phase_layer(nr_qubits, flags)

                    save_circ += block.qasm()
                    circuit.append(block, [*range(nr_qubits)])
                    which += 1

                # Execute the circuit on a simulator
                backend = Aer.get_backend('unitary_simulator')
                job = execute(circuit, backend)
                U = job.result().get_unitary()
                F_values_re.append(otoc(U, V, W, rho).real)
                F_values_im.append(otoc(U, V, W, rho).imag)

        # Saving the data for each circuit:
        with open(f'{location}/Circuit{circ}.txt', 'w') as circ_file:
            circ_file.write(save_circ)
        with open(f'{location}/Circuit{circ}-ReF.txt', 'w') as data_file:
            json.dump(F_values_re, data_file)
        with open(f'{location}/Circuit{circ}-ImF.txt', 'w') as data_file:
            json.dump(F_values_im, data_file)
        # Append F_values_re and _im into total
        F_total_re.append(F_values_re)
        F_total_im.append(F_values_im)

        # Plotting 
        fig = plt.figure(figsize=(12, 8))
        plt.title('Out-of-time-order correlator', fontsize=20)
        plt.xlabel(r"$t$", fontsize=14)
        plt.ylabel(r"$F(t)$", fontsize=14)
        plt.ylim(-1.05, 1.05)

        plt.plot(list(range(t_steps)),
                 F_values_re,
                 '-o',
                 color='blue',
                 label=r"$\mathfrak{R}(F(t))$")
        plt.plot(list(range(t_steps)),
                 F_values_im,
                 '-o',
                 color='red',
                 markerfacecolor='white',
                 label=r"$\mathfrak{I}(F(t))$")

        plt.legend(fontsize=14, loc="best")
        plt.tight_layout()
        fig.savefig(f'{location}/Plot-Circuit{circ}.pdf')
        plt.close()

    #Plotting the total plot
    fig = plt.figure(figsize=(12, 8))
    plt.title('Out-of-time-order correlator', fontsize=30)
    plt.xlabel(r"$\tau$", fontsize=30)
    plt.ylabel(r"$\mathfrak{R}[F(\tau)]$", fontsize=30)
    #plt.ylim(-0.3, 1.05)
    plt.ylim(-1.05, 1.05)
    # Set the size of the tick labels for both the x and y axes
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    highlight=0
    for F_values_re in F_total_re:
        if highlight == nr_circuits - 1:
            plt.plot(list(range(t_steps)),
                 F_values_re,
                 '-o',
                 color='blue')
        if highlight <= nr_circuits - 2:
            plt.plot(list(range(t_steps)),
                 F_values_re,
                 '-o',
                 color='#aec2d1')
            highlight+=1
    #plt.legend(fontsize=25, loc="best")
    plt.tight_layout()
    fig.savefig(f'{location}/Plot-Circuit-Total.pdf')
    plt.close()
