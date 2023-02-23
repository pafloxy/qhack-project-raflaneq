###########################################################################################
## IMPORTS ##
###########################################################################################
import numpy as np
from typing import Optional
from tqdm import tqdm
from collections import Counter
from .basic_utils import qsm, states, MCMCChain, MCMCState
# from .prob_dist import *
from .energy_models import IsingEnergyFunction
from .classical_mcmc_routines import test_accept, get_random_state
from qiskit import (
    QuantumCircuit,
    QuantumRegister,
    ClassicalRegister,
    execute,
)
from qiskit.extensions import UnitaryGate, XGate, ZGate, HamiltonianGate
################################################################################################
##  QUANTUM CIRCUIT CONSTRUCTION ##
################################################################################################


def initialise_qc(n_spins: int, bitstring: str) -> QuantumCircuit :
    """
    Initialises a quantum circuit with n_spins number of qubits in a state defined by "bitstring"
    (## NOTE : Qiskit's indexing convention for qubits (order of tensor product) is different from the conventional textbook one!)
    
    """

    spins = QuantumRegister(n_spins, name="spin")
    creg_final = ClassicalRegister(n_spins, name="creg_f")
    qc_in = QuantumCircuit(spins, creg_final)

    len_str_in = len(bitstring)
    assert len_str_in == len(
        qc_in.qubits
    ), "len(bitstring) should be equal to number_of_qubits/spins"

    # print("qc_in.qubits: ", qc_in.qubits)
    where_x_gate = [
        qc_in.qubits[len_str_in - 1 - i]
        for i in range(0, len(bitstring))
        if bitstring[i] == "1"
    ]
    if len(where_x_gate) != 0:
        qc_in.x(where_x_gate)
    return qc_in


def fn_qc_h1(num_spins: int, gamma, alpha, h, delta_time) -> QuantumCircuit :
    """
    Create a Quantum Circuit for time-evolution under
    hamiltonain H1 (described in the paper) 

    ARGS:
    ----
    num_spins: number of spins in the model
    gamma: float
    alpha: float
    h: list of field at each site
    delta_time: total evolution time time/num_trotter_steps
    """
    a = gamma
    # print("a:",a)
    b_list = [-(1 - gamma) * alpha * hj for hj in h]
    list_unitaries = [
        UnitaryGate(
            HamiltonianGate(
                a * XGate().to_matrix() + b_list[j] * ZGate().to_matrix(),
                time=delta_time,
            ).to_matrix(),
            label=f"exp(-ia{j}X+b{j}Z)",
        )
        for j in range(0, num_spins)
    ]
    qc = QuantumCircuit(num_spins)
    for j in range(0, num_spins):
        qc.append(list_unitaries[j], [num_spins - 1 - j])
    qc.barrier()
    # print("qc is:"); print(qc.draw())
    return qc


def fn_qc_h2(J:np.array, alpha:float, gamma:float, delta_time=0.8) -> QuantumCircuit :
    """
    Create a Quantum Circuit for time-evolution under
    hamiltonain H2 (described in the paper)

    ARGS:
    ----
    J: interaction matrix, interaction between different spins
    gamma: float
    alpha: float
    delta_time: (default=0.8, as suggested in the paper)total evolution time time/num_trotter_steps
    """
    num_spins = np.shape(J)[0]
    qc_for_evol_h2 = QuantumCircuit(num_spins)
    # calculating theta_jk
    upper_triag_without_diag=np.triu(J,k=1)
    theta_array=(-2*(1-gamma)*alpha*delta_time)*upper_triag_without_diag
    for j in range(0, num_spins - 1):
        for k in range(j+1,num_spins):
            angle=theta_array[j,k]
            qc_for_evol_h2.rzz(
                angle, qubit1=num_spins - 1 - j, qubit2=num_spins - 1 - k
            )
    # print("qc for fn_qc_h2 is:"); print(qc_for_evol_h2.draw())
    return qc_for_evol_h2


def trottered_qc_for_transition(num_spins: int, qc_h1: QuantumCircuit, qc_h2: QuantumCircuit, num_trotter_steps: int) -> QuantumCircuit:
    """ Returns a trotter circuit (evolution_under_h2 X evolution_under_h1)^(r-1) (evolution under h1)"""
    qc_combine = QuantumCircuit(num_spins)
    for i in range(0, num_trotter_steps - 1):
        qc_combine = qc_combine.compose(qc_h1)
        qc_combine = qc_combine.compose(qc_h2)
        qc_combine.barrier()
    qc_combine = qc_combine.compose(qc_h1)
    # print("trotter ckt:"); print(qc_combine.draw())
    return qc_combine


def combine_2_qc(init_qc: QuantumCircuit, trottered_qc: QuantumCircuit) -> QuantumCircuit:
    """ Function to combine 2 quantum ckts of compatible size.
        In this project, it is used to combine initialised quantum ckt and quant ckt meant for time evolution
    """
    num_spins = len(init_qc.qubits)
    qc = QuantumCircuit(num_spins, num_spins)
    qc = qc.compose(init_qc)
    qc.barrier()
    qc = qc.compose(trottered_qc)
    return qc


################################################################################################
##  QUANTUM MARKOV CHAIN CONSTRUCTION ##
################################################################################################

def run_qc_quantum_step(
    qc_initialised_to_s: QuantumCircuit, model: IsingEnergyFunction, alpha, n_spins: int
) -> str:

    """
    Takes in a qc initialized to some state "s". After performing unitary evolution U=exp(-iHt)
    , circuit is measured once. Function returns the bitstring s', the measured state .

    ARGS:
    ----
    qc_initialised_to_s:
    model:
    alpha:
    n_spins:
    
    """

    h = model.get_h# and not model.get_h() anymore
    J = model.get_J# and not model.get_J() anymore

    # init_qc=initialise_qc(n_spins=n_spins, bitstring='1'*n_spins)
    gamma = np.round(np.random.uniform(0.25, 0.6), decimals=2)
    time = np.random.choice(list(range(2, 12)))  # earlier I had [2,20]
    delta_time = 0.8
    num_trotter_steps = int(np.floor((time / delta_time)))
    # print(f"gamma:{gamma}, time: {time}, delta_time: {delta_time}, num_trotter_steps:{num_trotter_steps}")
    # print(f"num troter steps: {num_trotter_steps}")
    qc_evol_h1 = fn_qc_h1(n_spins, gamma, alpha, h, delta_time)
    qc_evol_h2 = fn_qc_h2(J, alpha, gamma, delta_time=delta_time)
    trotter_ckt = trottered_qc_for_transition(
        n_spins, qc_evol_h1, qc_evol_h2, num_trotter_steps=num_trotter_steps
    )
    qc_for_mcmc = combine_2_qc(qc_initialised_to_s, trotter_ckt)

    # run the circuit
    num_shots = 1
    quantum_registers_for_spins = qc_for_mcmc.qregs[0]
    classical_register = qc_for_mcmc.cregs[0]
    qc_for_mcmc.measure(quantum_registers_for_spins, classical_register)
    # print("qc_for_mcmc: ")
    # print( qc_for_mcmc.draw())
    state_obtained_dict = (
        execute(qc_for_mcmc, shots=num_shots, backend=qsm).result().get_counts()
    )
    state_obtained = list(state_obtained_dict.keys())[
        0
    ]  # since there is only one element
    return state_obtained


def quantum_enhanced_mcmc(
    n_hops: int,
    model: IsingEnergyFunction,
    # alpha,
    initial_state: Optional[str] = None,
    temperature=1,
):
    """
    version 0.2
    
    ARGS:
    ----
    Nhops: Number of time you want to run mcmc
    model:
    return_last_n_states:
    return_both:
    temp:

    RETURNS:
    -------
    Last 'return_last_n_states' elements of states so collected (default value=500). one can then deduce the distribution from it!
    
    """
    num_spins = model.num_spins

    if initial_state is None:
        initial_state = MCMCState(get_random_state(num_spins), accepted=True)
    else:
        initial_state = MCMCState(initial_state, accepted=True)
    
    current_state: MCMCState = initial_state
    energy_s = model.get_energy(current_state.bitstring)
    print("starting with: ", current_state.bitstring, "with energy:", energy_s)

    mcmc_chain = MCMCChain([current_state])

    print(mcmc_chain)
    for _ in tqdm(range(0, n_hops), desc='runnning quantum MCMC steps . ..' ):
        # get sprime
        qc_s = initialise_qc(n_spins= model.num_spins, bitstring=current_state.bitstring)
        s_prime = run_qc_quantum_step(
            qc_initialised_to_s=qc_s, model=model, alpha=model.alpha, n_spins= model.num_spins
        )
        # accept/reject s_prime
        energy_sprime = model.get_energy(s_prime)
        accepted = test_accept(
            energy_s, energy_sprime, temperature=temperature
        )
        mcmc_chain.add_state(MCMCState(s_prime, accepted))
        if accepted:
            current_state = mcmc_chain.current_state
            energy_s = model.get_energy(current_state.bitstring)

    return mcmc_chain 