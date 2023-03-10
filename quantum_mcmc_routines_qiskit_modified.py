###########################################################################################
## IMPORTS ##
###########################################################################################
import numpy as np
from typing import Optional
from tqdm import tqdm
from collections import Counter

from qiskit import IBMQ

# IBMQ.save_account(TOKEN)   #Input your token before using
IBMQ.load_account()
provider = IBMQ.get_provider(hub ='ibm-q')
backend = provider.get_backend('ibmq_lima')  ## !!!Remember to change provider later



## Removed just for testing
# from .basic_utils import qsm, states, MCMCChain, MCMCState
# # from .prob_dist import *
# from .energy_models import IsingEnergyFunction
# from .classical_mcmc_routines import test_accept, get_random_state

from qiskit import Aer
qsm = Aer.get_backend("qasm_simulator") ## Added for testing

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

## This function is not necessary now
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


def fn_qc_h2(J:np.array, alpha:float, gamma:float, delta_time=0.8, RZZ_twirl = False) -> QuantumCircuit :
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
    def is_commute(i,j):
        ## checks if twirling Paulis commute (c = 1) or anti-commute (c = -1)
        Z = 3
        c = 1
        if i == 0 or i == Z:
            c = c*(-1)
        
        if j == 0 or j == Z:
            c = c*(-1)
        
        return c ## = -1 if Paulis anti-commute, +1 is Paulis commute

    def pauli(qc:QuantumCircuit, index:int ,qubit:int):
        ## Checks which Pauli to apply for twirling
        if index == 1:
            qc.x(qubit)
        elif index == 2:
            qc.y(qubit)
        elif index == 3:
            qc.z(qubit)

    num_spins = np.shape(J)[0]
    qc_for_evol_h2 = QuantumCircuit(num_spins)
    # calculating theta_jk
    upper_triag_without_diag=np.triu(J,k=1)
    theta_array=(-2*(1-gamma)*alpha*delta_time)*upper_triag_without_diag

    if RZZ_twirl:
        for j in range(0, num_spins - 1):
            for k in range(j+1,num_spins):
                ## Pauli twirling for RZZ
                Pauli_a,Pauli_b = np.random.choice(range(3), size=2)    ## Choose two Pauli from {I,X,Y,Z}
                pauli(qc_for_evol_h2,Pauli_a,num_spins - 1 - j)
                pauli(qc_for_evol_h2,Pauli_b,num_spins - 1 - k)    

                ## Flip theta is Paulis anti-commute, else remains same
                check = is_commute(Pauli_a,Pauli_b)
                angle=check * theta_array[j,k]  
                qc_for_evol_h2.rzz(
                    angle, qubit1=num_spins - 1 - j, qubit2=num_spins - 1 - k
                )

                ## Finish twirling 
                pauli(qc_for_evol_h2,Pauli_a,num_spins - 1 - j)
                pauli(qc_for_evol_h2,Pauli_b,num_spins - 1 - k)  

    else:
        ## No twirling
        for j in range(0, num_spins - 1):
            for k in range(j+1,num_spins):
                angle= theta_array[j,k]  
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
    qc = QuantumCircuit(num_spins, num_spins*2)
    qc = qc.compose(init_qc)
    qc.barrier()

    qc = init_qc.compose(trottered_qc)
    return qc


################################################################################################
##  QUANTUM MARKOV CHAIN CONSTRUCTION ##
################################################################################################
# model: IsingEnergyFunction,


## Testing parameters
n_spins = 2
shape_of_J=(n_spins,n_spins)

## defining J matrix (mutual 1-1 interaction)
# J =  np.round(np.random.choice([+1, 0, -1], size=(n_spins, n_spins)), decimals=2) 
J =  np.random.uniform(low= -2, high= 2, size= shape_of_J )

J = 0. * (J + J.transpose() )
J = np.round( J - np.diag(np.diag(J)) , decimals= 3)

# defining h
h = np.round(0.5 * np.random.randn(n_spins), decimals=2)


def run_qc_quantum_step(h,J, alpha, n_spins: int,SPAM_twirl:bool = True, RZZ_twirl:bool = False
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
    if SPAM_twirl:
        c = np.random.choice([-1, 1], size=n_spins) ## For SPAM twirling, added into J,h
    else:
        c = np.ones(n_spins)    ## No spins are flipped

    c_ij = np.outer(c,c)
    flipping = int("0b" + "".join('1' if x == -1 else '0' for x in reversed(c)),2) ## Used later for flipping inital/final bits
                                                                                 ## Reversed c due to opposite encoding in qiskit
    h = c * h#* model.get_h# and not model.get_h() anymore
    J = c_ij * J #* model.get_J# and not model.get_J() anymore

    # init_qc=initialise_qc(n_spins=n_spins, bitstring='1'*n_spins)
    gamma = np.round(np.random.uniform(0.25, 0.6), decimals=2)
    time = np.random.choice(list(range(2, 12)))  # earlier I had [2,20]
    delta_time = 0.8
    num_trotter_steps = int(np.floor((time / delta_time)))
    # print(f"gamma:{gamma}, time: {time}, delta_time: {delta_time}, num_trotter_steps:{num_trotter_steps}")
    # print(f"num troter steps: {num_trotter_steps}")
    qc_evol_h1 = fn_qc_h1(n_spins, gamma, alpha, h, delta_time)
    qc_evol_h2 = fn_qc_h2(J, alpha, gamma, delta_time=delta_time,RZZ_twirl = RZZ_twirl)
    trotter_ckt = trottered_qc_for_transition(
        n_spins, qc_evol_h1, qc_evol_h2, num_trotter_steps=num_trotter_steps
    )

    initial_state = ClassicalRegister(n_spins)
    final_state = ClassicalRegister(n_spins)
    quantum_registers_for_spins = QuantumRegister(n_spins)

    qc_initialized = QuantumCircuit(quantum_registers_for_spins, initial_state,final_state)
    qc_initialized.h(range(n_spins))    ## Creating equal superpositon
    qc_initialized.measure(quantum_registers_for_spins,initial_state) ## Chooses initial state


    qc_for_mcmc = combine_2_qc(qc_initialized, trotter_ckt)


    # run the circuit
    num_shots = 2**(n_spins*6)  ## Chosen to take enough samples, can be modified later


    qc_for_mcmc.measure(qc_for_mcmc.qregs[0], final_state)
    # print("qc_for_mcmc: ")
    # print( qc_for_mcmc.draw())

    state_obtained_dict = (
        execute(qc_for_mcmc, shots=num_shots, backend=backend).result().get_counts()
    )

    proposal_matrix = np.zeros((2**n_spins,2**n_spins)) ## Stores proposal probability 

    for states in state_obtained_dict:
        input_state = int(states[n_spins:],2)^flipping  ## From first measurement
        output_state = int(states[:n_spins],2) ^flipping  ## From second measurement (Corrected for twirling using XOR operator ^)
        transition = state_obtained_dict[states]/num_shots

        proposal_matrix[input_state,output_state] = transition *(2**n_spins) ## Factor for normalizing

    return proposal_matrix


print(np.array_str(run_qc_quantum_step(h,J,0.8,n_spins,True,False), precision=3))

# def quantum_enhanced_mcmc(
#     n_hops: int,
#     model: IsingEnergyFunction,
#     # alpha,
#     initial_state: Optional[str] = None,
#     temperature=1,
# ):
#     """
#     version 0.2
    
#     ARGS:
#     ----
#     Nhops: Number of time you want to run mcmc
#     model:
#     return_last_n_states:
#     return_both:
#     temp:

#     RETURNS:
#     -------
#     Last 'return_last_n_states' elements of states so collected (default value=500). one can then deduce the distribution from it!
    
#     """
#     num_spins = model.num_spins

#     if initial_state is None:
#         initial_state = MCMCState(get_random_state(num_spins), accepted=True)
#     else:
#         initial_state = MCMCState(initial_state, accepted=True)
    
#     current_state: MCMCState = initial_state
#     energy_s = model.get_energy(current_state.bitstring)
#     print("starting with: ", current_state.bitstring, "with energy:", energy_s)

#     mcmc_chain = MCMCChain([current_state])

#     print(mcmc_chain)
#     for _ in tqdm(range(0, n_hops), desc='runnning quantum MCMC steps . ..' ):
#         # get sprime
#         qc_s = initialise_qc(n_spins= model.num_spins, bitstring=current_state.bitstring)
#         s_prime = run_qc_quantum_step(
#             qc_initialised_to_s=qc_s, model=model, alpha=model.alpha, n_spins= model.num_spins
#         )
#         # accept/reject s_prime
#         energy_sprime = model.get_energy(s_prime)
#         accepted = test_accept(
#             energy_s, energy_sprime, temperature=temperature
#         )
#         mcmc_chain.add_state(MCMCState(s_prime, accepted))
#         if accepted:
#             current_state = mcmc_chain.current_state
#             energy_s = model.get_energy(current_state.bitstring)

#     return mcmc_chain 