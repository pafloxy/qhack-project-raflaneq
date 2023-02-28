# qhack-project-raflaneq
Repo for your official project submission ot the QHACK '23 hackathon.

Our aim is to reproduce the results from the paper by **'Quantum-Enhanced Markov Chain Monte Carlo'** by *David Layden* et al [on *arxiv*](https://arxiv.org/abs/2203.12497), which we belive fits under the **Quantum Computing Today** challenge. Firstly our aim is to reproduce the results directly from a real quantum device, following which we would like to advance by introducing hardware efficient modifications to the algorithm and obtain plausible enhancements than the former.

To this end we will be extensively using the python project [**quMCMC**](https://github.com/pafloxy/quMCMC) which was partly authored by one of our team members previously. This project provides simulation based implementation of the algorithm mentioned and some tests for the same.

The following documents/files are uploaded:
1. PresentationNotebook.ipynb :
This is the main document containing all the results from our project.

2. data-collecting-notebook.ipynb
Contains code that we used to interface with ibmq hardware (ibmq_guadalupe) and store the data in DATA/raw-circuit-outputs/ as .pickle file

3. QuantumSamplingRoutines.py
Contains some functions used to create and run the quantum circuits

4. qumcmc
Folder containing .py files necessary for running the classical and quantum enchanced mcmc on simulator, and extracting associated data, which was partly authored by one our team members previously. A simple running example is given in minimal_example.ipynb

5. minimal_example.ipynb
Contains a simple example showing how the qumcmc packages can be used

6. DATA
Folder conatining all the data obtained from circuits,simulations and model parameters