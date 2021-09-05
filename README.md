# UndergraduateThesis
This is my undergraduate thesis project.

In this project, I will calculate the high-order harmonic generations(HHG) and above-threshold ionization spectra using C++.
I will also try parallel computing using OpenMp and MPI in this project.

(1).Source1.f90: Calculate the matrix of H0 and output to file "H0_3.txt".
(2).EignVal.f90: Take "H0_3.txt" as input, then using Jacob Iteration method to calculate its eigenvalues and eigenvectors, i.e. the eigen energies and eigen wave functions of field free hamiltonian H0. The eigenvalues are recorded in output file with user defined name.
(3).OperatorA_cmplx.f90: This file takes the result of "EignVal.f90" as input, then uses eigenvector-expansion method to get the matrix of A=exp(-i*H0*dt/2).
(4).OperatorB_cmplx.f90: This file caulculate the matrix of B=exp(-i*HI(t+dt/2)*dt).
(5).TimeProp_cmplx.f90: This file takes the result of 3&4 as input, then performs the time evolution of using groudstate of H0 as initial wavefunction. Specifically, using iteration formula: Psi(t+dt)=(A*B*A)*Psi(t) to evolve. The velocity-form dipole dv(t) is also calculated.
(6).GetHHG_cmplx.f90: Calculate the HHG spectrum by performing fourier transformation to dv(t).
