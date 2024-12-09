% This matlab code studies the Lindblad master equation defining
% interaction between atom and cavity in terms of the dissipative
% atom-field dynamics. 
% 
% Time evolution is obtained with the 2nd-order Runge-Kutta approach.
% 
% Atom and field hamiltonian is the Jaynes-Cummings (JC) model:
% H_JC = 0.5 * hbar * wa * sigma_z + hbar * wf * a^{+} * a + hbar * g * (sigma_minus * a^{+} + sigma_{+} * a); 
%
% The Lindblad master equation:
%   drho/dt = -(i/hbar)*[H_JC, rho] + L_cavity(rho) + L_atom(rho)
% where, 
%   L_cavity(rho) = 0.5 * kappa * (2 * a * rho * a^{+} - rho * a^{+} * a - a^{+} * a * rho)
%   L_atom(rho) = 0.5 * gamma * (2 * sigma_minus * rho * sigma_plus - rho * sigma_plus * sigma_minus - sigma_plus * sigma_minus * rho)
%
% An initial state of the field: a coherent state and atomic excited state.
% 
% Output: P_ex_st(t) = <es_st | rho(t) | ex_st >, where rho = |psi> <psi|. 
%
%  References: [1] N. Meher J. Phys. B: At. Mol. Opt. Phys. 57 073001 (2024); 
%              [2] H. J. Carmichael, Statistical Methods in Quantum Optics 1: 
%                  Master Equations and Fokker-Planck Equations (Springer, 1999).
%              [3] G. Lindblad, Communications in Mathematical Physics 48, 119 (1976). 
% 
% The atomic unit (au) is used in the calculation. 
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
% Dec 6, 2024 & University of North Dakota 
%
function [] = Lindblad_master_equation_1
%
clear; clc; format short
%
n_max = 20.;
%
a_ann = diag(sqrt(1:n_max-1),1);  % a          & annihilation operator for EM field
a_cr = a_ann';                    % a^{dagger} & creation operator for EM field
N_op = a_cr * a_ann;              % N = a^{dagger} * a
%
[Vec,En] = eig(N_op);                                     % Eigenvalue problem: N|n> = n|n>
En = diag(En);
[foo, ij] = sort(En);
En = En(ij);
[En(1),En(2),En(3),En(4),En(5)] % eigenvalues: 0, 1, 2, ...
% [En(1),En(2),En(3),En(4),En(5)]
% 0    1.0000    2.0000    3.0000    4.0000
% E0 = 0 + 0.5; 
%
Vacuum = Vec(:,1);
n1 = Vec(:,2);    % one-photon state & |1 >
n2 = Vec(:,3);    % two-photon state & |2 >
n3 = Vec(:,4);    % three-photon state & |3 >

Vec(:,1:5);
%%%

%%%%%%%%%%%%%%%
hbar = 1.;
wa = 1.0; % atomic frequency
wf = 1.0; % EM field frequency
g = 0.20; % coupling constant
kappa = 0.050; % cavity decay rate
gamma = 0.020; % atomic decay rate
%
sigma_z = [1, 0;
           0, -1];
%
sigma_plus = [0, 1;
              0, 0];
%
sigma_minus = [0, 0;
               1, 0];
%
gr_st = [0;
         1]; % ground state
ex_st = [1;
         0]; % excited state
%
Ia = eye(2); % identity matrix for the atom
If = eye(n_max); % identity matrix for the EM field
%
H_atom = 0.5 * hbar * wa * kron(sigma_z, If); % atomic Hamiltonian
H_field = hbar * wf * kron(Ia, a_cr * a_ann); % EM field Hamiltonian
H_interaction = hbar * g * (kron(sigma_plus, a_ann) + kron(sigma_minus, a_cr)); % interaction Hamiltonian
%
H_total_JC = H_atom + H_field + H_interaction; % total JC Hamiltonian
%
n_ph = 1.; % initial number of photons in the cavity
%
psi = kron(ex_st, If(:,n_ph)); % atom in excited state and cavity in vacuum
%
rho = psi * psi';  % initial density matrix 
%
dt = 0.10; % time step
ci = sqrt(-1.);
%
T_final = 250.;
T_time = 0:dt:T_final;
%
an_op = kron(Ia, a_ann);
acr_op = kron(Ia, a_cr);
%
sigma_minus_op = kron(If, sigma_minus);
sigma_plus_op = kron(If, sigma_plus);
%
for t = 1:length(T_time)
    %
    Prob_ex_st(t) = conj(kron(ex_st, If(:,n_ph)))' * rho * kron(ex_st, If(:,n_ph)); % probability of atom in |ex_st>
    K1 = -ci*(H_total_JC * rho - rho * H_total_JC) + ...
              0.5 * kappa * (2 * an_op * rho * acr_op - acr_op * an_op * rho - rho * acr_op * an_op ) + ...
              0.5 * gamma * (2 * sigma_minus_op * rho * sigma_plus_op - sigma_plus_op * sigma_minus_op * rho - ...
                             rho * sigma_plus_op * sigma_minus_op); 
    rho1 = rho + dt * K1;
    %
    K2 = -ci*(H_total_JC * rho1 - rho1 * H_total_JC) + ...
              0.5 * kappa * (2 * an_op * rho1 * acr_op - acr_op * an_op * rho1 - rho1 * acr_op * an_op ) + ...
              0.5 * gamma * (2 * sigma_minus_op * rho1 * sigma_plus_op - sigma_plus_op * sigma_minus_op * rho1 - ...
                             rho1 * sigma_plus_op * sigma_minus_op);     
    %
    rho = rho + 0.5 * dt * (K1 + K2); 
    %
end
%

figure(1)
plot(T_time, abs(Prob_ex_st), 'b-', 'LineWidth',1.5)
xlabel('Time (au)')
ylabel('Probaility') % 
set(gca,'FontSize',16)
%axis([0 n_max, 0 0.18])
box on


%%%
return
end
