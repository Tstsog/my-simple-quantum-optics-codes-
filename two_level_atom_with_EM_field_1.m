% This matlab code studies the two-level atom interacting with an external
% electromagnetic (EM) field  in terms of the Jaynes-Cummings (JC)
% Hamiltonian, defined as:
%
% H_JC = 0.5 * hbar * wa * sigma_z + hbar * wf * a^{+} * a + hbar * g * (sigma_minus * a^{+} + sigma_{+} * a); 
%
% Probability of finding the atom in ground (gr) state (st) and the field in (n+1)-photon state, 
% P_{gr}(t): P_{gr}(r) = |<gr,n+1|psi(t)>|^{2} =  sin(g*t*sqrt(n+1))^2.
%
% Probability of finding the atom in excited (ex) state (st) and the field in (n)-photon state, 
% P_{ex}(t): P_{ex}(t) = |<ex,n|psi(t)>|^{2} =  cos(g*t*sqrt(n+1))^2.
%
%  References: [1] N. Meher J. Phys. B: At. Mol. Opt. Phys. 57 073001 (2024); 
%              [2] M. O. Scully, M. S. Zubairy, Quantum Optics, Cambridge University Presse (1997) 
% 
% The atomic unit (au) is used in the calculation. 
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
% Dec 6, 2024 & University of North Dakota 
%
function [] = two_level_atom_with_EM_field_1
%
clear; clc; format short
%
n_max = 35.; % dimension of number operator N
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
Vacuum = Vec(:,1); % vacuum 
n1 = Vec(:,2);     % one-photon state & |1 >
n2 = Vec(:,3);     % two-photon state & |2 >
n3 = Vec(:,4);     % three-photon state & |3 >
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Jaynes-Cummingg model for two-level atom interacting with the EM field. 
%
hbar = 1.;
wa = 1.50; % atomic frequency
wf = 1.50; % EM field frequency
g  = 0.20; % coupling constant
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
n_ph = 6.; % initial number of photons in the cavity
%
gr_st_wf = kron(gr_st, If(:,n_ph+2)); % atom in ground state and field has n+1 photons
ex_st_wf = kron(ex_st, If(:,n_ph+1)); % atom in excited state and field has n photons
%
psi = ex_st_wf; % initial state 
dt = 0.1; % time step
ci = sqrt(-1.);
%
U_evol_op = expm(-ci * H_total_JC * dt/hbar); % unitary evolution operator
%
T_final = 15.;
T_time = 0:dt:T_final;
%
for t = 1:length(T_time)
    %
    Prob_gr_st(t) = abs(gr_st_wf' * psi)^2; % 
    Prob_ex_st(t) = abs(ex_st_wf' * psi)^2;
    %
    psi = U_evol_op * psi;
    psi = psi/norm(psi);
    %
end
%
figure(1)
hold on
%
plot(T_time, Prob_gr_st, 'b-', 'LineWidth',1.5) % ground state
plot(T_time, Prob_ex_st, 'g-', 'LineWidth',1.5) % excited state
%
plot((0:2:T_final), sin(g.*(0:2:T_final).*sqrt(n_ph+1)).^2, 'bo', 'LineWidth',1.5) % exact = sin(g*t*sqrt(n+1))^2 for ground state
plot((0:2:T_final), cos(g.*(0:2:T_final).*sqrt(n_ph+1)).^2, 'go', 'LineWidth',1.5) % exact = cos(g*t*sqrt(n+1))^2 for excited state
%
hold off
box on
xlabel('Time (au)')
ylabel('Probaility') % 
set(gca,'FontSize',16)
%axis([0 n_max, 0 0.18])
box on


%%%
return
end
