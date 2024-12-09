% This matlab code studies the two-level atom interacting with an external
% electromagnetic (EM) field  in terms of the Jaynes-Cummings (JC)
% Hamiltonian, defined as:
%
% H_JC = 0.5 * hbar * wa * sigma_z + hbar * wf * a^{+} * a + hbar * g * (sigma_minus * a^{+} + sigma_{+} * a); 
%
% An initial state of the field: a coherent state and atomic excited state.
% 
% W(t) = <sigma_z(t)> = exp(-|alpha|^{2}) * sum_{n=0} (|alpha|^(2*n)/n!) * cos(2*g*t*sqrt(n+1))
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
function [] = two_level_atom_with_EM_field_2
%
clear; clc; format short
%
n_max = 35.;
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
g = 0.50; % coupling constant
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
% coherent states
alpha = 2.750;
%
coh_st = 0.;
for n = 0:n_max-1
    coh_st = coh_st + exp(-0.5*abs(alpha)^2) * (abs(alpha)^(n)./sqrt(factorial(n))) * Vec(:,n+1);
end
%coh_st;
%
psi = kron(ex_st,coh_st); % initial state is coherent state 
dt = 0.10; % time step
ci = sqrt(-1.);
%
U_evol_op = expm(-ci * H_total_JC * dt/hbar); % unitary evolution operator
%
T_final = 60.;
T_time = 0:dt:T_final;
%
for t = 1:length(T_time)
    %
    W(t) = psi' * kron(sigma_z, If) * psi; % 
    %
    psi = U_evol_op * psi;
    psi = psi/norm(psi);
    %
end
%
%%%
dtt = 1.0;
T_time_ex = 0:dtt:T_final;
W_exact = zeros(1,length(T_time_ex));
%sm(1)
for tt = 1:length(T_time_ex)
    W_exact(tt) = exp(-abs(alpha)^2).* sum( ((alpha).^(2*(0:1:n_max-1))./factorial(0:1:n_max-1)).* cos(2.*g.*T_time_ex(tt).*sqrt((0:1:n_max-1)+1)) ) ;
end

%
figure(1)
hold on
plot(T_time, W, 'b-', 'LineWidth',1.5)
plot(T_time_ex, W_exact, 'ro', 'LineWidth',1.5)
%
xlabel('time (au)')
ylabel('W(t)')
hold off
box on
xlabel('Time (au)')
ylabel('W(t)') % 
set(gca,'FontSize',16)
%axis([0 n_max, 0 0.18])
box on




%%%
return
end
