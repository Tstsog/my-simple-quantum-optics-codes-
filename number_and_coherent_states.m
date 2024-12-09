% This matlab code defines the number operator in terms of the
% annihiliation (a) and creation (a^{dagger}) operator, 
% whose eigenvalue and eigenvector correspond to number of photon in the quantized electromagnetic (EM) field and 
% number state (function) with same photon number, respectively; a coherent states representing superposition of 
% the number states. Also computes a probability of occupying n photon in a
% coherent state. 
%
% Coherent states: |alpha > = exp(-|alpha|^{2}) * sum_{n=0}(alpha^{2}/sqrt(n!))*|n >.
%
% Probability, P_{n}: P_{n} = |<n|alpha>|^{2} = exp(-|alpha|^{2}) *|alpha|^{2*n}/n! .
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
function [] = number_and_coherent_states
%
clear; clc; format short
%
n_max = 25.; % dimension of number operator N
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coherent state (coh_st)
%
alpha = 2.5;
%
coh_st = 0.;
for n = 0:n_max-1
    coh_st = coh_st + exp(-0.5*abs(alpha)^2) * (alpha^(n)./sqrt(factorial(n))) * Vec(:,n+1);
end
%coh_st;
%
norm(coh_st); % = 1.0000
%%%

for n = 0:n_max-1
    Prob_n(n+1) = abs(Vec(:,n+1)' * coh_st)^2; % P_{n} = |<n|alpha>|^2
end
%Prob_n;

% an exact formula for P_{n}
Pn_exact = exp(-abs(alpha)^2).*abs(alpha).^(2.*(0:1:n_max-1))./factorial((0:1:n_max-1)); % P_n = e^(|alpha|^2)*|alpha|^{2n}/n!

%
figure(1)
hold on
plot((0:1:n_max-1), Prob_n, 'b-', 'LineWidth',1.5)
plot((0:1:n_max-1), Pn_exact, 'ro', 'LineWidth',1.5)
hold off
xlabel('number of photon')
ylabel('Probaility') % of number state |n> in coherent state
set(gca,'FontSize',16)
axis([0 n_max, 0 0.18])
box on


%%%
return
end
