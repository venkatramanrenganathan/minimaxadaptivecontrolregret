% Have a fresh start
clear all; close all; clc;

% Define system parameters - Difference of values should be greater than 2
a = [2 5]; 

% Define finite set of system matrices
B = [1 0]';
AMatrices = {};
BMatrices = {};
for i = 1:length(a)
    AMatrices{i} = [a(i) -1; 1 0];
    BMatrices{i} = B;
end
C = [1 0];
D = 0;

% Infer the system dimensions
[n,m] = size(B);
           
% Define design parameters.
startGamma = 500;
Q = eye(n);
R = eye(m);

% Get minimax adaptive control (MAC) policies
[status, gamma, PMatrices, KMatrices] = Approach3(AMatrices, BMatrices, Q, R, startGamma)
