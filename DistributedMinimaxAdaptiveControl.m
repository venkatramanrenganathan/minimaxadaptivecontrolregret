%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code simulates the control of a distributed implementation of the 
% minimax adaptive control algorithm.
%
% Copyrights Authors: 1) Venkatraman Renganathan - Lund University, Sweden.
%                     2) Anders Rantzer - Lund University, Sweden.
%
% Email: venkatraman.renganathan@control.lth.se
%
% Courtesy: 1.) Daniel Cedelberg - Linköping University, Sweden.
%           2.) Anders Hansson - Linköping University, Sweden.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a fresh start
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the network data using graph structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Form Incidence Matrix
I = [-1  0  0  0  0  0 
      1 -1 -1  0  0  0
      0  1  0 -1  0  0
      0  0  1  0 -1 -1
      0  0  0  1  0  0 
      0  0  0  0  1  0 
      0  0  0  0  0  1];
 
% Infer state and input dimensions
[n,m] = size(I);

% Get the concatenated dimension for forming Zbar = [-v x u] [.]'
concatDim = n+n+m;
  
% Form Dynamics b and B, C, D such that B/b = I
b = 0.1;
B = b*I;
C = eye(n);
D = zeros(n,m);

% Specify the adjacency matrix
adjMatrix = zeros(n,n);
adjMatrix(1,2) = 1;
adjMatrix(2,1) = 1;
adjMatrix(2,3) = 1;
adjMatrix(2,4) = 1;
adjMatrix(3,2) = 1;
adjMatrix(3,5) = 1;
adjMatrix(4,2) = 1;
adjMatrix(4,6) = 1;
adjMatrix(4,7) = 1;
adjMatrix(5,3) = 1;
adjMatrix(6,4) = 1;
adjMatrix(7,4) = 1;

% Get the degree vector
degreeVec = sum(adjMatrix, 1)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare the dynamics matrix A for each number of Models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the number of possible plant network dynamics
numModels = 2;

% Populate all the possible plant models
for dIter = 1:numModels      
    % Placeholders for a matrix for each component
    ai = zeros(n, 1);
    % Flag for stopping while loop
    dynPrepFlag = 1;
    % Counter for populating ai vector
    iter = 1;
    % Loop through until ai^2 + bb^{T} < ai is satisfied
    while dynPrepFlag
        % Generate a random ai
        ai(iter, 1) = rand;
        % Check condition
        if(ai(iter, 1)^2 - ai(iter, 1) + 2*b^2*degreeVec(iter, 1) < 0)        
            iter = iter + 1;    
        end
        % if counter exceeds limit, form the A matrix & stop while loop
        if(iter > n)
            A = diag(ai);
            dynPrepFlag = 0;
        end
    end
    % Prepare A and B matrices
    AMatrices{dIter} = A;
    BMatrices{dIter} = B; 
    % Check for correctness (condition on connectivity & speed of propagation)
    if(all(eig(AMatrices{dIter}) < 1))
        fprintf('# %d Dynamics agrees the condition A is Schur \n', dIter);
    end
    if(eig(AMatrices{dIter}^2 + BMatrices{dIter}*BMatrices{dIter}' - AMatrices{dIter}) < 1)
        fprintf('# %d Dynamics agrees the condition A^2 + BB^{T} < A \n', dIter);
    end    
    % Compute H_infinity control as u = Kx, with K = B'*(A-I)^{-1}
    Kinfinity{dIter} =  BMatrices{dIter}'*inv(AMatrices{dIter} - eye(n));
    if(any(eig(AMatrices{dIter} + BMatrices{dIter}*Kinfinity{dIter}) > 1))
        disp('Closed loop Eigenvalue > 1. So, there exists atleast an unstable state.');
    else
        disp('Closed loop Hinfinity controlled system is stable')
    end
end

disp('Finished Computing H_infinity Control Gains');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate the trajectory with Hinfinity & Minimax Adaptive Control (MAC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set time horizon
T = 50;

% Select disturbance. 0: none, 1: white noise, 2: sinusoidal, 3: constant step 
disturbFlag = 2;
if(disturbFlag == 0)
    w_advers = zeros(n, T+1);       % Zero Adversarial disturbance
elseif(disturbFlag == 1)
    w_advers = 0.1*randn(n, T+1);   % White Noise Adversarial disturbance
elseif(disturbFlag == 2)
    w_advers = zeros(n, T+1);       % Sinusoidal Adversarial disturbance
    for i = 1:n
        w_advers(i, :) = sin(1:T+1);
    end
elseif(disturbFlag == 3)
    w_advers = 0.5*ones(n, T+1);    % Step Adversarial disturbance
end


% Define Placeholders for states and inputs for MAC & H_infty policies
x_minmax = zeros(n, T+1);                 % MAC States  
u_minmax = zeros(m, T);                   % MAC Inputs
x_hinfty = zeros(n, T+1);                 % H_infty States 
u_hinfty = zeros(m, T);                   % H_infty Inputs
Z_models = zeros(numModels,1);            % Z for each numModel
Z_Matrix = zeros(T, numModels);           % Z vector of size T for each models
ZbarMatrix = zeros(concatDim, concatDim); % Zbar matrices 

% Generate & populate the same initial condition for both MAC & H_infty
x_0 = 15*ones(n,1);
x_minmax(:, 1) = x_0;
x_hinfty(:, 1) = x_0;

% Loop through and simulate until the end of time horizon
for t = 2:T+1  
    
    % Define true dynamics to be used from all the possible models in simulation
    modelNum = 2;
%     if(t<T/2)
%         modelNum = 2;
%     else
%         modelNum = 1;
%     end

    % Based on true modelNum specified, get the true dynamics & Gain matrices
    A_true = AMatrices{modelNum};
    B_true = BMatrices{modelNum};
    K_true = Kinfinity{modelNum};
    
    % Simulate the Hinfty system
    u_hinfty(:,t-1) = K_true*x_hinfty(:,t-1);  % H_infty full state feedback control
    x_hinfty(:,t) = A_true*x_hinfty(:,t-1) + B_true*u_hinfty(:,t-1) + w_advers(:,t-1);  % H_infty state update    
    
    % Simulate the MAC system    
    if(t > 2)
        wSums = zeros(numModels,1);
        for i = 1:numModels
            sumMatrix  = [eye(n) AMatrices{i} BMatrices{i}]';
            wSums(i,1) = trace(sumMatrix'*ZbarMatrix*sumMatrix);
        end
        [~, minIdx] = min(wSums);     % Get minIdx corresponding to least wSums 
    else
        [~, minIdx] = min(Z_models);  % Get minIdx corresponding to least Z  
    end
    u_minmax(:, t-1) = Kinfinity{minIdx}*x_minmax(:, t-1);  % MAC full state feedback control
    x_minmax(:, t) = A_true*x_minmax(:, t-1) + B_true*u_minmax(:, t-1) + w_advers(:,t-1); % MAC state update
    if(any(eig(A_true + B_true*Kinfinity{minIdx}) > 1))
        disp('Closed loop Eigenvalue > 1. So, there exists atleast an unstable state.');
    end
    
    % Form the Z for each models
    for i = 1:numModels
        Ai = AMatrices{i};
        Bi = BMatrices{i};
        Z_models(i) = Z_models(i) + norm(Ai*x_minmax(:, t-1) + Bi*u_minmax(:, t-1) - x_minmax(:, t))^2;
        Z_Matrix(t,i) = Z_models(i);        
    end
    
%     % Should we do this instead below by stacking node wise ??
%     ZbarUpdate = [];
%     for j = 1:n
%         ZbarUpdate = [ZbarUpdate; -x_minmax(j, t); x_minmax(j, t-1); u_minmax(j, t-1)];
%     end

    ZbarUpdate = [-x_minmax(:, t); x_minmax(:, t-1); u_minmax(:, t-1)];
    ZbarMatrix = ZbarMatrix + ZbarUpdate*ZbarUpdate';
    
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot minimax and Hinfinity trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Visualize the network by plotting it
figure(1);
plot(graph(adjMatrix),'layout','force');

figure(2);
tvec = 0:T; % Form the time vector for plotting
subplot(2,1,1);
plot(tvec, x_minmax', '*-');
xlabel('Time');
ylabel('States');
title('Minimax Adaptive Control States');
subplot(2,1,2); 
plot(tvec, (x_minmax - x_hinfty)', '*-');
xlabel('Time');
ylabel('States');
title('Difference of Minimax and H infinity States');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 40);