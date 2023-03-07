%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code simulates the control of a double integrator model using the 
% minimax adaptive control algorithm.
%
% Copyrights Authors: 1) Venkatraman Renganathan - Lund University, Sweden.
%                     2) Anders Rantzer - Lund University, Sweden.
%                     3) Andrea Ianelli - University of Stuttgart
%
% Email: venkatraman.renganathan@control.lth.se
%
% Courtesy: 1.) Daniel Cedelberg - Linköping University, Sweden.
%           2.) Anders Hansson - Linköping University, Sweden.
%           3.) Olle Kjellkvist - Lund University, Sweden.
%
%
% RUNNING WITH MOSEK - if it causes trouble
% 1. Go to terminal in mac
% 2. type the following and press enter: xattr -dr com.apple.quarantine /Users/venkat/Documents/MATLAB/mosek/10.0/tools/platform/osx64x86/bin
% 3. Include Mosek into path
% 4. In matlab command window type the following and press enter: setenv('PATH', [getenv('PATH') ';/Users/venkat/Documents/MATLAB/mosek/10.0/tools/platform/osx64x86/bin']);
% 5. In matlab command window type the following and press enter: mosekdiag
% 6. MATLAB will now give the message that "mosekopt works OK. You can use MOSEK in MATLAB."
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a fresh start
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set time horizon
T = 20;
tvec = 0:T; 

% Define system parameters of double integrator
A = [2 -1  1; 
     1  0  0; 
     0  0  0];
B = [0  0  1]';
C = [1  0  0];
D = 0;

% Define system parameters of vinnicombe example
% A = 1;
% B = 1;
% C = 1;
% D = 0;

% % F16 Aircraft
% A = [1 0.1025 0.2080 -0.0502 -0.0057
%      0 1.1175 4.1534 -0.8000 -0.1010
%      0 0.0955 1.0722 -0.0541 -0.0153
%      0 0      0       0.1353  0
%      0 0      0       0       0.1353];
% B = [-0.0377  -0.0040
%      -1.0042  -0.1131
%      -0.0453  -0.0175
%      0.8647    0
%      0         0.8647];
% C = [1 0 0 0 0
%      0 1 0 0 0];
% D = zeros(size(C, 1));


% Define the finite set of linear system models
AMatrices = {A, A};
BMatrices = {B, -B};

% Infer the system dimensions
[n,m] = size(B);

% Get to know the number of linear system models available
numModels = length(BMatrices);

% Form the state-space with difference linear models
Plants = cell(numModels, 1);
for i = 1:numModels
    Plants{i, 1} = ss(AMatrices{i}, BMatrices{i}, C, D);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Minimax Adaptive Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           


% Define design parameters.
Q = eye(n);         % State penalty matrix 
R = eye(m);         % Input penalty matrix 
startGamma = 100;   % Starting gain from disturbance to error

computeFlag = 0;

if(computeFlag == 1)
    % Get minimax adaptive control (MAC) policies solving LMIs
    disp('Started Computing Minimax Adaptive Control Gains');
    [MAC_status, MAC_gamma, MAC_PMatrices, MAC_KMatrices] = Approach3(AMatrices, BMatrices, Q, R, startGamma);    
    save('doubleIntegratorData.mat','MAC_gamma','MAC_PMatrices', 'MAC_KMatrices');
    disp('Finished Computing Minimax Control Gains');
else
    % Load the data fromt the file
    load('doubleIntegratorData.mat');
    disp('Finished Loading Minimax Control Gains');
end

% Infer the square of gamma
gammaSqr = (MAC_gamma*MAC_gamma);

% Check the answers
P11Logic = prod(MAC_PMatrices{1,1}, 'all'); 
P12Logic = prod(MAC_PMatrices{1,2}, 'all'); 
P21Logic = prod(MAC_PMatrices{2,1}, 'all'); 
P22Logic = prod(MAC_PMatrices{2,2}, 'all'); 

correctnessFlag = 0;
if(all(eig(MAC_PMatrices{1,2} - MAC_PMatrices{1,1})) >= 0)
    disp('P_{12} > P_{11} Satisfied');
    correctnessFlag = 1;
end
correctnessFlag = 0;
if(all(eig(MAC_PMatrices{2,1} - MAC_PMatrices{2,2})) >= 0)
    disp('P_{21} > P_{22} Satisfied');
    correctnessFlag = 1;
end
correctnessFlag = 0;
if(all(eig(gammaSqr*eye(n) - MAC_PMatrices{2,1})) >= 0)
    disp('P_{21} < gammaSqr*I Satisfied');
    correctnessFlag = 1;
end
if(correctnessFlag == 1)
    disp('Computed Minimax Adaptive Control Gains are Correct');
else
    disp('Computed Minimax Adaptive Control Gains has Problems');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute hinf controller and disturbance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Started Computing H_infinity Control Gains');

% Initialize the control gains, disturbance gains and worst case gamma 
K_Gains = cell(numModels, 1);
F_Gains = cell(numModels, 1);
Gammas  = zeros(numModels, 1);

% Get H-infinity Optimal Control & adversarial disturbance via Dynamic Games
for i = 1:numModels
    [K_Gains{i,1}, F_Gains{i,1}, Gammas(i,1)] = HinfDynamicGame(Plants{i,1}, Q, R);
end
disp('Finished Computing H_infinity Control Gains');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the system to compute the regret.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Placeholders for states and inputs for MAC & H_inf policies
x_minmax = zeros(n, T+1);  % States of MAC
u_minmax = zeros(m, T);    % MAC inputs
x_hinfty = zeros(n, T+1, numModels);  % States of H_infty  
u_hinfty = zeros(m, T, numModels);    % H_infty inputs
w_advers = zeros(n, T+1, numModels);  % Adversarial disturbance

% Populate the initial condition
% x_0 = -5;
x_0 = [-5 1 10]';
% x_0 = [-5 1 10 3 -2]';
for i = 1:numModels
    x_minmax(:, 1) = x_0;
    x_hinfty(:, 1, i) = x_0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Simulate the trajectory with the Hinfinity control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through all finite set of linear system models
for i = 1:numModels    
    % Rename model matrices for notational brevity
    Ai = AMatrices{i};
    Bi = BMatrices{i};
    Ki = K_Gains{i,1};
    Fi = F_Gains{i,1};    
    % Loop through the entire time horizon
    for t = 2:T+1        
        % Compute the H_infty full state feedback control
        u_hinfty(:,t-1,i) = Ki*x_hinfty(:,t-1,i);        
        % Compute the H_infty worst case adversarial disturbance
        w_advers(:,t-1,i) = Fi*x_hinfty(:,t-1,i); % sin(t-1)*ones(n, 1); % Fi*x_hinfty(:,t-1,i);        
        % Perform the state update with the Hinfty control & Hinfty noise
        x_hinfty(:,t,i) = Ai*x_hinfty(:,t-1,i) + Bi*u_hinfty(:,t-1,i) + w_advers(:,t-1,i);        
    end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the trajectory with the Minimax Adaptive control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define true dynamics to be used from all the possible models in simulation
modelNum = 2;

% Adversarial disturbance
w_advers = zeros(n, T+1);  
  
% Rename model matrices for notational brevity
Ai = AMatrices{modelNum};
Bi = BMatrices{modelNum};
Fi = F_Gains{modelNum,1};
zplus = 0;
zminus = 0;
zpp = zeros(T,1);
zmm = zeros(T,1);

% Loop through the entire time horizon
for t = 2:T+1         

    % Compute the minimax full state feedback control   
    if(zplus<=zminus)
        % (A, B) fits disturbance model well
        disp('Controller selects model 1 as the best fit for disturbance trajectory');
        u_minmax(:,t-1) = -MAC_KMatrices{1}*x_minmax(:,t-1); 
    else
        % (A, -B) fits disturbance model well
        disp('Controller selects model 2 as the best fit for disturbance trajectory');
        u_minmax(:,t-1) = -MAC_KMatrices{2}*x_minmax(:,t-1); 
    end     
    
    % Compute the H_infty worst case adversarial disturbance
%     w_advers(:,t-1) = Fi*x_minmax(:,t-1); % sin(t-1)*ones(n, 1); 
    
    % Construct Confusing disturbance
    modelsCombination = rand(1, numModels);
    if(mod(t,2)==0)
        modelsCombination(modelNum) = -1;
    else
        modelsCombination(modelNum-1) = -1;
    end
    w_advers_sum = 0;
    for i = 1:numModels
        w_advers_sum = w_advers_sum + modelsCombination(i)*(AMatrices{i}*x_minmax(:,t-1) + BMatrices{i}*u_minmax(:,t-1));
    end
    
    % Set the flag if you want to add randomness to the disturbance
    addrandomFlag = 0;
    if (addrandomFlag == 1)    
        addRandomness = 0.1*rand(n, 1); % 
        w_advers(:,t-1) = w_advers_sum + addRandomness;
    else        
        w_advers(:,t-1) = w_advers_sum;
    end

    % Perform the state update with the minimax control & Hinfty noise
    x_minmax(:,t) = Ai*x_minmax(:,t-1) + Bi*u_minmax(:,t-1) + w_advers(:,t-1);  

    % Update the Z variable 
    zplus = zplus + gammaSqr*norm(A*x_minmax(:,t-1) + B*u_minmax(:,t-1) - x_minmax(:,t))^2;
    zminus = zminus + gammaSqr*norm(A*x_minmax(:,t-1) - B*u_minmax(:,t-1) - x_minmax(:,t))^2;

    % Store history of z values
    zpp(t-1) = zplus;
    zmm(t-1) = zminus;
end    





% % Loop through all finite set of linear system models
% for i = 1:numModels    
%     % Rename model matrices for notational brevity
%     Ai = AMatrices{i};
%     Bi = BMatrices{i};
%     Fi = F_Gains{i,1};
%     Ki = MAC_KMatrices{i};
%     Z  = zeros(2*n+m, 2*n+m);
%     % Loop through the entire time horizon
%     for t = 2:T+1       
%         % Compute the H_infty worst case adversarial disturbance
%         w_advers(:,t-1,i) = Fi*x_minmax(:,t-1,i);
%         % Variable storing state update difference for all j models
%         jModeldiff = zeros(numModels, 1);
%         for j = 1:numModels                
%             % Rename model matrices for notational brevity
%             Aj = AMatrices{j};
%             Bj = BMatrices{j};
%             Kj = MAC_KMatrices{j};
%              
%         end
%         % Get controlIndex based on a minimization
%         [~, controlIndex] = min(jModeldiff);
%         % Compute the minimax full state feedback control   
%         if(t > 2)
%             u_minmax(:,t-1,i) = -MAC_KMatrices{controlIndex}*x_minmax(:,t-1,i);                
%         else
%             u_minmax(:,t-1,i) = -Ki*x_minmax(:,t-1,i);        
%         end        
%         % Perform the state update with the minimax control & Hinfty noise
%         x_minmax(:,t,i) = Ai*x_minmax(:,t-1,i) + Bi*u_minmax(:,t-1,i) + w_advers(:,t-1,i);  
%         % Update the Z variable as Z = Z + [-v x u][-v x u]^T
%         Z_update = [-x_minmax(:,t,i); x_minmax(:,t-1,i); u_minmax(:,t-1,i)];
%         Z = Z + Z_update*Z_update';
%     end    
% end


%% Plot minimax and Hinfinity trajectories
for i = 1:numModels
    figure(i);
    plot(tvec, (x_minmax(:,:)-x_hinfty(:,:,i))');
    hold on;
    xlabel('Time');
    ylabel('States');
    title('Difference States');
    hold off;
    a = findobj(gcf, 'type', 'axes');
    h = findobj(gcf, 'type', 'line');
    set(h, 'linewidth', 4);
    set(a, 'linewidth', 4);
    set(a, 'FontSize', 40);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the Model-based Regret & the Total Regret
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiate model based regret to zero for each model
modelbasedRegret = zeros(numModels, 1);
modelRegrets = zeros(T+1,numModels);

% Loop through all finite set of linear system models
for i = 1:numModels   
    minmaxStateSum = 0;
    minmaxInputSum = 0;
    hinftyStateSum = 0;
    hinftyInputSum = 0;
    % Loop through the entire time horizon
    for t = 1:T+1
        % compute u'Ru and add it to control sum
        if(t < T+1)
            minmaxInputSum = minmaxInputSum + u_minmax(:,t)'*R*u_minmax(:,t);
            hinftyInputSum = hinftyInputSum + u_hinfty(:,t,i)'*R*u_hinfty(:,t,i);
        end
        % compute x'Qx and add it to state sum
        minmaxStateSum = minmaxStateSum + x_minmax(:,t)'*Q*x_minmax(:,t);
        hinftyStateSum = hinftyStateSum + x_hinfty(:,t,i)'*Q*x_hinfty(:,t,i);  
        
        % Record the regret incurred upto that time
        modelRegrets(t, i) = (minmaxStateSum + minmaxInputSum) - (hinftyStateSum + hinftyInputSum); 
    end
    minmaxStateInputSum = minmaxStateSum + minmaxInputSum;
    hinftyStateInputSum = hinftyStateSum + hinftyInputSum;
    modelbasedRegret(i,1) = minmaxStateInputSum - hinftyStateInputSum;
end

% Compute totalRegret as max of modelbasedRegret
totalRegret = max(modelbasedRegret);
fprintf('Total regret: %.3f \n', totalRegret);

%% Plot the regret vs time
Tvec = 0:T;
figure();
plot(Tvec, modelRegrets(:, 1), '-ob');
xlabel('Time');
ylabel('Total Regret');
hold off;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 40);
