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

% set properties for plotting
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 20;     % Set time horizon
tvec = 0:T; % A vector of times till time horizon 

% Define system parameters of double integrator
A = [2 -1  1; 
     1  0  0; 
     0  0  0];
B = [0  0  1]';


% Define the finite set of linear system models
AMatrices = {A, A};
BMatrices = {B, -B};
A1 = A;
A2 = A;
B1 = B;
B2 = -B;

% Infer the system dimensions
[n,m] = size(B);

% Get to know the number of linear system models available
numModels = length(BMatrices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Minimax Adaptive Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
% Define design parameters.
Q = eye(n);         % State penalty matrix 
R = eye(m);         % Input penalty matrix 
startGamma = 100;   % Starting gain from disturbance to error
computeFlag = 0;    % Flag to decide if we want to compute the gains or load the precomputed gains 

if(computeFlag == 1)
    % Get minimax adaptive control (MAC) policies solving LMIs
    disp('Started Computing Minimax Adaptive Control Gains');
    [MAC_status, MAC_gamma, MAC_PMatrices, MAC_KMatrices] = Approach2(AMatrices, BMatrices, Q, R, startGamma);    
    save('doubleIntegratorData.mat','MAC_gamma','MAC_PMatrices', 'MAC_KMatrices');
    disp('Finished Computing Minimax Control Gains');
else
    % Load the data fromt the file
    disp('Loading Precomputed Minimax Control Gains');
    load('doubleIntegratorData.mat');
    disp('Finished Loading Minimax Control Gains');
end

% Infer the square of gamma
gammaSqr = (MAC_gamma*MAC_gamma);

% Check the answers
correctnessFlag = 0;
if(all(eig(MAC_PMatrices{1,2} - MAC_PMatrices{1,1})) >= 0)
    disp('P_{12} > P_{11} Satisfied');
    correctnessFlag = 1;
end
if(all(eig(MAC_PMatrices{2,1} - MAC_PMatrices{2,2})) >= 0)
    disp('P_{21} > P_{22} Satisfied');
    correctnessFlag = 1;
end
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

% Flag to optimize Hinfty gamma. 0: Use MAC_Gamma, 1: Get Bisection gamma
gammaOptFlag = 0; 

% Get H-infinity Optimal Control & adversarial disturbance via Dynamic
% Game Approach - USE MINIMAX GAMMA
for i = 1:numModels
    [K_Gains{i,1}, F_Gains{i,1}, Gammas(i,1)] = HinfDynamicGame(AMatrices{i}, BMatrices{i}, Q, R, MAC_gamma, gammaOptFlag);
end
disp('Finished Computing H_infinity Control Gains');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the system to compute the regret.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select the disturbance type: 1: adverse, 2: confuse, 3: sine
disturbanceSelect = 1;   

% Define Placeholders for states and inputs for MAC & H_inf policies
x_minmax  = zeros(n, T+1);             % States of MAC
u_minmax  = zeros(m, T);               % MAC inputs
x_hinfty  = zeros(n, T+1, numModels);  % States of H_infty  
u_hinfty  = zeros(m, T, numModels);    % H_infty inputs
w_advers  = zeros(n, T, numModels);    % Adversarial disturbance

% Populate the initial condition
x_0 = [-5 1 3]';
for i = 1:numModels
    x_minmax(:, 1) = x_0;
    x_hinfty(:, 1, i) = x_0;
end

% Frequency response parameters
Ts = 0.01;
z = tf('z', Ts);

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
        % Generate the disturbance according to the selection made
        if(disturbanceSelect == 1)
            % Generate the H_infty worst case adversarial disturbance
            w_advers(:,t-1,i) = Fi*x_hinfty(:,t-1,i); 
        elseif(disturbanceSelect == 2)
            % Generate confusing adversarial disturbance
            modelsCombination = ones(numModels,1);
            if(mod(t,2)==0)
                modelsCombination(i) = -1;             
            end
            % modelsCombination(i) = -1; % Subtract contribution of i^th model            
            for j = 1:numModels
                w_advers(:,t-1,i) = w_advers(:,t-1,i) + modelsCombination(j)*(AMatrices{j}*x_hinfty(:,t-1,i) + BMatrices{j}*u_hinfty(:,t-1,i));
            end
        elseif(disturbanceSelect == 3)
            % Find the transfer function from disturbance to [x,u]
            T_d_to_z_i = [eye(n); Ki]*inv(z*eye(n) - (Ai + Bi*Ki));
            % Find the frequency where worst gain occurs
            [gpeaki, fpeaki] = getPeakGain(T_d_to_z_i);
            % Evaluate the transfer function at the frequency
            sysi_at_fmax = evalfr(T_d_to_z_i, exp(j*fpeaki*Ts));
            % Perform a Singular value decomposition
            [Ui, Si, Vi] = svd(sysi_at_fmax);
            % Get the angle of Vi complex vector
            viAngle = angle(Vi);
            % Generate sinusoidal adversarial disturbance at that frequency
            for j = 1:n
                w_advers(j,t-1,i) = sin(fpeaki*(t-1) + viAngle(j))*real(Vi(j,1));
            end
        end    
        % Update the state with the Hinfty control & generated disturbance
        x_hinfty(:,t,i) = Ai*x_hinfty(:,t-1,i) + Bi*u_hinfty(:,t-1,i) + w_advers(:,t-1,i);        
    end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the trajectory with the Minimax Adaptive control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define true dynamics to be used from all the possible models in simulation
modelNum = 2;

% Define the sign of K matrix based on the model selected
if(modelNum == 1)
    sgn = 1;
elseif(modelNum == 2)
    sgn = -1;
end    

% Adversarial disturbance
w_minimax = zeros(n, T);  
  
% Rename model matrices for notational brevity
Fi = F_Gains{modelNum,1};
zplus = 0;
zminus = 0;

% Loop through the entire time horizon
for t = 2:T+1         

    % Compute the minimax full state feedback control   
    if(zplus <= zminus)
        % (A, B) fits disturbance model well
        modelSelected = 1;
        Ki = MAC_KMatrices{1};
        disp('Controller selects model 1 as the best fit for disturbance trajectory');
        u_minmax(:,t-1) = -MAC_KMatrices{1}*x_minmax(:,t-1); 
    else
        % (A, -B) fits disturbance model well
        modelSelected = 2;
        Ki = MAC_KMatrices{2};
        disp('Controller selects model 2 as the best fit for disturbance trajectory');
        u_minmax(:,t-1) = -MAC_KMatrices{2}*x_minmax(:,t-1); 
    end     
    
    % Generate the disturbance according to the selection made
    if(disturbanceSelect == 1)
        % Use the worst case disturbance policy instead
        w_minimax(:,t-1) = Fi*x_minmax(:,t-1); 
    elseif(disturbanceSelect == 2)
        % Generate Confusing disturbance
        modelsCombination = ones(numModels,1);
        if(mod(t,2)==0)
            modelsCombination(modelNum) = -1; % Subtract the contribution of original        
        end
        for i = 1:numModels
            w_minimax(:,t-1) = w_minimax(:,t-1) + modelsCombination(i)*(AMatrices{i}*x_minmax(:,t-1) + BMatrices{i}*u_minmax(:,t-1));
        end
    elseif(disturbanceSelect == 3)
        % Find the transfer function from disturbance to [x,u]
        T_d_to_z_minimax = [eye(n); Ki]*inv(z*eye(n) - (A + B*Ki));
        % Find the frequency where worst gain occurs
        [gpeak_minmax, fpeak_minmax] = getPeakGain(T_d_to_z_minimax);
        % Evaluate the transfer function at the frequency
        sys_minmax_at_fmax = evalfr(T_d_to_z_minimax, exp(j*fpeak_minmax*Ts));
        % Perform a Singular value decomposition
        [U_minmax, S_minmax, V_minmax] = svd(sys_minmax_at_fmax);
        % Get the angle of V_minmax complex vector
        vminmax_Angle = angle(V_minmax);
        % Generate sinusoidal adversarial disturbance at that frequency
        for j = 1:n
            w_minimax(j,t-1) = sin(fpeak_minmax*(t-1) + vminmax_Angle(j))*real(V_minmax(j,1));
        end
    end    

    % Perform the state update with the minimax control & Hinfty noise
    x_minmax(:,t) = A*x_minmax(:,t-1) + sgn*B*u_minmax(:,t-1) + w_minimax(:,t-1);  

    % Update the Z variable 
    zplus = zplus + gammaSqr*norm(A*x_minmax(:,t-1) + B*u_minmax(:,t-1) - x_minmax(:,t))^2;
    zminus = zminus + gammaSqr*norm(A*x_minmax(:,t-1) - B*u_minmax(:,t-1) - x_minmax(:,t))^2;

end    

%% Plot the minimax trajectories vs time
figNum = 1;
figure(figNum);
plot(tvec, x_minmax(:,:)');
hold on;
xlabel('Time');
ylabel('$x^{\dagger}$', 'interpreter', 'latex');
hold off;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 40);

%% Plot the Hinfty trajectories vs time
figNum = figNum + 1;
figure(figNum);
plot(tvec, x_hinfty(:,:,2)');
hold on;
xlabel('Time');
ylabel('$x^{\star}_{2}$', 'interpreter', 'latex');
hold off;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 40);

%% Plot the Difference trajectories vs time
figNum = figNum + 1;
figure(figNum);
plot(tvec, (x_minmax(:,:) - x_hinfty(:,:,2))');
hold on;
xlabel('Time');
ylabel('$x^{\dagger} - x^{\star}_{2}$', 'interpreter', 'latex');
hold off;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 40);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the Model-based Regret & the Total Regret
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Record regret values at each time and Initialize to zero for each model
modelRegrets = zeros(numModels, T+1);    

% Loop through all finite set of linear system models
for i = 1:numModels   
    % Loop through the entire time horizon
    for t = 1:T+1
        minmaxStateSum = 0;
        minmaxInputSum = 0;
        hinftyStateSum = 0;
        hinftyInputSum = 0;
        disturbanceSum = 0;
        % compute u'Ru and add it to control sum
        if(t < T+1)
            minmaxInputSum = u_minmax(:,t)'*R*u_minmax(:,t);
            hinftyInputSum = u_hinfty(:,t,i)'*R*u_hinfty(:,t,i);
            disturbanceSum = w_minimax(:,t)'*w_minimax(:,t);
        end
        % compute x'Qx and add it to state sum
        minmaxStateSum = x_minmax(:,t)'*Q*x_minmax(:,t);
        hinftyStateSum = x_hinfty(:,t,i)'*Q*x_hinfty(:,t,i);  
        
        % Record the regret incurred at time t
        modelRegrets(i, t) = (minmaxStateSum + minmaxInputSum) - (hinftyStateSum + hinftyInputSum) - gammaSqr*disturbanceSum;         
    end
end

% Compute modelbased regret as cumulative sum upto T+1
modelbasedRegret = cumsum(modelRegrets, 2);

% Compute totalRegret as max of modelbasedRegret
totalRegret = max(modelbasedRegret(:, end));
fprintf('Total regret: %.3f \n', totalRegret);

%% Plot the regret vs time
Tvec = 0:T;
figNum = figNum+1;
figure(figNum);
plot(Tvec, modelbasedRegret(2, :), '-ob');
xlabel('Time');
ylabel('$\mathcal{R}(\pi^{\dagger}, \pi^{\star}_{2})$');
hold off;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 40);

%% Plot the Hinfty control inputs vs time
Tvec = 0:T-1;
figNum = figNum + 1;
figure(figNum);
for i = 1:numModels
    if(i==1)
        plot(Tvec, u_hinfty(:,:,i), '-ob');
    else
        plot(Tvec, u_hinfty(:,:,i), '-or');
    end
    xlabel('Time');
    ylabel('$u^{\star}_{1}, u^{\star}_{2}$');
    hold on;
    a = findobj(gcf, 'type', 'axes');
    h = findobj(gcf, 'type', 'line');
    set(h, 'linewidth', 4);
    set(a, 'linewidth', 4);
    set(a, 'FontSize', 40);
end

%% Plot the Minimax Adaptive control inputs vs time
Tvec = 0:T-1;
figNum = figNum + 1;
figure(figNum);
plot(Tvec, u_minmax(:,:), '-or');
xlabel('Time');
ylabel('$u^{\dagger}$');
hold on;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 40);

%% Plot the Hinfinity disturbance inputs vs time
Tvec = 0:T-1;
figNum = figNum + 1;
figure(figNum);
for i = 1:numModels
    if(i==1)
        plot(Tvec, w_advers(:,:,i)');
    else
        plot(Tvec, w_advers(:,:,i)');
    end
    xlabel('Time');
    ylabel('$w^{\star}_{1}, w^{\star}_{2}$');
    hold on;
    a = findobj(gcf, 'type', 'axes');
    h = findobj(gcf, 'type', 'line');
    set(h, 'linewidth', 4);
    set(a, 'linewidth', 4);
    set(a, 'FontSize', 40);
end

%% Plot the Minimax disturbance inputs vs time
Tvec = 0:T-1;
figNum = figNum + 1;
figure(figNum);
plot(Tvec, w_minimax(:,:)');
xlabel('Time');
ylabel('$w^{\dagger}$');
hold on;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 40);

