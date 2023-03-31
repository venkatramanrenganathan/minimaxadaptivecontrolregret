%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code simulates the control of a random system model using the 
% minimax adaptive control algorithm.
%
% Copyrights Authors: 1) Venkatraman Renganathan - Lund University, Sweden.
%                     2) Anders Rantzer - Lund University, Sweden.
%                     3) Andrea Iannelli - University of Stuttgart, Germany
%
% Email: venkatraman.renganathan@control.lth.se
%
% Courtesy: 1.) Daniel Cedelberg - Linköping University, Sweden.
%           2.) Anders Hansson - Linköping University, Sweden.
%           3.) Olle Kjellkvist - Lund University, Sweden.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a fresh start
clear all; close all; clc;

% set properties for plotting
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Get current working folder info
currentFolder = pwd;
currentFolder = strcat(currentFolder, '/CDC_2023_Figures/');

% Set figure number
figNum = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 100;     % Set time horizon
tvec = 0:T; % A vector of times till time horizon 

% Set the number of models needed
numModels = 4;

% Set the system dimension
n = 3; % Number of states
m = 1; % Number of controls

% Define the finite set of linear system models
AMatrices = {};
BMatrices = {};

% Set flag to prepare or load the precomputed random system matrices
computeFlag = 0;

if(computeFlag == 1)
    disp('Started preparing random system matrices');
    % Prepare random systems of numModels count
    for i = 1:numModels
        % Set a flag to denote if the random system is controllable
        notOkayFlag = 1;
        while(notOkayFlag == 1)
            % generate random system
            A = 2*rand(n, n);
            % A = (A + A')/2;
            B = 2*rand(n, 1);
            % find the rank of controllability matrix
            if(rank(ctrb(A,B)) == n && all(eig(A) < 1))
                notOkayFlag = 0;
                AMatrices{i} = A;
                BMatrices{i} = B;
            end            
        end
    end
    save('RandomSystemData.mat','AMatrices','BMatrices');
    disp('Finished preparing random system matrices');
else
    disp('Started loading precomputed random system matrices');
    load('RandomSystemData.mat');
    disp('Finished loading precomputed random system matrices');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Minimax Adaptive Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
% Define design parameters.
Q = eye(n);         % State penalty matrix 
R = eye(m);         % Input penalty matrix 
startGamma = 100;   % Starting gain from disturbance to error

if(computeFlag == 1)
    % Get minimax adaptive control (MAC) policies solving LMIs
    disp('Started Computing Minimax Adaptive Control Gains');
    [MAC_status, MAC_gamma, MAC_PMatrices, MAC_KMatrices] = Approach3(AMatrices, BMatrices, Q, R, startGamma);    
    save('RandomSystemMinimaxData.mat','MAC_gamma','MAC_PMatrices', 'MAC_KMatrices');
    disp('Finished Computing Minimax Control Gains');
else
    % Load the data fromt the file
    disp('Loading Precomputed Minimax Control Gains');
    load('RandomSystemMinimaxData.mat');
    disp('Finished Loading Minimax Control Gains');
end

% Infer the square of gamma
gammaSqr = (MAC_gamma*MAC_gamma);

% Check the answers - Need to do it for all nxn matrices
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
if(all(eig(gammaSqr*eye(n) - MAC_PMatrices{3,1})) >= 0)
    disp('P_{31} < gammaSqr*I Satisfied');
    correctnessFlag = 1;
end
if(all(eig(gammaSqr*eye(n) - MAC_PMatrices{2,3})) >= 0)
    disp('P_{23} < gammaSqr*I Satisfied');
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

% Initialize the control gains, disturbance gains and worst case gamma 
K_Gains = cell(numModels, 1);
F_Gains = cell(numModels, 1);
Gammas  = zeros(numModels, 1);

if(computeFlag == 1)
    disp('Started Computing H_infinity Control Gains');
    % Flag to optimize Hinfty gamma. 0: Use MAC_Gamma, 1: Get Bisection gamma
    gammaOptFlag = 0; 

    % Get H-infinity Control & adversarial disturbance via Dynamic Games Approach - USE MINIMAX GAMMA
    for i = 1:numModels
        [K_Gains{i,1}, F_Gains{i,1}, Gammas(i,1)] = HinfDynamicGame(AMatrices{i}, BMatrices{i}, Q, R, MAC_gamma, gammaOptFlag);
    end
    save('RandomSystemHinftyData.mat','Gammas','K_Gains', 'F_Gains');
    disp('Finished Computing H_infinity Control Gains');
else
    % Load the data fromt the file
    disp('Loading Precomputed H_infinity Control Gains');
    load('RandomSystemHinftyData.mat');
    disp('Finished Loading H_infinity Control Gains');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the system to compute the regret.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select the disturbance type: 1: adverse, 2: confuse, 3: sinusoidal
disturbanceSelect = 2;   

% Define true dynamics to be used from all the possible models in simulation
modelNum = 2; 

% Set the index of dynamics that confusing w should trick policy to choose
cheatModelNum = 3;

% Define Placeholders for states and inputs for MAC & H_inf policies
x_minmax  = zeros(n, T+1);  % States of MAC
u_minmax  = zeros(m, T);    % MAC inputs
x_hinfty  = zeros(n, T+1);  % States of H_infty  
u_hinfty  = zeros(m, T);    % H_infty inputs
w_advers  = zeros(n, T);    % Adversarial disturbance

% Populate the initial condition
x_0 = [-5 1 3]'; % will also work for rand(n,1); 
x_minmax(:, 1) = x_0;
x_hinfty(:, 1) = x_0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the trajectory with the Minimax Adaptive control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the true dynamics
ATrue = AMatrices{modelNum};
BTrue = BMatrices{modelNum};

% Adversarial disturbance
w_minimax = zeros(n, T);  

% Set placeholder for zData
zData = zeros(numModels, 1);

% Loop through the entire time horizon
for t = 2:T+1         

    % Get model that best fits disturbance trajectory in least squares sense
    [zMin, zID] = min(zData);
    fprintf('Controller selects model %d as the best fit for disturbance trajectory \n', zID);
    
    % Compute the minimax full state feedback control   
    u_minmax(:,t-1) = -MAC_KMatrices{zID}*x_minmax(:,t-1);
    
    % Generate the disturbance according to the selection made
    if(disturbanceSelect == 1 || disturbanceSelect == 3)        
        % Use same worst-case adversarial/sinusoidal hinfinity disturbance
        w_minimax(:,t-1) = w_advers(:,t-1); 
    elseif(disturbanceSelect == 2)
        % Generate Confusing disturbance
%         modelsCombination = 0.01*rand(numModels,1);        
        modelsCombination = zeros(numModels,1);        
        % Subtract the contribution of original 
        modelsCombination(modelNum) = -1;                
        % Choose the i = cheatModelNum \neq nodelNum which you need policy to choose
        modelsCombination(cheatModelNum) = 1;
        for i = 1:numModels
%             if(mod(t,2) == 0 || mod(t,3) == 0)
%                 modelsCombination(i) = 1;        
%             end
            w_minimax(:,t-1) = w_minimax(:,t-1) + modelsCombination(i)*(AMatrices{i}*x_minmax(:,t-1) + BMatrices{i}*u_minmax(:,t-1));
        end
    end    
    
    % Perform the state update with the minimax control & Hinfty noise
    x_minmax(:,t) = ATrue*x_minmax(:,t-1) + BTrue*u_minmax(:,t-1) + w_minimax(:,t-1);  

    % Update the Z variable     
    for i = 1:numModels
        zData(i, 1) = zData(i, 1) + gammaSqr*norm(AMatrices{i}*x_minmax(:,t-1) + BMatrices{i}*u_minmax(:,t-1) - x_minmax(:,t))^2;
    end

end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Simulate the trajectory with the Hinfinity control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get Hinfty control for i = modelNum
i = modelNum;
% Rename model matrices for notational brevity
Ai = AMatrices{i};
Bi = BMatrices{i};
Ki = K_Gains{i,1};
Fi = F_Gains{i,1};    

% Frequency response parameters
Ts = 1;
z = tf('z', Ts);
% Find the transfer function from disturbance to [x,u]
T_d_to_z_i = [eye(n); Ki]*inv(z*eye(n) - (Ai + Bi*Ki));
% Find & record the frequency where worst gain occurs
[peakGain, peakFreq] = getPeakGain(T_d_to_z_i);
% [peakGain, peakFreq] = hinfnorm(T_d_to_z_i);
% Evaluate T_d_to_z_i at z = e^{jwn}
T_d_to_z_i_eval = [eye(n); Ki]*inv(exp(j*peakFreq*Ts)*eye(n) - (Ai + Bi*Ki));
% Evaluate the transfer function at the peak frequency
% sysi_at_fmax = evalfr(T_d_to_z_i, exp(j*peakFreq*Ts));
% sysi_at_fmax = evalfr(T_d_to_z_i, peakFreq);
% Perform a Singular value decomposition
[Ui, Si, Vi] = svd(T_d_to_z_i_eval);
% Get the angle of Vi complex vector
viAngle = angle(Vi);
% Get the magnitude of Vi complex vector
viMagnitude = abs(Vi);

% Loop through the entire time horizon
for t = 2:T+1  
    % Compute the H_infty full state feedback control
    u_hinfty(:,t-1) = Ki*x_hinfty(:,t-1); 
    % Generate the disturbance according to the selection made
    if(disturbanceSelect == 1)
        % Generate the H_infty worst case adversarial disturbance
        w_advers(:,t-1) = Fi*x_hinfty(:,t-1); 
    elseif(disturbanceSelect == 2)
        % Use the same confusing disturbance from Minimax control
        w_advers(:,t-1) = w_minimax(:,t-1);
    elseif(disturbanceSelect == 3)            
        % Generate sinusoidal adversarial disturbance at peak frequency 
        % Use Vi(:, 1) for worst case input direction
        for j = 1:n
            w_advers(j,t-1) = viMagnitude(j,1)*cos(peakFreq*(t-1) + viAngle(j, 1));
        end
    end    
    % Update the state with the Hinfty control & generated disturbance
    x_hinfty(:,t) = Ai*x_hinfty(:,t-1) + Bi*u_hinfty(:,t-1) + w_advers(:,t-1);        
end   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the Model-based Regret & the Total Regret
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Record Cost Difference Regret values at each time and initialize to zero for each model
modelCostDiffRegrets = zeros(1,T+1);    

% Record state & control difference regret values at each time and initialize to zero for each model
modelStateCtrlRegrets = zeros(1,T+1);   

% Loop through all finite set of linear system models
for i = 1:numModels   
    % Loop through the entire time horizon
    for t = 1:T+1   
        % Initialize all contributions to zero
        minmaxInputSum = 0;
        hinftyInputSum = 0;
        minmaxStateSum = 0;
        hinftyStateSum = 0;
        stateDifference = 0;
        controlDifference = 0;
        
        % Compute Regret as Cost Difference
        % compute u'Ru and add it to control sum
        if(t < T+1)
            minmaxInputSum = u_minmax(:,t)'*R*u_minmax(:,t);
            hinftyInputSum = u_hinfty(:,t)'*R*u_hinfty(:,t);            
        end
        % compute x'Qx and add it to state sum
        minmaxStateSum = x_minmax(:,t)'*Q*x_minmax(:,t);
        hinftyStateSum = x_hinfty(:,t)'*Q*x_hinfty(:,t);  
        
        % Record the Cost Difference regret incurred at time t        
        modelCostDiffRegrets(t) = (minmaxStateSum + minmaxInputSum) - (hinftyStateSum + hinftyInputSum);  
        
        % Compute State-Control Difference Regret
        % Compute Control Difference 
        if(t<T+1)
            controlDifference = (u_minmax(:,t) - u_hinfty(:,t))'*R*(u_minmax(:,t) - u_hinfty(:,t));
        end
        % Compute State Difference 
        stateDifference = (x_minmax(:,t) - x_hinfty(:,t))'*Q*(x_minmax(:,t) - x_hinfty(:,t));
        % Obtain State-Control Difference Regret
        modelStateCtrlRegrets(t) = stateDifference + controlDifference;
        
    end
end

% Compute modelbased regret as cumulative sum upto T+1
modelbasedCostDiffRegret = cumsum(modelCostDiffRegrets, 2);

% Compute modelbased regret as cumulative sum upto T+1
modelbasedStateCtrlDiffRegret = cumsum(modelStateCtrlRegrets, 2);

% Compute totalRegret as max of modelbasedRegret
totalCostRegret = max(modelbasedCostDiffRegret(:, end));
fprintf('Difference of sum of squares total regret: %.3f \n', totalCostRegret);

% Compute totalRegret as max of modelbasedRegret
totalStateControlRegret = max(modelbasedStateCtrlDiffRegret(:, end));
fprintf('Sum of square of differences total regret: %.3f \n', totalStateControlRegret);


%% Plot the State-Control Difference regret vs time
Tvec = 0:T;
TDivVec = Tvec + 1;
plotRegrets = zeros(T+1, 1);
tpower = 1;
for t = 1:T+1
    plotRegrets(t,1) = modelbasedStateCtrlDiffRegret(1,t)/(t)^(tpower);    % tpower (0.8, 1) decreases to 0, (0, 0.79)-increases
end
figNum = figNum+1;
figure(figNum);
set(gcf, 'Position', get(0, 'Screensize'));
plot(Tvec, modelbasedStateCtrlDiffRegret, '-.b');
hold on;
plot(Tvec, plotRegrets, '-.r');
xlabel('time');
ylabel('regret');
hold off;
legend('$\mathcal{R}(\bar{\pi}^{\dagger}, \pi^{\star}_{2}, T)$', '$\frac{\mathcal{R}(\bar{\pi}^{\dagger}, \pi^{\star}_{2}, T)}{T}$', 'Location','southeast');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 6);
set(a, 'linewidth', 6);
set(a, 'FontSize', 60);
if(disturbanceSelect == 1)
    saveas(figure(figNum), [currentFolder 'AdverseWRegret.png']);
elseif(disturbanceSelect == 2)
    saveas(figure(figNum), [currentFolder 'ConfuseWRegret.png']);
elseif(disturbanceSelect == 3)
    saveas(figure(figNum), [currentFolder 'SinusoidWRegret.png']);
end

%% Plot the minimax trajectories vs time
figNum = figNum + 1;
figure(figNum);
set(gcf, 'Position', get(0, 'Screensize'));
plot(tvec, x_minmax(:,:)');
hold on;
xlabel('time');
% ylabel('$x^{\bar{\pi}^{\dagger}}_{k}$', 'interpreter', 'latex');
ylabel('states');
hold off;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 6);
set(a, 'linewidth', 6);
set(a, 'FontSize', 60);
if(disturbanceSelect == 1)
    saveas(figure(figNum), [currentFolder 'MACAdverseWStates.png']);
elseif(disturbanceSelect == 2)
    saveas(figure(figNum), [currentFolder 'MACConfuseWStates.png']);
elseif(disturbanceSelect == 3)
    saveas(figure(figNum), [currentFolder 'MACSinusoidWStates.png']);
end
% 
%% Plot the Hinfty trajectories vs time
figNum = figNum + 1;
figure(figNum);
set(gcf, 'Position', get(0, 'Screensize'));
plot(tvec, x_hinfty(:,:)');
hold on;
xlabel('time');
% legendText1 = '$x^{\pi^{\star}_{';
% legendText2 = num2str(modelNum);
% legendText3 = '}}_{k}$'; 
% legendText = strcat(legendText1, strcat(legendText2, legendText3)); 
% ylabel(legendText, 'interpreter', 'latex');
ylabel('states');
hold off;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 6);
set(a, 'linewidth', 6);
set(a, 'FontSize', 60);
if(disturbanceSelect == 1)
    saveas(figure(figNum), [currentFolder 'HinfAdverseWStates.png']);
elseif(disturbanceSelect == 2)
    saveas(figure(figNum), [currentFolder 'HinfConfuseWStates.png']);
elseif(disturbanceSelect == 3)
    saveas(figure(figNum), [currentFolder 'HinfSinusoidWStates.png']);
end

% %% Plot the Difference trajectories and difference of control inputs (both Hinfty and Minimax) vs time
% figNum = figNum + 1;
% Tvec = 0:T-1;
% figure(figNum);
% set(gcf, 'Position', get(0, 'Screensize'));
% plot(Tvec, u_minmax(:,:) - u_hinfty(:,:), '-.k');
% hold on;
% plot(tvec, x_minmax(:,:)' - x_hinfty(:,:)');
% hold on;
% xlabel('time');
% legendText1 = '$x^{\bar{\pi}^{\dagger}}_{k} - x^{\pi^{\star}_{';
% legendText2 = num2str(modelNum);
% legendText3 = '}}_{k}$'; 
% legend1Text = strcat(legendText1, strcat(legendText2, legendText3)); 
% legendText4 = '$u^{\bar{\pi}^{\dagger}}_{k} - u^{\pi^{\star}_{';
% legendText5 = num2str(modelNum);
% legendText6 = '}}_{k}$'; 
% legend2Text = strcat(legendText4, strcat(legendText5, legendText6));  
% legend(legend2Text, '$[x^{\bar{\pi}^{\dagger}}_{k}]_{1} - [x^{\pi^{\star}_{2}}_{k}]_{1}$', '$[x^{\bar{\pi}^{\dagger}}_{k}]_{2} - [x^{\pi^{\star}_{2}}_{k}]_{2}$', '$[x^{\bar{\pi}^{\dagger}}_{k}]_{3} - [x^{\pi^{\star}_{2}}_{k}]_{3}$', 'Location','southeast');
% hold off;
% a = findobj(gcf, 'type', 'axes');
% h = findobj(gcf, 'type', 'line');
% set(h, 'linewidth', 6);
% set(a, 'linewidth', 6);
% set(a, 'FontSize', 60);
% if(disturbanceSelect == 1)
%     filename = 'AdverseWStatesControlsDiff.png';    
% elseif(disturbanceSelect == 2)
%     filename = 'ConfuseWStatesControlsDiff.png';    
% elseif(disturbanceSelect == 3)
%     filename = 'SinusoidWStatesControlsDiff.png';    
% end
% file = fullfile(currentFolder, filename);
% exportgraphics(figure(figNum), file);
%
%% Plot the Hinfinity disturbance inputs vs time
% Tvec = 0:T-1;
% figNum = figNum + 1;
% figure(figNum);
% set(gcf, 'Position', get(0, 'Screensize'));
% plot(Tvec, w_advers(:,:)');   
% xlabel('Time');
% legendText1 = '$w^{\pi^{\star}_{';
% legendText2 = num2str(modelNum);
% legendText3 = '}}_{k}$'; 
% legendText = strcat(legendText1, strcat(legendText2, legendText3)); 
% ylabel(legendText, 'interpreter', 'latex');
% hold on;
% a = findobj(gcf, 'type', 'axes');
% h = findobj(gcf, 'type', 'line');
% set(h, 'linewidth', 6);
% set(a, 'linewidth', 6);
% set(a, 'FontSize', 60);
% if(disturbanceSelect == 1)
%     filename = 'HinfAdverseW.png';
% elseif(disturbanceSelect == 2)
%     filename = 'HinfConfuseW.png';    
% elseif(disturbanceSelect == 3)
%     filename = 'HinfSinusoidW.png';    
% end
% file = fullfile(currentFolder, filename);
% exportgraphics(figure(figNum), file);
%
%
% %% Plot the difference of control inputs (both Hinfty and Minimax) vs time
% Tvec = 0:T-1;
% figNum = figNum + 1;
% figure(figNum);
% set(gcf, 'Position', get(0, 'Screensize'));
% plot(Tvec, u_minmax(:,:) - u_hinfty(:,:,modelNum), '-r');
% hold on;
% xlabel('time');
% legendText1 = '$u^{\bar{\pi}^{\dagger}} - u^{\pi^{\star}_{';
% legendText2 = num2str(modelNum);
% legendText3 = '}}$'; 
% legendText = strcat(legendText1, strcat(legendText2, legendText3));  
% ylabel(legendText);
% a = findobj(gcf, 'type', 'axes');
% h = findobj(gcf, 'type', 'line');
% set(h, 'linewidth', 6);
% set(a, 'linewidth', 6);
% set(a, 'FontSize', 60);
% if(disturbanceSelect == 1)
%     saveas(figure(figNum), [currentFolder 'AdverseWControl.png']);
% elseif(disturbanceSelect == 2)
%     saveas(figure(figNum), [currentFolder 'ConfuseWControl.png']);
% elseif(disturbanceSelect == 3)
%     saveas(figure(figNum), [currentFolder 'SinusoidWControl.png']);
% end

% %% Plot the Cost Difference regret vs time
% Tvec = 0:T;
% plotRegrets = zeros(T+1, 1);
% for i = 1:T+1
%     plotRegrets(i,1) = modelbasedCostDiffRegret(modelNum, i)/i;
% end
% figNum = figNum+1;
% figure(figNum);
% set(gcf, 'Position', get(0, 'Screensize'));
% plot(Tvec, modelbasedCostDiffRegret(modelNum, :), '-.b');
% hold on;
% plot(Tvec, plotRegrets, '-.r');
% xlabel('time');
% ylabel('regret');
% hold off;
% a = findobj(gcf, 'type', 'axes');
% h = findobj(gcf, 'type', 'line');
% legend('$\bar{\mathcal{R}}(\pi^{\dagger}, \pi^{\star}_{2}, T)$', '$\frac{\bar{\mathcal{R}}(\pi^{\dagger}, \pi^{\star}_{2}, T)}{T}$', 'Location','southeast');
% set(h, 'linewidth', 6);
% set(a, 'linewidth', 6);
% set(a, 'FontSize', 60);
% if(disturbanceSelect == 1)
%     saveas(figure(figNum), [currentFolder 'RegretAdverseW.png']);
% elseif(disturbanceSelect == 2)
%     saveas(figure(figNum), [currentFolder 'RegretConfuseW.png']);
% elseif(disturbanceSelect == 3)
%     saveas(figure(figNum), [currentFolder 'RegretSinusoidW.png']);
% end
%
% 
% %% Plot the Minimax disturbance inputs vs time
% Tvec = 0:T-1;
% figNum = figNum + 1;
% figure(figNum);
% set(gcf, 'Position', get(0, 'Screensize'));
% plot(Tvec, w_minimax(:,:)');
% xlabel('Time');
% ylabel('$w^{\bar{\pi}^{\dagger}}_{k}$');
% hold on;
% a = findobj(gcf, 'type', 'axes');
% h = findobj(gcf, 'type', 'line');
% set(h, 'linewidth', 6);
% set(a, 'linewidth', 6);
% set(a, 'FontSize', 60);
% if(disturbanceSelect == 1)
%     saveas(figure(figNum), [currentFolder 'MACAdverseW.png']);
% elseif(disturbanceSelect == 2)
%     saveas(figure(figNum), [currentFolder 'MACConfuseW.png']);
% elseif(disturbanceSelect == 3)
%     saveas(figure(figNum), [currentFolder 'MACSinusoidW.png']);
% end