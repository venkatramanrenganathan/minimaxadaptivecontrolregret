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



load('RandomSystemData.mat');
load('RandomSystemMinimaxData.mat');
Q = eye(3);         % State penalty matrix 
R = eye(1);         % Input penalty matrix 
startGamma = 100;   % Starting gain from disturbance to error
gammaOptFlag = 0;

% Initialize the control gains, disturbance gains and worst case gamma 
K_Gains = cell(4, 1);
F_Gains = cell(4, 1);
Gammas  = zeros(4, 1);

Hinf_gammas = [1.266 4.544 2.913 2.298];
MBSubGaps = MAC_gamma - Hinf_gammas

minGap = MAC_gamma - max(Hinf_gammas)
maxGap = MAC_gamma - min(Hinf_gammas)

% Get H-infinity Control & adversarial disturbance via Dynamic Games Approach - USE MINIMAX GAMMA
% i = 4;
% [K_Gains{i,1}, F_Gains{i,1}, Gammas(i,1)] = HinfDynamicGame(AMatrices{i}, BMatrices{i}, Q, R, MAC_gamma, gammaOptFlag);