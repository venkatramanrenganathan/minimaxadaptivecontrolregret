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