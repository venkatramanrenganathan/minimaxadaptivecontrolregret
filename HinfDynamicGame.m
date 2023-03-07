function [Kmatrix, Fmatrix, gamma] = HinfDynamicGame(statespace, Q, R)
% HinfDynamicGame solves the H-infinity Optimal Control via Dynamic Games
% given the state space system and the penalty matrices - Q, R.


A = statespace.A;
B = statespace.B;
C = statespace.C;
D = statespace.D;

n = size(A, 1);
m = size(B, 2);
T = 1000;
invR = inv(R);
gamma = 50;


% Define UB and LB used in the bisection.  UB will always be the largest
% feasible value of gamma. 
UB = gamma;
LB = 0;
tol = 1e-8;
gammaFail = 0;

% Solve the infinite horizon H_infty Dynamic Game
while (UB - LB >= tol)
    
    % bisect
    gamma = (UB+LB)/2;
    
    Rext = [R, zeros(m, n); zeros(n, m), -gamma^2*eye(n)];
    Bext = [B, eye(n)];
    [P, K, ~] = idare(A, Bext, Q, Rext, [],[]);
    
    % 
    if isequal(K, [])
        UB = gamma;
    else
        LB = gamma;
    end
end


Kmatrix = -inv(R + B'*P*B + B'*P*inv(gamma^2*eye(n) - P)*P*B)*(B'*P*A + B'*P*inv(gamma^2*eye(n) - P)*P*A);
Fmatrix = inv(P - gamma^2*eye(n) - P*B*inv(R+B'*P*B)*B'*P)*(P*A - P*B*inv(R+B'*P*B)*B'*P*A);

% % Solve the finite horizon H_infty Dynamic Game
% while (UB - LB >= tol)
%     
%     % perform gamma bisection
%     gamma = (UB+LB)/2;
%     
%     % placeholders for recursion matrices
%     Kmatrices = zeros(m, n, T);    % control gain matrices
%     Fmatrices = zeros(n, n, T+1);  % disturbance gain matrices
%     Lmatrices = zeros(n, n, T);    % lambda matrices
%     Mmatrices = zeros(n, n, T+1);  % M matrices  
%     
%     % populate the terminal condition
%     Mmatrices(:,:, T+1) = Q;
%     
%     % Loop backwards in time
%     for k = T:-1:1 
%         GammaSqr = gamma*gamma;
%         invGammaSqr = 1/(GammaSqr);        
%         Lmatrices(:,:, k) = eye(n) + (B*invR*B' - invGammaSqr*eye(n))*Mmatrices(:,:, k+1);
%         Mmatrices(:,:, k) = Q + A'*Mmatrices(:,:, k+1)*inv(Lmatrices(:,:, k))*A;        
%         try chol(GammaSqr*eye(n) - Mmatrices(:,:, k));
%             gammaFail = 0;        
%         catch ME
%             disp('Matrix is not symmetric positive definite');
%             gammaFail = 1;
%             break;
%         end
%         % Record the control and disturbance gain matrices
%         Kmatrices(:,:, k) = -invR*B'*Mmatrices(:,:, k+1)*inv(Lmatrices(:,:, k))*A;
%         Fmatrices(:,:, k) = invGammaSqr* Mmatrices(:,:, k+1)*inv(Lmatrices(:,:, k))*A;        
%     end
%     
%     % Change upperbound and lowerbound accordingly for bisection
%     if(gammaFail == 0)
%         UB = gamma;
%     else
%         LB = gamma;
%     end
%     
% end
% 
% % Report the converged first K and L matrices
% Kmatrix = Kmatrices(:,:, 1);
% Fmatrix = Fmatrices(:,:, 1);

end

