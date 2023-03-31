function [Kmatrix, Fmatrix, gamma] = HinfDynamicGame(A, B, Q, R, gamma, gammaOptFlag)
% HinfDynamicGame solves the H-infinity Optimal Control via Dynamic Games
% given the state space matrices A, B and the penalty matrices - Q, R.

% Infer the dimensions
n = size(A, 1);
m = size(B, 2);

if(gammaOptFlag == 1)
    % Compute the bisection gamma and control gains
    % Initialize a starting disturbance attenuation level
    % gamma = 100;

    % Define UB and LB used in the bisection.  
    UB = gamma;
    LB = 0;

    % Set the tolerance for bisection to converge
    tol = 1e-8;

    % Solve the infinite horizon H_infty Dynamic Game
    while (UB - LB >= tol)

        % Bisect the gamma
        gamma = (UB+LB)/2;

        % Augmented Input penalty matrix for both u and w players
        Rext = [R, zeros(m, n); zeros(n, m), -gamma^2*eye(n)];

        % Augmented Inout matrix
        Bext = [B, eye(n)];

        % Solve the Hinfinity Riccatti
        [P, K, ~] = idare(A, Bext, Q, Rext, [],[]);

        % Set the upper and lower limits for bisection
        if isequal(K, [])
            disp('K is Empty');
            disp(gamma)
            UB = gamma;
        else
            disp('K is not Empty');
            disp(gamma)
            LB = gamma;
        end
    end
    

else
    % Just compute the gains for the given gamma
    % Augmented Input penalty matrix for both u and w players
    Rext = [R, zeros(m, n); zeros(n, m), -gamma^2*eye(n)];

    % Augmented Inout matrix
    Bext = [B, eye(n)];

    % Solve the Hinfinity Riccatti
    [P, K, ~] = idare(A, Bext, Q, Rext, [],[]);    
end


% Extract the feedback gains for u and w players from the converged gamma
Kmatrix = -inv(R + B'*P*B + B'*P*inv(gamma^2*eye(n) - P)*P*B)*(B'*P*A + B'*P*inv(gamma^2*eye(n) - P)*P*A);
Fmatrix = -inv(P - gamma^2*eye(n) - P*B*inv(R+B'*P*B)*B'*P)*(P*A - P*B*inv(R+B'*P*B)*B'*P*A);


% % Solve the finite horizon H_infty Dynamic Game
% T = 1000;
% invR = inv(R);
% gammaFail = 0;
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

