function [status, PMatrices, KMatrices] = checkFeasibilityApproach2(AMatrices, BMatrices, Q, R, gamma)
% checkFeasibilityApproach2: Checks if the feasibility problem is feasible using
% approach 2. We only fix the K-matrices by solving a Riccati equation, and
% find all P-matrices using optimization.

     N = length(AMatrices);
     [n, ~] = size(BMatrices{1});
     PMatrices = [];
     KMatrices = [];
      
     % Calculate Ki for i = 1,..., N.
     [status, ~, KiMatrices] = SolveRiccatiEquation(AMatrices, BMatrices, Q, R, gamma);
     
     % Check if the Ricatti equation was not succesfully solved.
     if status == 1
        disp(['No solution to the Ricatti equation for gamma equal to ', num2str(gamma)]);
        return;
     end
     
     % Define variables. We store the matrices Pij in a 2D cell array
     % denoted by 'P'.
     P = reshape(sdpvar(n*ones(N^2), n*ones(N^2)), N, N);
     
     constraints = [];
     % Add 'simple' constraints.
     for i = 1:N
         for j = 1:N
                constraints = [constraints, gamma^2*eye(n) - P{i, j} >= 0, P{i, j} >= 0];
         end
     end
     
     % Add 'complicated' constraints.
     for i = 1:N
         for j = 1:N
             for k = 1:N
              
                 Kk = KiMatrices{k};
                 Aik = AMatrices{i} - BMatrices{i}*Kk;
                 Ajk = AMatrices{j} - BMatrices{j}*Kk;
                 Pik = P{i, k};
                 Pij = P{i, j};  
                 schurMatrix = [Pik - Q - Kk'*R*Kk + gamma^2/2*(Aik'*Aik + Ajk'*Ajk),      gamma^2/2*(Aik + Ajk)';
                                gamma^2/2*(Aik + Ajk),                                     gamma^2*eye(n) - Pij];
                            
                 constraints = [constraints, schurMatrix >= 0];                                 
             end
         end
     end
     
     
     
      % Solve feasibility problem
      ops = sdpsettings('solver','mosek','verbose', 0);
      sol = optimize(constraints, [], ops);

      % Yalmip returns 0 if the problem is feasible, 1 if the problem 
      % is infeasible etc.
      status = sol.problem;
     if ((status ~= 0) && (status ~= 1))
        disp(['Numerical issues in approach 2. The Yalmip status code is: ', yalmiperror(status)])
     end
        
        
      % Convert P matrices from sdpvar-objects to actual numerical
      % matrices.
      PMatrices = cellfun(@value, P, 'UniformOutput', false);
      KMatrices = KiMatrices;
    
end