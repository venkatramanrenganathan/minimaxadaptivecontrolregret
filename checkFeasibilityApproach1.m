function [status, PMatrices, KMatrices] = checkFeasibilityApproach1(AMatrices, BMatrices, Q, R, gamma)
% checkFeasibilityApproach1: Checks if the feasibility problem is feasible using
% approach 1.
% Input:   AMatrices - Cell array of size 1xN with N system matrices Ai.
%          BMatrices - Cell array of size 1xN with N state-input matrices Bi.
% Output:  status - Consistent with Yalmip status codes.
%          PMatrices - Cell array of size NxN with the matrices Pij.
%          KMatrices - Cell array of size 1xN with the matrices Ki. 
    
     N = length(AMatrices);
     [n, ~] = size(BMatrices{1});
     PMatrices = [];
     KMatrices = [];
      
     % Calculate Pii and Ki for i = 1,..., N.
     [status, PiiMatrices, KiMatrices] = SolveRiccatiEquation(AMatrices, BMatrices, Q, R, gamma);
     
     % Check if the Ricatti equation was not succesfully solved.
     if status == 1
        disp(['No solution to the Ricatti equation for gamma equal to ', num2str(gamma)]);
        return;
     end
     
     % Define variables. We store the matrices Pij in a 2D cell array
     % denoted by 'P'. We have already calculated the value of the diagonal
     % elements in this cell array.
     P = reshape(sdpvar(n*ones(N^2), n*ones(N^2)), N, N);
     P(logical(eye(N, N))) = PiiMatrices;
     
     % Add 'complicated' constraints.
     constraints = [];
    
     for i = 1:N
         for j = 1:N
             for k = 1:N
              
                  % For i = j = k we shall have no constraints since this
                  % constraint is automatically satisfied when solving the
                  % Ricatti equation.
                  if ~((k == i) && (j == i))  
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
     end
     
     % Add 'simple' constraints.
     for i = 1:N
         for j = 1:N
            if i ~= j
                constraints = [constraints, gamma^2*eye(n) - P{i, j} >= 0, P{i, j} >= 0];
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
        disp(['Numerical issues in approach 1. The Yalmip status code is: ', yalmiperror(status)])
     end
        
        
      % Convert P matrices from sdpvar-objects to actual numerical
      % matrices.
      PMatrices = cellfun(@value, P, 'UniformOutput', false);
      KMatrices = KiMatrices;
    
end