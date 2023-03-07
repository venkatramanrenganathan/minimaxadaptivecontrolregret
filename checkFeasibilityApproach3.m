function [status, PMatrices, KMatrices]= checkFeasibilityApproach3(AMatrices, BMatrices, Q, R, linearizationPointKMatrices, gamma)
% CHECKFEASIBILITYAPPORACH3: Checks if the feasibility problem is feasible using
% approach 3.
% Input:   AMatrices - Cell array of size 1xN with N system matrices Ai.
%          BMatrices - Cell array of size 1xN with N state-input matrices Bi.
%          linearizationPointKMatrices - The linearization point. 
% Output:  status - Integer. Contains Yalmip status code.
%          PMatrices - Cell array of size NxN with the matrices Pij.
%          KMatrices - Cell array of size 1xN with the matrices Ki. 


     N = length(AMatrices);
     [n, m] = size(BMatrices{1});
      
     % Define variables. We store the matrices Pij in a 2D cell array
     % denoted by 'P'. 
     if N ~= 1
        P = reshape(sdpvar(n*ones(N^2), n*ones(N^2)), N, N);
        K = sdpvar(m*ones(N), n*ones(N));
     else
        P = cell(1, 1);
        K = cell(1, 1);
        P{1, 1} = sdpvar(n, n);
        K{1, 1} = sdpvar(m, n);
     end
     
     % Add simple constraints
     constraints = [];
     for i = 1:N
         for j = 1:N
             constraints = [constraints, gamma^2*eye(n) - P{i, j} >= 0, P{i, j} >= 0];
         end
     end
     
     % Add complicated constraints.
      for i = 1:N
         Ai = AMatrices{i};
         Bi = BMatrices{i};
         for j = 1:N
             Aj = AMatrices{j};
             Bj = BMatrices{j};
             for k = 1:N
                Aik = Ai - Bi*K{k};
                Ajk = Aj - Bj*K{k};
                S = Ai'*Ai + Aj'*Aj - (Ai'*Bi + Aj'*Bj)*K{k} - K{k}'*(Bj'*Aj + Bi'*Ai);
                
                Yijk = sqrtm(Bi'*Bi+Bj'*Bj)*linearizationPointKMatrices{k};
                Dsqrt = sqrtm(Bi'*Bi + Bj'*Bj);
                Lijk = K{k}'*Dsqrt*Yijk + Yijk'*Dsqrt*K{k} - Yijk'*Yijk;
                
                schurMatrix = [P{i, k} - Q + gamma^2/2*S + gamma^2/2*Lijk,     gamma^2/2*(Aik + Ajk)',        K{k}';
                               gamma^2/2*(Aik + Ajk),         gamma^2*eye(n) - P{i, j},                        zeros(n, m);
                               K{k},                          zeros(m, n),                                     inv(R)];
                constraints = [constraints, schurMatrix >= 0];
             end
         end
      end
      
    
     % Solve the feasibility problem.
    ops = sdpsettings('solver','mosek','verbose',0);
    sol = optimize(constraints, [], ops);
    
    % Yalmip returns zero if the problem is feasible.
    status = sol.problem;
    
    if ((status ~= 0) && (status ~= 1))
        disp(['Numerical issues in approach 3. The Yalmip status code is: ', yalmiperror(status)])
    end
     
    % Save the calculated matrices.
    PMatrices = cellfun(@value, P, 'UniformOutput', false);
    KMatrices = cellfun(@value, K, 'UniformOutput', false);

end