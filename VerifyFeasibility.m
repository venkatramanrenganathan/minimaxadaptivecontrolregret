function [feasible] = VerifyFeasibility(PMatrices, KMatrices, AMatrices, BMatrices, Q, R, gamma)
% Checks if the proposed solution is feasible. 
% Output:     feasible - 1 if the solution is feasible, 0 otherwise. 

     N = length(AMatrices);
     [n, ~] = size(BMatrices{1});
     feasible = true;
     feasibilityTol = 1e-9;
     
     % Check if 0 < Pij < gamma^2 I
     for i = 1:N
         for j = 1:N
            evals = eig(PMatrices{i, j});
            if ((min(evals) <= 0 - feasibilityTol) || (max(evals) >= gamma^2 - feasibilityTol))
               feasible = false; 
               disp('FALSE SIMPLE CONSTRAINTS')
            end
         end
     end
     
     % Check difficult constraints 
     for i = 1:N
         for j = 1:N
             for k = 1:N
                 Pik = PMatrices{i, k};
                 Pij = PMatrices{i, j};
                 Kk = KMatrices{k};
                 Aik = AMatrices{i} - BMatrices{i}*Kk;
                 Ajk = AMatrices{j} - BMatrices{j}*Kk;
                 LHS = Pik - Q - Kk'*R*Kk - 1/4*(Aik + Ajk)'*((inv(Pij) - gamma^(-2)*eye(n))\(Aik + Ajk)) + gamma^2/4*(Aik - Ajk)'*(Aik - Ajk);
                 evals = eig(LHS);
                 if min(evals) <= 0 - feasibilityTol
                    feasible = false;
                    disp('FALSE DIFFICULT CONSTRAINTS')
                 end
             end
         end
     end
end


