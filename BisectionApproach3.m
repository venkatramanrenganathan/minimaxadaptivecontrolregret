function [status, gamma, PMatrices, KMatrices] = BisectionApproach3(AMatrices, BMatrices, Q, R, InitialK, startGamma)
% BISECTIONAPPROACH3: Applies bisection to solve the linearized problem. 
% Input:     AMatrices, BMatrices - The different possible (A, B)-pairs.
%            Q, R                 - Matrices used to define the cost function.
%            InitialK             - K-matrices used in the linearization.
%            startGamma           - Initial value of gamma used in the
%                                   bisection.
% Output:    status - Consistent with Yalmip status codes.
%            gamma, PMatrices, KMatrices - Found value of gamma and
%                                          corresponding matrices.

   PMatrices = [];
   KMatrices = [];
   
   % Define tolerance used in the bisection.  When UB - LB < tol we are
   % done.
   tol = 0.01;
   
   % Define UB and LB used in the bisection.  UB will always be the largest
   % feasible value of gamma. 
   UB = startGamma;
   LB = 0;
   
   % Check that the start value of gamma actually gives a feasible problem.
   [status, PMatricesTemp, KMatricesTemp] = checkFeasibilityApproach3(AMatrices, BMatrices, Q, R, InitialK, UB);
   
   % Gamma equals 1 indicates that approach 3 found no feasible value of
   % gamma. This scenario has never occurred in any of our experiments.
   if status ~= 0
      gamma = 1;
      return; 
   end
   PMatrices = PMatricesTemp;
   KMatrices = KMatricesTemp;
   
   while (UB - LB >= tol)
            mid = (UB+LB)/2;
            [statusTemp, PMatricesTemp, KMatricesTemp] = checkFeasibilityApproach3(AMatrices, BMatrices, Q, R, InitialK, mid);
            
            % If feasible problem.
            if statusTemp == 0
               UB = mid;
               PMatrices = PMatricesTemp;
               KMatrices = KMatricesTemp;
            else
               LB = mid;
            end
   end
   
   % We know that UB is a valid value of gamma.
   gamma = UB;
end