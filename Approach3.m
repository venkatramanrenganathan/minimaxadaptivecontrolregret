function [status, gamma, PMatrices, KMatrices] = Approach3(AMatrices, BMatrices, Q, R, startGamma)
% APPROACH3: See paper for details of how approach 3 works.
% Input:     AMatrices, BMatrices - The different possible (A, B)-pairs.
%            Q, R                 - Matrices used to define the cost function.
%            startGamma           - Initial value of gamma used in the
%                                   bisection.
% Output:    status - Yalmip status code.
%            gamma, PMatrices, KMatrices - Found value of gamma and
%                                          corresponding matrices.

      disp('Running approach 3')
      [~, newGamma, ~, KMatrices1] = Approach2(AMatrices, BMatrices, Q, R, startGamma);
    
     
      % Linearize and solve the linearization using bisection. Repeat until
      % a new linearization does not lower the value of gamma.
      KMatrices = KMatrices1;
      oldGamma = newGamma + 10;  % Just so it enters the loop the first time.
      
      while abs(newGamma - oldGamma) > 0.05
          oldGamma = newGamma;
         [status, newGamma, PMatrices, KMatrices] = BisectionApproach3(AMatrices, BMatrices, Q, R, KMatrices, oldGamma);
      end
      gamma = newGamma;
end

