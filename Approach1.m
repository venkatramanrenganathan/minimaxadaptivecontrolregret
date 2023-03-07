function [status, gamma, PMatrices, KMatrices] = Approach1(AMatrices, BMatrices, Q, R, startGamma)
% APPROACH1: See paper for details of how approach 1 works.
% Input:     AMatrices, BMatrices - The different possible (A, B)-pairs.
%            Q, R                 - Matrices used to define the cost function.
%            startGamma           - Initial value of gamma used in the
%                                   bisection.
% Output:    status - Consistent with Yalmip status code.
%            gamma, PMatrices, KMatrices - Found value of gamma and
%                                          corresponding matrices.

    disp('Running approach 1');
    
    % Initialize UB and LB for bisection. UB will always be the largest
    % feasible value of gamma.
    UB = startGamma;
    LB = 0;
    
    % Define the tolerance used in the bisection. When UB - LB < tol we are
    % done.
    tol = 0.01;
    
    % Define containers in which to store P and K matrices.
    PMatrices = [];
    KMatrices = [];
    
    % If the initial value of gamma is "too large", the feasibility problem
    % becomes ill-conditioned which causes numerical issues for Mosek. We
    % reduce the initial value until Mosek returns feasible/infeasible
    % problem instead of numerical issues. Status code 4 in Yalmip
    % corresponds to numerical issues. Status code 9 in Yalmip corresponds
    % to unknown problems.
    status = 4;
    while (status == 4) || (status == 9) 
        [status, PMatricesTemp, KMatricesTemp] = checkFeasibilityApproach1(AMatrices, BMatrices, Q, R, UB);
        UB = UB*0.98;
    end
    UB = UB/0.98;
    
    disp(['First value of gamma that did not cause numerical issues: ', num2str(UB)]);
     
    % Status equals 1 corresponds to 'Infeasible problem'.
    if status == 1  
       gamma = -1;
       return; 
    % Status equals 9 corresponds to 'Unknown problem in solver'.
    elseif status == 9
       gamma = -9;
       return;
    % Status equals 4 corresponds to 'Numerical problems'.
    elseif status == 4
       gamma = -4;
       return;
    % Status equals 0 corresponds to feasible problem.
    elseif status == 0
            PMatrices = PMatricesTemp;
            KMatrices = KMatricesTemp;
    end
    
    
    % Apply bisection.
    while (UB - LB >= tol)
            mid = (UB+LB)/2;
            
            % Try approach 1 for gamma equals mid.
            [statusTemp, PMatricesTemp, KMatricesTemp] = checkFeasibilityApproach1(AMatrices, BMatrices, Q, R, mid);
            %disp(['Current value of gamma: ', num2str(mid)]);
            %disp(['Status: ', num2str(statusTemp)]);
            
            % If feasible optimization problem.
            if statusTemp == 0
               UB = mid;
               PMatrices = PMatricesTemp;
               KMatrices = KMatricesTemp;
               
            % If Mosek returns something else than feasible optimization
            % problem. This could be infeasible problem/numerical issues/unknown problem.  
            else
               LB = mid;
            end
    end
    
    % After the bisection we know that UB is a feasible value of gamma. 
    gamma = UB;
end