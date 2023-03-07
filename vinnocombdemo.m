%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code simulates the control of a doible integrator model using the 
% minimax adaptive control algorithm.
%
% Copyrights Authors: 1) Venkatraman Renganathan - Lund University, Sweden.
%                     2) Anders Rantzer - Lund University, Sweden.
%
% Email: venkatraman.renganathan@control.lth.se
%
% Courtesy: 1.) Daniel Cedelberg - Linköping University, Sweden.
%           2.) Anders Hansson - Linköping University, Sweden.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a fresh start
clear all; close all; clc;

% Define system parameters of vinnicombe example
A = 1;
B = 1;
C = 1;
D = 0;

% Define the finite set of linear system models
AMatrices = {A, A};
BMatrices = {B, -B};

% Infer the system dimensions
[n,m] = size(B);

% Get to know the number of linear system models available
numModels = length(BMatrices);

% Form the state-space with difference linear models
Plants = cell(numModels, 1);
for i = 1:numModels
    Plants{i, 1} = ss(AMatrices{i}, BMatrices{i}, C, D);
end

% Compute and display the Vinnicomb gain
vinnicomb_gain = A + sqrt(A*A + 1);
fprintf('Vinnicomb Gain = %.4f \n', vinnicomb_gain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Minimax Adaptive Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
% Define design parameters.
Q = 0.1;         % State penalty matrix 
R = 0.1;         % Input penalty matrix 
startGamma = 2.5;   % Starting gain from disturbance to error

% Get minimax adaptive control (MAC) policies
%[MAC_status, MAC_gamma, MAC_PMatrices, MAC_KMatrices] = Approach1(AMatrices, BMatrices, Q, R, startGamma);
[MAC_status, MAC_gamma, MAC_PMatrices, MAC_KMatrices] = Approach3(AMatrices, BMatrices, Q, R, startGamma);
disp(MAC_PMatrices{1} < (MAC_gamma*MAC_gamma));
disp(MAC_PMatrices{2} < (MAC_gamma*MAC_gamma));
fprintf('Minimax Gain = %.4f \n', MAC_gamma);
disp('Finished Computing Minimax Adaptive Control Gains');


%% Helper Functions

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [status, gamma, PMatrices, KMatrices] = Approach2(AMatrices, BMatrices, Q, R, startGamma)
% APPROACH2:  See paper for details of how approach 2 works.
% Input:     AMatrices, BMatrices - The different possible (A, B)-pairs.
%            Q, R                 - Matrices used to define the cost function.
%            startGamma           - Initial value of gamma used in the
%                                   bisection.
% Output:    status - Yalmip status code.
%            gamma, PMatrices, KMatrices - Found value of gamma and
%                                          corresponding matrices.


    disp('Running approach 2');
    
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
        [status, PMatricesTemp, KMatricesTemp] = checkFeasibilityApproach2(AMatrices, BMatrices, Q, R, UB);
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
            
            % Try approach 2 for gamma equals mid.
            [statusTemp, PMatricesTemp, KMatricesTemp] = checkFeasibilityApproach2(AMatrices, BMatrices, Q, R, mid);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [AMatrices, BMatrices] = PerturbPartsOfMatrices(ANominal, BNominal, epsilon, N)
% PERTURBPARTSOFMATRICES: Perturbs all nonzero elements of the nominal
% system description.
% Input:       ANominal, BNominal - Nominal system description.
%              epsilon - The perturbation level. 
%              N - Number of (A,B)-pairs that should be generated.
% Output:      AMatrices, BMatrices - N different (A,B)-pairs.

    [n, m] = size(BNominal);
    
    AMatrices = cell(1, N);
    BMatrices = cell(1, N);
    
    zeroElementsANominal = (ANominal == 0);
    zeroElementsBNominal = (BNominal == 0);
    
    AMatrices{1} = ANominal;
    BMatrices{1} = BNominal;
    
    for i = 2:N
       % Generate random matrices to perturb with. Note the scaling factor.
       APerturbation = (max(abs(ANominal), [], 'all')*(rand(n, n)*2 - ones(n, n)));
       BPerturbation = (max(abs(BNominal), [], 'all')*(rand(n, m)*2 - ones(n, m)));
       
       % Zero elements in A and B should not be perturbed.
       APerturbation(zeroElementsANominal) = 0;
       BPerturbation(zeroElementsBNominal) = 0;
       
       AMatrices{i} = ANominal + epsilon*APerturbation ;
       BMatrices{i} = BNominal + epsilon*BPerturbation;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [status, PiiMatrices, KiMatrices] = SolveRiccatiEquation(AMatrices, BMatrices, Q, R, gamma)
% SOLVERICCATIEQUATION: Solves the minmax equation to calculate Pii
% and Ki for every i. 
% Input:   AMatrices  - Cell array, each cell corresponding to one Ai. 
%          BMatrices  - Cell array, each cell corresponding to one Bi.
% Output:  PiiMatrices - Cell array, each cell corresponding to one Pi.
%          KiMatrices  - Cell array, each cell corresponding to one Ki. 

        N = length(AMatrices);
        [n, m] = size(BMatrices{1});
        
        PiiMatrices = cell(1, N);
        KiMatrices = cell(1, N);

        Rext = [R, zeros(m, n); zeros(n, m), -gamma^2*eye(n)];
        
        
        for i = 1:N
           Bext = [BMatrices{i}, eye(n)];
           [Pii, Kext, ~] = idare(AMatrices{i}, Bext, Q, Rext, [],[]);
           % If the Riccati equations don't have a solution we represent it 
           % by letting the variable 'status' equal 1.
           if isequal(Kext, [])
             status = 1;
             return;
           end
           PiiMatrices{i} = Pii;
           KiMatrices{i} = Kext(1:m, :); 
        end
        
        % Status equals 0 corresponds to that the Ricatti equations were
        % succesfully solved.
        status = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%