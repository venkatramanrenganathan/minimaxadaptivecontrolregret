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