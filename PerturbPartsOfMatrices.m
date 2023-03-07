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