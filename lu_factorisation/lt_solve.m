function [b] = lt_solve(LU,b)
    n = size(LU, 1);        
    %L is contained in the LU matrix therefore, it does not need to be
    %extracted
    % Forward Substitution algorithm based on algorithm 4.2 (pg26 in notes)
    % Normalisation of b(i) is skipped since diagonal of the L matrix is
    %always 1
    for i = 2:n
        for j = 1:i-1
            b(i) = b(i) - LU(i,j)*b(j);
        end
    end
end
