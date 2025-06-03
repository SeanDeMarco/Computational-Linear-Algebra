function [b] = ut_solve(LU, b)
    n = size(LU, 1);    
    %U matrix is contained in the LU matrix therefore it does not need to
    %be explicitally extracted
    % Backward substitution algorithm based on algorithm 4.4 (pg29 in
    % notes)
    b(n) = b(n) / LU(n, n);
    for i = n-1:-1:1
        for j = i+1:n
            b(i) = b(i) - LU(i, j) * b(j);
        end
        b(i) = b(i) / LU(i, i);
    end
end