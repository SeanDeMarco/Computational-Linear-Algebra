function R = cholesky_tridiagonal(A)
    main_diag = diag(A); 
    subdiag = diag(A, -1); 

    n = length(main_diag); 
    R = zeros(n, n); 

    for i = 1:n
        % Computing diagonal entry 
        if i == 1
            R(i, i) = sqrt(main_diag(i));
        else
            R(i, i) = sqrt(main_diag(i) - R(i-1, i)^2);
        end

        % Computing superdiagonal entry 
        if i < n
            R(i, i+1) = subdiag(i) / R(i, i);
        end
    end
end