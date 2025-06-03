function[H,U] = hessenberg(A)
    H = A;
    n = size(A,1);
    U = eye(n);
    for k = 1:n-2
        x = H(k+1:n, k);
        if x(1) >= 0
            v = x + norm(x) * eye(length(x), 1);
        else
            v = x - norm(x) * eye(length(x), 1);
        end
        v = v/norm(v, 2);
        R = eye(length(x)) - 2 * v * transpose(v);
        H(k+1:n, k:n) = transpose(R) * H(k+1:n, k:n);
        H(1:n, k+1:n) = H(1:n, k+1:n) * R;
        U(k+1:n, :) = R * U(k+1:n, :);
    end
end