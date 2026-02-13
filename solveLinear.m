function x = solveLinear(y, n, k)
    % Solves Log Linear system to recover base singular values/ eigenvalues from
    % their compound. Solve the linear system y = Lx
    C = nchoosek(1:n,k);
    L = full(sparse(repmat(1:size(C,1),1,k),C(:),1,size(C,1),n));
    y = log(y);
    x = L \ y;
    x = exp(x);
end

