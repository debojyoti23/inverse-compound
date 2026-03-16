function x = solveLinear_complex(y, n, k)
    % Solves Log Linear system to recover base singular values/ eigenvalues from
    % their compound. Solve the linear system y = Lx
    C = nchoosek(1:n,k);
    L = full(sparse(repmat(1:size(C,1),1,k),C(:),1,size(C,1),n));
    % magnitude
    y_abs = abs(y);
    y_abs = log(y_abs);
    x_abs = L \ y_abs;
    x_abs = exp(x_abs);
    % phase
    y_angle = angle(y);
    x_angle = L \ y_angle;
    x = x_abs .* exp(1i * x_angle);
end

