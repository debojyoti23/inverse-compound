function [Uk_rcp, k_out] = compute_reciprocal(Uk,n,k)
    % this function computes C_k(U) from given C_{n-k}(U) where U is
    % non-singular square matrix and full rank using the formula C_k(A) =
    % (detA) SP(C_{n-k}^T)^{-1} PS

    C = nchoosek(1:n,k);
    P = fliplr(eye(size(Uk)));
    S = diag((-1).^sum(C,2));
    detU = nthroot(det(Uk),nchoosek(n-1,k-1));
    Uk_rcp = detU*S*P*(Uk'\eye(size(Uk')))*P*S;
    k_out = n-k;
end

