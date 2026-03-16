clear;
clc;

%-------------------- Input Square base matrix A------------------
n = 6; 
k = 4;
% A = diag([2 2 3 1])
% A = rand(n);
% Create A with real eigenvalues and eigenvector. A full rank and
% diagonalizable. Non-distinct eigenvalues might fail
U0 = rand(n);
U0 = normalize(U0,'norm',2);
% D0 = diag(rand(n,1));
D0 = diag([5 1 1 -1 2 6]);
A = U0 * D0 * inv(U0)
%-----------------------------------------------------------------

tol = 1e-10;
[n,~] = size(A);
[U, D1] = eig(A);

CA = compound(A,k);

[V, D2] = eig(CA);
[~,idx] = sort(abs(diag(D2)),'descend');
V = V(:,idx);
d2 = diag(D2);
D2 = diag(d2(idx));

if ~ (k>1 & k<n)
    error("k must be > 1 and < %d",n);
end

% Check if A full rank
if min(abs(diag(D2))) < tol
    error('Error! Base matrix A is not full rank')
end

% [U_rec, W_Mat] = decompose_Comp_EigVec(V, n, k);
[U_rec, W_Mat] = decompose_Comp_EigVec_1(V, n, k, n);

% Recover a correct order of the compound eigenvectors
P = zeros(1,size(V,2));
V_rec = normalize(compound(U_rec, k),'norm',2);
d2_rec = diag((CA * V_rec)' * V_rec);
D2_rec = diag(d2_rec);

% Recover absolute eigenvalues by solving a linear system in the 
% log scale
y = abs(diag(D2_rec));
% Handle zeros and non-zeros separately
if min(y) > 1e-8
    x = solveLinear(y,n,k);
    D1_rec = diag(x);
    % Recover the eigenvalue signs (which makes sense when real)
    if min(real(diag(D2_rec))) < 0
        y_sgn = sign(real(diag(D2_rec)));
        y_sgn = (y_sgn + 1) / 2;
        [x_sgn, ~] = solve_mod2(y_sgn,n,k);
        x_sgn = x_sgn * 2 - 1;
        D1_rec = diag(x_sgn) * D1_rec;
    end
else
    y_nz = y(y > 1e-8);
    if length(y_nz) <= 1
        % happens if x has >= n-k zeros
        error('Linear system not solvable due to <= %d non zero eigenvalues!', k)
    end
    for z = 1 : n-k-1
        if length(y_nz) == nchoosek(n-z,k)
            break
        end
    end
    x_nz = solveLinear(y_nz,n-z,k);

    % Recover the eigenvalue signs (which makes sense when real)
    if min(real(diag(D2_rec))) < 0
        y_sgn = sign(real(diag(D2_rec)));
        y_sgn = y_sgn(y > tol);
        y_sgn = (y_sgn + 1) / 2;
        [x_sgn, ~] = solve_mod2(y_sgn,n-z,k);
        x_sgn = x_sgn * 2 - 1;
        x_nz = x_nz .* x_sgn;
    end

    % Insert zeros back into x_nz
    C1 = nchoosek(1:n,n-z);
    for i=1:size(C1,1)
        d1 = zeros(n,1);
        d1(C1(i,:)) = x_nz;
        D1_rec = diag(d1);
        D12 = compound(D1_rec,k);
        if norm(diag(D2_rec) - diag(D12), 2) < 1e-8
            break
        end
    end
end

%--------- Finally, recover the original matrix ------------------
A_rec = U_rec * D1_rec * inv(U_rec)
%-----------------------------------------------------------------