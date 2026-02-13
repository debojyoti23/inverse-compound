clear;
clc;

%----------------- Input Rectangular Base Matrix A ----------------
n = 7; m = 6;
k = 3; % Keep k <= min(m,n)/2
A = rand(n,m)

% U0 = rand(n);
% U0 = normalize(U0,'norm',2);
% D0 = diag([5 0 -9 -6 8 1]);
% A = U0 * D0 * inv(U0)
%-----------------------------------------------------------------


[n,m] = size(A);
[U1, S1, W1] = svd(A);

CA = compound(A,k);
[U2, S2, W2] = svd(CA);

% Compute the tolerance to account for the k-compound operation
% det(A) is computed by LU decomposition: PA = LU for permutation matrix P
% det(A) = det(P)det(L)det(U) where det(P) is +1 or -1, det(L) is 1 and
% det(U) = \prod U_ii. Thus the error propagated in det is approximately
% upper bounded by dim(A) * \eps * max(A(:))
tol = 1e-15 * k * max(CA(:));


% Recover original singular vectors upto (possibly inconsistent) 
% sign flips
U1_rec = decompose_Comp_EigVec(U2,n,k);
W1_rec = decompose_Comp_EigVec(W2,m,k);

% Recover a correct order of the compound eigenvectors
U2_rec = normalize(compound(U1_rec, k),'norm',2);
W2_rec = normalize(compound(W1_rec, k),'norm',2);
SU_rec = sqrt((CA * CA' * U2_rec)' * U2_rec);
SW_rec = sqrt((CA' * CA * W2_rec)' * W2_rec);

% recompute the tolerance due to the above matrix operation
% C=AB propagates the error \eps in A and B as: n\eps\max(A(:))\max(B(:)),
% when A in mxn and B is nxm matrix. So normalized matrices does not add up
% to scaling error

% tol_new = tol * (max(size(CA)) + max(size(U2_rec)) + max(size(W2_rec))) * max(CA(:))^2;
% tol_new = min(tol_new, 1e-4); % restrict error precision bound
tol_new = 1e-4;

% Keep only the singular vectors corresponding to positive singular values
U2_rec = U2_rec(:,diag(SU_rec) > tol_new);
W2_rec = W2_rec(:,diag(SW_rec) > tol_new);
S2_rec = diag(SU_rec);
S2_rec = diag(S2_rec(S2_rec > tol_new));

% Determine rank(A)
rankA = min(m,n);
try
    for i=min(m,n):-1:1
        if size(S2_rec,1) == nchoosek(i,k)
            rankA = i;
            break
        end
    end
catch
    error('Floating point precision error while checking for zero. Reset tolerance.')
end

% Recover consistent sign of U1_rec and W1_rec by comparing 
% against U2, W2.
[~, perm] = ismembertol(diag(S2_rec),diag(S2(1:size(S2_rec,1),1:size(S2_rec,1))),tol);
[U1_rec, W1_rec] = fixSign(U1_rec(:,1:rankA), U2(:,perm), U2_rec, W1_rec(:,1:rankA), W2(:,perm), W2_rec, k, tol_new);

% Recover singular values, which are positive
S1_rec = solveLinear(diag(S2_rec), rankA, k);
S1_rec = diag(S1_rec);

%----------------Finally, recover the original matrix ---------------
A_rec = U1_rec * S1_rec * W1_rec'
%--------------------------------------------------------------------