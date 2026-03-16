clear;
clc;

%----------------- Input Rectangular Base Matrix A ----------------
n = 5; m = 5;
k = 1; % Keep k <= min(m,n)/2
A = rand(n,m)
A(:,3) = A(:,4)+ 0.5 * A(:,5);
% A = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 0 0];
% A = diag([4,1,3,2,1]);

% Create a perturbation matrix of same dimension of A
Q = rand(size(A,1),size(A,1));
Q1 = compound(Q,k);


% U0 = rand(n);
% U0 = normalize(U0,'norm',2);
% D0 = diag([5 0 -9 -6 8 1]);
% A = U0 * D0 * inv(U0)
%-----------------------------------------------------------------

[n,m] = size(A);
[U1, S1, W1] = svd(A);

CA = compound(A,k);

% Perturb CA : in order to make non-zero singular values distinct
CA = Q1 * CA;

[U2, S2, W2] = svd(CA);


% tol = 1e-15 * k * max(CA(:));
tol = 1e-10;

U2 = U2(:,diag(S2)>tol);
W2 = W2(:,diag(S2)>tol);
S2 = S2(diag(S2)>tol, diag(S2)>tol);

% compute rank(A)
rankA = min(m,n);
for i=min(m,n):-1:1
    if size(S2,1) == nchoosek(i,k)
        rankA = i;
        break
    end
end

% Ensure k < rank(A)
if k >= rankA
    error('Error! k must be less than rank of the original matrix');
end
% trivial case k = 1
if k == 1
    A_rec = inv(Q1) * CA
    return
end

% Recover original singular vectors upto (possibly inconsistent) 
% sign flips
U1_rec = decompose_Comp_EigVec_1(U2,n,k,rankA);
W1_rec = decompose_Comp_EigVec_1(W2,m,k,rankA);

% Recover a correct order of the compound eigenvectors
U2_rec = normalize(compound(U1_rec, k),'norm',2);
W2_rec = normalize(compound(W1_rec, k),'norm',2);
S2_rec = sqrt((CA * CA' * U2_rec)' * U2_rec);

% Recover consistent sign of U1_rec and W1_rec by comparing 
% against U2, W2.
[~, perm] = ismembertol(diag(S2_rec),diag(S2),tol);
[U1_rec, W1_rec] = fixSign(U1_rec, U2(:,perm), U2_rec, W1_rec, W2(:,perm), W2_rec, k, tol);

% Recover singular values, which are positive
S1_rec = solveLinear(diag(S2_rec), rankA, k);
S1_rec = diag(S1_rec);

%----------------Finally, recover the original matrix ---------------
A_rec = U1_rec * S1_rec * W1_rec';
%--------------------------------------------------------------------

% Remove the perturbation
A_rec = inv(Q) * A_rec