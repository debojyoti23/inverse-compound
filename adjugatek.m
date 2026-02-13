function Adj = adjugatek(A,k)
% computes kth order adjugate adj_k(A)
if nargin < 2 
    k = 1;
end
[n,m] = size(A);
if n ~= m
    error("Error: Input must be a square matrix!")
end

idx = nchoosek(1:n,k);
Adj = zeros(size(idx,1));

for i = 1:size(idx,1)
    for j = 1:size(idx,1)
        p = sum(idx(i,:))  + sum(idx(j,:));
        rows = setdiff(1:n,idx(i,:));
        cols = setdiff(1:n,idx(j,:));
        Adj(i,j) = (-1)^p * det(A(cols,rows));
    end
end



