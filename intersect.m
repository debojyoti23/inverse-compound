function B_out = intersect(B1,B2)
    % computes the basis of the intersection of two subspaces S1 
    % and S2 with basis B1 and B2 respectively
    N = null([B1 -B2], 1e-10);
    B_out = B1 * N(1:size(B1,2),:);
    B_out = normalize(B_out, 'norm',2);
end