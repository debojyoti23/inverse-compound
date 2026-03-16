function [V1_out, W1_out] = fixSign(V1, V2_true, V2, W1, W2_true, W2, k, tol)
    % Recover consistent sign of the singular vectors V1 and W1
    
    n = size(V1,2);
    row_idx_firstnonzero_V2 = ones(1,size(V2,2));
    row_idx_firstnonzero_W2 = ones(1,size(W2,2));
    for j = 1:size(V2,2)
        row_idx_firstnonzero_V2(j) = find(V2(:,j),1);
    end
    for j = 1:size(W2,2)
        row_idx_firstnonzero_W2(j) = find(W2(:,j),1);
    end
    % Compare sign of the first non-zero element of each column
    lin_idx_nzero_V2 = sub2ind(size(V2),row_idx_firstnonzero_V2,1:size(V2,2));
    sign_flip_V2 = abs(V2(lin_idx_nzero_V2) + V2_true(lin_idx_nzero_V2)) < tol;
    sign_flip_V2 = (1 - sign_flip_V2) * 2 - 1; % convert to {+1,-1} array, -1 for flip
    lin_idx_nzero_W2 = sub2ind(size(W2),row_idx_firstnonzero_W2,1:size(W2,2));
    sign_flip_W2 = abs(W2(lin_idx_nzero_W2) + W2_true(lin_idx_nzero_W2)) < tol;
    sign_flip_W2 = (1 - sign_flip_W2) * 2 - 1;
    
    % Propagate the mismatch information from compound to base
    sign_flips = dec2bin(0:2^n-1) - '0';
    sign_flips = -(2 * sign_flips - 1);
    
    V1_out = V1; % Fix V1
    sign_flip_W2 = sign_flip_W2 .* sign_flip_V2;
    for j = 1:size(sign_flips, 1)
        res_flip_W2 = diag(compound(diag(sign_flips(j,:)), k))';
        if all(res_flip_W2 .* sign_flip_W2 == 1)
            % W1 sign-consistent to V1 
            W1_out = W1 * diag(sign_flips(j,:));
            return
        end
    end
end

