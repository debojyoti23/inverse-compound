function [V1_out, W1_out] = fixSign(V1, V2_true, V2, W1, W2_true, W2, k, tol)
    % Recover consistent sign of the singular vectors V1 by comparing 
    % their kth compound V2 against the certificate V2_true
    
    n = size(V1,2);
    sign_flip_comp_V = abs(V2(1,:) + V2_true(1,:)) < tol;
    sign_flip_comp_V = (1 - sign_flip_comp_V) * 2 - 1; % convert to {+1,-1} array, -1 for flip
    sign_flip_comp_W = abs(W2(1,:) + W2_true(1,:)) < tol;
    sign_flip_comp_W = (1 - sign_flip_comp_W) * 2 - 1;
    % Propagate the mismatch information from compound to base
    sign_flips = dec2bin(0:2^n-1) - '0';
    sign_flips = -(2 * sign_flips - 1);
    % res_flips = zeros(size(sign_flips,1),nchoosek(n,k));
    for i=1:size(sign_flips,1)
        res_flip_1 = diag(compound(diag(sign_flips(i,:)), k))';
        idx_pairflip = find(res_flip_1 .* sign_flip_comp_V < 0);
        temp_flip_W = sign_flip_comp_W;
        temp_flip_W(idx_pairflip) = temp_flip_W(idx_pairflip) * -1;
        for j = 1:size(sign_flips, 1)
            res_flip_2 = diag(compound(diag(sign_flips(j,:)), k))';
            if all(res_flip_2 .* temp_flip_W == 1)
                V1_out = V1 * diag(sign_flips(i,:));
                W1_out = W1 * diag(sign_flips(j,:));
                % fprintf('Iteration count %d %d',i,j)
                return
            end
        end
    end
end

