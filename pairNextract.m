function [basis_out]  = pairNextract(basis_in, n, d_out)
    % basis_in: each column is linearized basis
    % d_out: dimension of the intersects to be kept
    % basis_out: each column is linearized basis to the intersects given by
    % d_out basis vectors
    d_in = size(basis_in,1)/n;
    basis_out = [];
    for i = 1:size(basis_in,2)-1
        B1 = reshape(basis_in(:,i),n,d_in);
        for j = i+1 : size(basis_in,2)
            B2 = reshape(basis_in(:,j),n,d_in);
            B = intersect(B1,B2);
            if size(B,2) == d_out && d_out > 1
                basis_out = [basis_out B(:)];
            elseif size(B,2) == d_out && d_out == 1
                if rank([basis_out B(:)], 1e-10) > size(basis_out,2)
                    basis_out = [basis_out B(:)];
                end
            end
        end
    end
end

