function [U, W_Mat] = decompose_Comp_EigVec_1(V,n,k,r)
    % this function recovers U where V = compound(U,k), U has r independent
    % columns representing eigenvectors / singular vectors. Ensure k < r
    
    tol = 1e-6;
    
    % Construct wedge matrix for each column of V
    W_Mat = zeros(nchoosek(n,k+1),n,size(V,2));
    for i=1:size(V,2)
        W_Mat(:,:,i) = wedge2matrix(V(:,i), n, k);
    end


    % Save nullspaces
    U = zeros(n,r);
    kbasis = zeros(n*k,size(V,2));
    for i=1:size(V,2)
        N = null(W_Mat(:,:,i), tol);
        kbasis(:,i) = N(:);
    end
    if k<= ceil(r/2)
        U = pairNextract(kbasis, n, 1);
    else
        overlap_window = 2*k - r;
        basis_overlap = kbasis;
        while overlap_window > 1
            basis_overlap = pairNextract(basis_overlap, n, overlap_window);
            overlap_window = 2 * overlap_window - r;
        end
        U = pairNextract(basis_overlap, n, 1);
    end
end

