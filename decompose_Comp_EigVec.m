function U = decompose_Comp_EigVec(V,n,k)
    % this function decomposes the k-compound V of n dimensional eigenvectors
    % to recover the base Eigenvectors U
    
    tol = 1e-8;

    % Construct a matrix for each vector vi = comp(ui1,...,uik) such
    % that its nullspace is span(ui1,...,uik)
    W_Mat = zeros(nchoosek(n,k+1),n,size(V,2));
    for i=1:size(V,2)
        W_Mat(:,:,i) = wedge2matrix(V(:,i), n, k);
    end
    
    % Recover Eigenvectors
    U = zeros(n,n);
    counter = 0;
    for i=1:size(W_Mat,3)
        N1 = null(W_Mat(:,:,i), tol);
        if isempty(N1)
            % fprintf('null!\n')
            % skip the degenerate eigenvectors
            continue
        end
        for j=i+1:size(W_Mat,3)
            N2 = null(W_Mat(:,:,j), tol);
            N = intersect(N1,N2);
            if ~isempty(N) & rank(N)==1
                if rank([U(:,1:counter) N(:,1)], 1e-10) == counter + 1
                    U(:,counter+1) = N(:,1);
                    counter = counter + 1;
                else
                    continue
                end
            end
        end
    end

end

