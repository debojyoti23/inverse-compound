function P = wedge2matrix(v, n, k)
    % given the wedge v = v1 \hat v2 \hat...\hat vk, construct a matrix P
    % having nullspace as span(v1,v2,...vk).
    C = nchoosek(1:n, k+1);
    Val = zeros(size(C));
    for i = 1:size(C,1)
        seq = C(i,:);
        for j = 1:length(seq)
            sseq = seq(seq ~= seq(j));
            Val(i,j) = (-1)^(j-1) * v(find_index(sseq, n));
        end
    end
    P = full(sparse(repmat(1:size(C,1),1,k+1), C(:), Val(:), size(C,1), n));
end

