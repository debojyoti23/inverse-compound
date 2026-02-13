function indx = find_index(u, n)
    % Find the index of the k-length sequence u in the n choose k list
    k = length(u);
    indx = 0;
    n_i = 1;
    u_i = 1;
    while u_i <= k
        if n_i < u(u_i)
            indx = indx + nchoosek(n-n_i,k-u_i);
            n_i = n_i + 1;
        else
            n_i = n_i + 1;
            u_i = u_i + 1;
        end
    end
    indx = indx + 1;
end

