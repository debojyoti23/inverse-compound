function Ck = compound(A,k)
    % Creates kth compound of A
    [n,m]=size(A);
    idx_row = nchoosek(1:n,k);
    idx_col = nchoosek(1:m,k);
    
    if isa(A, 'sym')
        Ck = sym(zeros(size(idx_row,1),size(idx_col,1)));
    else
        Ck = zeros(size(idx_row,1),size(idx_col,1));
    end
    
    
    for i = 1:size(idx_row,1)
        for j = 1:size(idx_col,1)
            Ck(i,j) = det(A(idx_row(i,:),idx_col(j,:)));
        end
    end
end

