function [x,hasSolution] = solve_mod2(y,n,k)
% y is a 0-1 vector, 1 for positive, 0 for negative sign
% Solve by Gaussian elimination with respect to xor operations
C = nchoosek(1:n,k);
L = full(sparse(repmat(1:size(C,1),1,k),C(:),1,size(C,1),n));
m = size(L,1);
y = 1 - y; % complement y, abiding by xor logic
L = logical(L);
y = logical(y);
Aug = [L y];

row = 1;

for col = 1:n
    % Find pivot
    pivot = find(Aug(row:m, col), 1) + row - 1;
    if isempty(pivot)
        break
        % error('Matrix does not have independent columns');
    end

    % Swap rows
    Aug([row pivot], :) = Aug([pivot row], :);

    % Eliminate other rows
    for r = 1:m
        if r~= row && Aug(r,col)
            Aug(r,:) = xor(Aug(r,:), Aug(row,:));
        end
    end

    row = row + 1;
end

% Check consistency in remaining rows
for r = n+1:m
    if Aug(r,end)
        hasSolution = false;
        x = [];
        return
    end
end

% Extract solution
x = double(Aug(1:n,end));
hasSolution = true;
x = 1 - x; % complement such that 0 - negative, 1 - positive

end

