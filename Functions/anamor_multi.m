function [Z] = anamor_multi(Y, nc)
% Function to transform multiple variables into normal scores (Gaussian values)
% Input:
%   Y: n x (nc + nv) matrix, where nc = # coordinates, nv = # variables
% Output:
%   Z: n x nv matrix of transformed (Gaussian) variables

[n, p] = size(Y);
Z = zeros(n, p);

data = Y;

% Small noise for breaking ties
for i = 1:p
    yi = data(:,i);
    if length(unique(yi)) < n
        noise = randn(n,1) * 1e-9;
    else
        noise = zeros(n,1);
    end

    [ysort, id] = sortrows([yi, noise], [1 2]);
    [~, id2] = sort(id);
    rank = (1:n)';
    zi = norminv(rank / (n + 1), 0, 1);  % Gaussian values
    Z(:,i) = zi(id2);
end
