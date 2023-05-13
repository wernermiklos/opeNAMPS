function outvecs = gram_schmidt_orth(invecs,eps)
%GRAM_SCHMIDT_ORTH Summary of this function goes here
%   Detailed explanation goes here
    if nargin == 1
        eps = 1.e-14;
    end
    overlaps = (invecs(:,1)') * invecs(:,2:end);
    [x,i] = sort(overlaps,'ascend');
    V = invecs(:, [1,i+1]);



    n = size(V, 1);
    k = size(V, 2);
    U = zeros(n, k);
    U(:, 1) = V(:, 1) ;
    for i = 2:k
        U(:, i) = V(:, i);
        for j = 1:i - 1
            U(:, i) = U(:, i) - (U(:, j)'*U(:,i) )/( U(:,j)' * U(:, j)) * U(:, j);
        end
        if sqrt(U(:, i)'*U(:,i)) > eps
            U(:, i) = U(:, i) / sqrt(U(:, i)'*U(:,i));
        end
    end

    outvecs = zeros(n,k);
    for i = 1:k
        if sqrt(U(:,i)' * U(:,i)) > eps
            outvecs(:,i) = U(:,i) / sqrt(U(:,i)' * U(:,i));
        end
    end
    outvecs'*outvecs

end

