 function C = khatrirao(A, B)

%   Computes the Khatri-Rao product of matrices A and B.
%   C = KHATRIRAO(A, B) returns the column-wise Kronecker product.
%   A and B must have the same number of columns.

    % Check column compatibility
    if size(A, 2) ~= size(B, 2)
        error('A and B must have the same number of columns.');
    end

    % Preallocate output
    m = size(A, 1);
    p = size(B, 1);
    n = size(A, 2);
    C = zeros(m * p, n);

    % Compute column-wise Kronecker product
    for i = 1:n
        C(:, i) = kron(A(:, i), B(:, i));
    end
end