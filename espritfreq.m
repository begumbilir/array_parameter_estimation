function estimated_freq = espritfreq(X, d)
% X: data matrix
% d: number of sources

[M, n] = size(X);

X_2 = X(:, 1:n-1);
Y = X(:, 2:n);

Z = [X_2; Y];


[U_z, ~, ~]  = svd(Z);

U_x = U_z( 1:M , 1:d ); % principal d left singular vectors
U_y = U_z( M+1:end, 1: d);

eigvals = eig(pinv(U_x) * U_y);

normalized_eigvals = eigvals ./ abs(eigvals); % map it to the unit circle

angles = angle(normalized_eigvals);

estimated_freq = angles/(2 * pi); %T=1
