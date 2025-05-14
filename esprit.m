function estimated_thetas = esprit(X,d)
% X: data matrix
% d: number of sources

Delta = 1/2;
[M, ~] = size(X);
X_2 = X(1:M-1, :);
Y = X(2:M, :);

% Stack them in Z
Z = [X_2; Y];

% Perform singular value decomposition
[U_z, ~, ~]  = svd(Z);

%Take only principal d (number of sources) components
U_x = U_z( 1:M-1 , 1:d ); % principal d left singular vectors
U_y = U_z( M:end, 1: d);

U_x_inv_U_y = pinv(U_x) * U_y;
eigvals = eig(U_x_inv_U_y);

% Normalize eigenvalues to map it to the unit circle
normalized_eigvals = eigvals ./ abs(eigvals); 

angles = angle(normalized_eigvals);

estimated_thetas = asin(angles/(2 * pi * Delta)) * 180/pi;

