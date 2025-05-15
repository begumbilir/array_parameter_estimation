%% Estimation of directions - 1. Implement the ESPRIT algorithm

% INPUTS:
% X: data matrix
% d: number of sources

% OUTPUT:
% estimated_thetas: estimated DOA angle in degrees


function estimated_thetas = esprit(X,d)

Delta = 1/2;
[M, ~] = size(X);
X_2 = X(1:M-1, :); %% slicing the initial data matrix 
Y = X(2:M, :);

% Stack them in Z
Z = [X_2; Y];

% Perform singular value decomposition
[U_z, ~, ~]  = svd(Z);

%Take only principal d (number of sources) components
U_x = U_z( 1:M-1 , 1:d ); % principal d left singular vectors
U_y = U_z( M:end, 1:d);

U_x_inv_U_y = pinv(U_x) * U_y;
eigvals = eig(U_x_inv_U_y);

% Normalize eigenvalues to map it to the unit circle
normalized_eigvals = eigvals ./ abs(eigvals); 

% Extract the exponent
phase_shifts = angle(normalized_eigvals); 

% Convert to degrees
estimated_thetas = asin(phase_shifts/(2 * pi * Delta)) * 180/pi; 

