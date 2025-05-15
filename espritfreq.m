%% Estimation of frequencies: 1&2. Use the structure of S

% INPUTS:
% X: data matrix
% d: number of sources

% OUTPUTS:
% estimated_freq: the estimated normalized frequencies of the source signals


function estimated_freq = espritfreq(X, d)


[M, N] = size(X);

X_2 = X(:, 1:N-1);
Y = X(:, 2:N);

% Stack them in Z
Z = [X_2; Y];

% Perform singular value decomposition
[U_z, ~, ~]  = svd(Z);

%Take only principal d (number of sources) components
U_x = U_z( 1:M , 1:d ); % principal d left singular vectors
U_y = U_z( M+1:end, 1: d);

eigvals = eig(pinv(U_x) * U_y);

% Normalize eigenvalues to map it to the unit circle
normalized_eigvals = eigvals ./ abs(eigvals); % map it to the unit circle

% Extract the exponent
angles = angle(normalized_eigvals);

% Convert to normalized frequency
estimated_freq = angles/(2 * pi); % T=1
