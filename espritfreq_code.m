% Parameters
M = 3; % number of antennas
Delta = 0.5; % distance between elements (in m)
num_sources = 2; % number of sources
theta_true = [-20; 30]; % true directions of arrival in degrees
f = [0.1; 0.12]; %normalized frequencies of the sources
SNR = 20; % signal-to-noise ratio in dB
N = 20; % number of snapshots
%m = floor(N/3); %number of rows of A matrix, ideal case
n = N - M + 1; %number of samples of each source

% Generate signals
S = exp(1i * 2 * pi * f * (1:n));



% Steering vectors for the true directions
element_positions = (0:M-1) * Delta; % positions of the array elements
A = zeros(M, num_sources);
for i = 1:num_sources
    A(:, i) = exp(1i * 2 * pi * element_positions' * sind(theta_true(i)));
end


Phi = diag(exp(1i * 2 * pi * f));
% Received signal at the array
X = A * S;
Y = A * Phi * S;

% Add noise
noise_power = 10^(-SNR/10);
noise = sqrt(noise_power/2) * (randn(2*M, n) + 1i * randn(2*M, n));

Z = [X; Y];
Z = Z + noise;

[U_z, ~, ~]  = svd(Z);

U_x = U_z( 1:M , 1:num_sources ); % principal d left singular vectors
U_y = U_z( M+1:end, 1: num_sources);

eigvals = eig(pinv(U_x) * U_y);

normalized_eigvals = eigvals ./ abs(eigvals); % map it to the unit circle

angles = angle(diag(normalized_eigvals));

estimated_freq = angles/(2 * pi); %T=1
