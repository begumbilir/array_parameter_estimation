function [X,A,S] = gendata(M, N, Delta, theta, f, SNR)
%for SNR greater than 40, no noise is added

element_positions = (0:M-1) * Delta; % Positions of the array elements

% Generate signals
S = exp(1i * 2 * pi * f * (0:N-1));

num_sources = size(f, 1);

% Steering vectors for the true directions
A = zeros(M, num_sources);
for i = 1:num_sources
    A(:, i) = exp(1i * 2 * pi * element_positions' * sind(theta(i)));
end

% Received signal at the array
X = A * S;

% Add noise
noise_power = 10^(-SNR/10);
noise = sqrt(noise_power/2) * (randn(M, N) + 1i * randn(M, N));

if SNR < 40 % for very big SNR's do not add noise
    X = X + noise;
end