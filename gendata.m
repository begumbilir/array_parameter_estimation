%% Signal model - 1. data, array response and signal matrix generation

% INPUTS:
% M = no of antennas
% N = no of samples
% Delta = normalized antenna spacing: 0.5 - scalar
% theta = direction of the sources in degrees
% f = normalized frequencies of the sources 
% SNR_dB = SNR_dB per source

% OUTPUTS:
% X = Generated data matrix
% A = Generated array response matrix
% S = Generated source matrix


function [X,A,S] = gendata(M, N, Delta, theta, f, SNR_dB)

% Generate source signals
S = exp(1i * 2 * pi * f * (0:N-1)); % num_sources x N samples

num_sources = size(f, 1); % return no of rows -> f is a column vector 

% Steering vectors for the true directions
A = zeros(M, num_sources);

% Positions of the array elements
element_positions = (0:M-1) * Delta;  %n*Delta -> without 2pi

for i = 1:num_sources
    A(:, i) = exp(1i * 2 * pi * element_positions' * sind(theta(i))); % sind: theta in degrees
    % Put the steering vectors together in columns
end

% Received signal at the array
X = A * S;

% Add noise
SNR_lin = 10^(SNR_dB/10); % linear scale

% Signal power = num_sources 
noise_power_spectral_density = sqrt(N) / SNR_lin;  % SNR defined per source

% Each noise sample is complex -> add normalization factor to preserve the scaling 
Noise = sqrt(noise_power_spectral_density / 2) * (randn(M, N) + 1i * randn(M, N));

%for SNR greater than 40, no noise is added for data generation

if SNR_dB < 40 % for very large SNR, do not add noise
    X = X + Noise;
end