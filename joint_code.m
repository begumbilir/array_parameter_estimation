% Parameters
M = 3; % number of antennas
Delta = 0.5; % distance between elements (in m)
num_sources = 2; % number of sources
theta_true = [-20; 30]; % true directions of arrival in degrees
f = [0.1; 0.3]; %normalized frequencies of the sources
SNR = 20; % signal-to-noise ratio in dB
N = 20; % number of snapshots
%m = floor(N/3); %number of rows of A matrix, ideal case (smoothing factor
%in time)
n = N - M + 1; %number of samples of each source
%n = N - m + 1; %number of samples of each source


% Generate signals
S = zeros(num_sources, n);
sample_idx = 1:n;
for i = 1:num_sources
    S(i, :) = exp(1i * 2 * pi * f(i) * sample_idx);
end


% Steering vectors for the true directions
element_positions = (0:M-1) * Delta; % positions of the array elements
A = zeros(M, num_sources);
for i = 1:num_sources
    A(:, i) = exp(1i * 2 * pi * element_positions' * sind(theta_true(i)));
end


Phi = exp(1i * 2 * pi * f');
Theta = exp(1i * 2 * pi * Delta * sind(theta_true)');

F = [ones(1, num_sources); Phi; Theta];
% Received signal at the array
K = khatrirao(F,A) * S;

% Add noise
noise_power = 10^(-SNR/10);
noise = sqrt(noise_power/2) * (randn(size(K)) + 1i * randn(size(K)));


K = K + noise;


[U, ~, ~]  = svd(K);

U_x = U( 1:M , 1:num_sources ); % principal d left singular vectors
U_y = U( M+1:2*M, 1: num_sources);
U_z = U( 2*M+1:end, 1: num_sources);

% Assume M >= num sources (which is the case) for psuedo inv of U_x
M_y = pinv(U_x) * U_y;
M_z = pinv(U_x) * U_z;

[T_inv, eigvals_Phi] = eig(M_y);

%normalized_eigvals = eigvals_Phi ./ abs(eigvals_Phi); % map it to the unit circle

% Estimate freq
angles_f = angle(diag(eigvals_Phi));
estimated_freq = angles_f/(2 * pi); %T=1

% Estimate theta (DoA's)
eigvals_Theta = inv(T_inv) * M_z * T_inv; %Scaling of T has no effect, since signal power is 1
angles = angle(diag(eigvals_Theta));
estimated_thetas = asin(angles/(2 * pi * Delta)) * 180/pi;

% Sort and Display Results
estimated_freqs = sort(real(estimated_freq));
estimated_thetas = sort(real(estimated_thetas));

fprintf('True DOAs:       %s\n', mat2str(f));
fprintf('Estimated DOAs:  %s\n', mat2str(estimated_freq));
fprintf('True DOAs:       %s\n', mat2str(theta_true));
fprintf('Estimated DOAs:  %s\n', mat2str(estimated_thetas));
