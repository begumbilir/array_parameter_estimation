% Parameters
M = 3; % number of antennas
Delta = 0.5; % distance between elements (in m)
num_sources = 2; % number of sources
theta_true = [-20; 30]; % true directions of arrival in degrees
f = [0.1; 0.3]; %normalized frequencies of the sources
SNR = 20; % signal-to-noise ratio in dB
N = 20; % number of snapshots

% Generate signals
S = exp(1i * 2 * pi * f * (1:N));



% % Steering vectors for the true directions
% A = zeros(M, num_sources);
% element_positions = (0:M-1) * Delta; % positions of the array elements
% for i = 1:num_sources
%     A(:, i) = exp(1i * 2 * pi * element_positions' * sind(theta_true(i)));
% end

% Steering vectors for the true directions
A = zeros(M-1, num_sources);
element_positions = (0:M-2) * Delta; % positions of the array elements
for i = 1:num_sources
    A(:, i) = exp(1i * 2 * pi * element_positions' * sind(theta_true(i)));
end

Phi = diag(exp(1i * 2 * pi * Delta * sind(theta_true)));

% Received signal at the array
X = A * S;
Y = A * Phi * S;

% Add noise
noise_power = 10^(-SNR/10);
noise = sqrt(noise_power/2) * (randn(2*(M-1), N) + 1i * randn(2*(M-1), N));
% noise2 = sqrt(noise_power/2) * (randn(M, N) + 1i * randn(M, N));
% %X = X + noise;
% %Y = Y + noise2;

Z = [X; Y];
Z = Z + noise;


[U_z, e, V_z]  = svd(Z, 'econ');
%e(1:2, 1:2) * V_z(:, 1:2)' == U_z(:, 1:2)' * Z

% U_x = U_z( 1:M , 1:num_sources ); % principal d left singular vectors
% U_y = U_z( M+1:end, 1: num_sources);

U_x = U_z( 1:M-1 , 1:num_sources ); % principal d left singular vectors
U_y = U_z( M:end, 1: num_sources);

%U_x_psuedo_inv = inv(U_x'* U_x)*U_x';

U_x_psuedo_inv2 = (U_x' * U_x) \ U_x'; %pinv(U_x)

U_x_inv_U_y = pinv(U_x) * U_y;
[T_inv, eigvals, T_H] = eig(U_x_inv_U_y);

normalized_eigvals = eigvals ./ abs(eigvals); % map it to the unit circle

angles = angle(diag(normalized_eigvals));

estimated_thetas = asin(angles/(2 * pi * Delta)) * 180/pi;

% Sort and Display Results
estimated_thetas = sort(real(estimated_thetas));

fprintf('True DOAs:       %s\n', mat2str(theta_true));
fprintf('Estimated DOAs:  %s\n', mat2str(estimated_thetas));



% Zero-forcing beamformer
% T = T_inv \ eye(num_sources);
T = inv(T_inv);
W_H = T * inv(U_x); 

S_rec = W_H * X;



W =U_z(:,1:2)* T';
S_rec2 = W' * Z;

% Display rank
S_mtrx = [S; S_rec];
fprintf('Rank of S matrix:       %s\n', mat2str(rank(S_mtrx)));

% Compare with true S
disp('Norm of error between estimated and true S:')
disp(norm(S_rec - S, 'fro') / norm(S, 'fro'))

% Plot true vs estimated real part for 1st source
figure;
plot(real(S(1,:)), 'b', 'LineWidth', 1.5); hold on;
plot(real(S_rec(1,:)), 'r--', 'LineWidth', 1.5);
legend('True S_1', 'Estimated S_1');
xlabel('Snapshot');
ylabel('Real(S)');
title('Source signal recovery with ESPRIT (no noise)');
grid on;
