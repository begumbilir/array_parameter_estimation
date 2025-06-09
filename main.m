%% Signal model - 2. Plot singular values of X

M = 5; % number of antennas
Delta = 0.5; %  normalized distance between antennas 
theta = [-20; 30]; % true directions of arrival in degrees
f = [0.1; 0.3]; % normalized frequencies of the sources w.r.t to the carrier -> for narrowband model
SNR_dB = 20; % signal-to-noise ratio in dB
N = 20; % number of snapshots/samples

% Data generation
[X_gen, A_gen, S_gen] = gendata(M, N, Delta, theta, f, SNR_dB);

% Singular value decomposition
singular_vals = svd(X_gen); 

% Plotting results
figure()
plot(1:length(singular_vals), singular_vals, 'x', 'MarkerSize', 10, 'LineWidth', 4); 
xlabel('Index', 'FontSize', 14);
ylabel('Singular Value', 'FontSize', 14);
title('Singular Values of X', 'FontSize', 16);
set(gca, 'FontSize', 12);
grid on;


%% Estimation of directions  - 2. Test correctness with high SNR

M = 5; % number of antennas
Delta = 0.5; %  normalized distance between antennas 
theta = [-20; 30]; % true directions of arrival in degrees
f = [0.1; 0.3]; % normalized frequencies of the sources w.r.t to the carrier -> for narrowband model
N = 20; % number of snapshots/samples
SNR_dB = 50; % perfect reconstruction, but with order ambiguity

% Data generation with high SNR (noiseless)
[X_gen, A_gen, S_gen] = gendata(M, N, Delta, theta, f, SNR_dB);

% DOA estimation
estimated_thetas = esprit(X_gen, size(f,1));


%% Estimation of frequencies  - 3. Test correctness with high SNR

M = 5; % number of antennas
Delta = 0.5; %  normalized distance between antennas 
theta = [-20; 30]; % true directions of arrival in degrees
f = [0.1; 0.3]; % normalized frequencies of the sources w.r.t to the carrier -> for narrowband model
N = 20; % number of snapshots/samples
SNR_dB = 50; % perfect reconstruction, but with order ambiguity

% Data generation with high SNR (noiseless)
[X_gen, A_gen, S_gen] = gendata(M, N, Delta, theta, f, SNR_dB);

% Frequency estimation
estimated_freq = espritfreq(X_gen, size(f,1));


%% Joint estimation of directions and frequencies - 3. Test correctness with high SNR

M = 5; % number of antennas
Delta = 0.5; %  normalized distance between antennas 
theta = [-20; 30]; % true directions of arrival in degrees
f = [0.1; 0.12]; % normalized frequencies of the sources w.r.t to the carrier -> for narrowband model
N = 20; % number of available snapshots/samples
SNR_dB = 50; % perfect reconstruction, but with order ambiguity
m = 3; % smoothing factor -> maximum the number of antennas
d = size(f,1);
Ns = N-m+1; % number of samples used after smoothing 

% Array response and signal generation 
[~, A_gen, S_gen] = gendata(M, Ns, Delta, theta, f, SNR_dB);

% Define the vector phi and matrix Phi
phi = exp(1j * 2 * pi * f);  
Phi = diag(phi);             % Construct diagonal matrix

% % Create F_phi vector composed of [Phi^0 Phi Phi^2 ... Phi^(m-1)]-> size (m*d) x (d)
% F_phi = zeros(m*d, d); 
% for k = 0:(m-1)
%     F_phi(k*d + (1:d), :) = Phi^k; % Stack each Phi^k
% end

% Received signal at the array by constructing the Khatri Rao structure
% Apply blocks kronecker product
K_kron = zeros(M*m, d);        
for k = 0:(m-1)
    row_start = k*M + 1;
    row_end = (k+1)*M;
    K_kron(row_start:row_end, :) = A_gen * (Phi^k);
end

K = K_kron * S_gen; %S_gen size d x (N-m+1)


% Add noise
SNR_lin = 10^(SNR_dB/10); % linear scale

% Calculate signal power
source_power = norm(S_gen(1,:)); % norm(S_gen(1,:)) = norm(S_gen(2,:)), both sources have save signal power
noise_power_spectral_density = source_power / SNR_lin;  % SNR defined per source

% Each noise sample is complex -> add normalization factor to preserve the scaling 
Noise = sqrt(noise_power_spectral_density / 2) * (randn(size(K)) + 1i * randn(size(K)));

K = K + Noise;


% Perform joint diagonalization

[estimated_thetas_joint, estimated_freq_joint ] = joint(K,d,m); 


%% Comparison - 1.  Make a plot of the estimation performance of the three algorithms
% Parameters
M = 3; % number of antennas
Delta = 0.5; % distance between elements (in m)
d = 2; % number of sources
theta_true = [-20; 30]; % true directions of arrival in degrees
f_true = [0.1; 0.12]; %normalized frequencies of the sources
SNR_db = 0:4:20; % signal-to-noise ratio in dB
N = 20; % number of snapshots
m = M; % smoothing factor in time
Ns = N - m + 1; %number of samples of each source for freq estimation (smoothing in time by factor M)
num_runs = 1000;


% Preallocate
theta_estimates = zeros(length(SNR_db), d, num_runs);
freq_estimates = zeros(length(SNR_db), d, num_runs);
theta_joint_estimates = zeros(length(SNR_db), d, num_runs);
freq_joint_estimates = zeros(length(SNR_db), d, num_runs);

% Define the vector phi and matrix Phi
phi = exp(1j * 2 * pi * f_true);  
Phi = diag(phi);             % Construct diagonal matrix

element_positions = (0:M-1) * Delta; % positions of the array elements
for idx = 1:length(SNR_db)
    SNR = SNR_db(idx);
    fprintf('Running SNR = %d dB\n', SNR);
    for run = 1:num_runs        
        [X, A, S] = gendata(M, N, Delta, theta, f_true, SNR);
        % Data matrix construction for joint estimation
        K_kron = zeros(M*m, d);        
        for k = 0:(m-1)
            row_start = k*M + 1;
            row_end = (k+1)*M;
            K_kron(row_start:row_end, :) = A * (Phi^k);
        end
        
        K = K_kron * S(:, 1:Ns);
        
        % Add noise
        SNR_lin = 10^(SNR_dB/10); % linear scale
        source_power = norm(S_gen(1,:));
        noise_power_spectral_density = source_power / SNR_lin;  % SNR defined per source
        Noise = sqrt(noise_power_spectral_density / 2) * (randn(size(K)) + 1i * randn(size(K)));
        K = K + Noise;

        % Apply ESPRIT
        theta_hat = esprit(X, d);
        % Apply Espritfreq
        freq_hat = espritfreq(X, d);
        % Apply joint
        [theta_joint, freq_joint] = joint(K,d,m);

        % Sort estimates to match true order
        theta_hat = sort(real(theta_hat));
        freq_hat = sort(freq_hat);
        theta_joint = sort(real(theta_joint));
        freq_joint = sort(freq_joint);
        
        theta_estimates(idx, :, run) = theta_hat;
        freq_estimates(idx, :, run) = freq_hat;
        theta_joint_estimates(idx, :, run) = theta_joint;
        freq_joint_estimates(idx, :, run) = freq_joint;
    end
end


% Compute mean and std dev
mean_theta = squeeze(mean(theta_estimates, 3));
std_theta = squeeze(std(theta_estimates, 0, 3));
mean_freq = squeeze(mean(freq_estimates, 3));
std_freq = squeeze(std(freq_estimates, 0, 3));

mean_theta_joint = squeeze(mean(theta_joint_estimates, 3));
std_theta_joint = squeeze(std(theta_joint_estimates, 0, 3));
mean_freq_joint = squeeze(mean(freq_joint_estimates, 3));
std_freq_joint = squeeze(std(freq_joint_estimates, 0, 3));

% Plot results
figure;
% Plot of esprit
subplot(2,1,1);
errorbar(SNR_db, mean_theta(:,1), std_theta(:,1), '-o'); hold on;
errorbar(SNR_db, mean_theta(:,2), std_theta(:,2), '-x');
xlabel('SNR (dB)'); ylabel('Estimated Angle (degrees)');
legend('θ1', 'θ2'); title('Angle Estimation Performance of esprit');

% Plot of espritfreq
subplot(2,1,2);
errorbar(SNR_db, mean_freq(:,1), std_freq(:,1), '-o'); hold on;
errorbar(SNR_db, mean_freq(:,2), std_freq(:,2), '-x');
xlabel('SNR (dB)'); ylabel('Estimated Frequency');
legend('f1', 'f2'); title('Frequency Estimation Performance of espritfreq');


% Plot results
figure;
% Plot of joint theta estimation
subplot(2,1,1);
errorbar(SNR_db, mean_theta_joint(:,1), std_theta(:,1), '-o'); hold on;
errorbar(SNR_db, mean_theta_joint(:,2), std_theta(:,2), '-x');
xlabel('SNR (dB)'); ylabel('Estimated Joint Angle (degrees)');
legend('θ1', 'θ2'); title('Joint Angle Estimation Performance');

% Plot of joint freq. estimation
subplot(2,1,2);
errorbar(SNR_db, mean_freq_joint(:,1), std_freq(:,1), '-o'); hold on;
errorbar(SNR_db, mean_freq_joint(:,2), std_freq(:,2), '-x');
xlabel('SNR (dB)'); ylabel('Estimated Joint Frequency');
legend('f1', 'f2'); title('Joint Frequency Estimation Performance');

%%  Comparison - 2. Compute two zero-forcing beamformers
M = 5; % number of antennas
Delta = 0.5; % distance between elements (in m)
theta = [-20; 30]; % true directions of arrival in degrees
f = [0.1; 0.12]; %normalized frequencies of the sources
SNR = 100; % for SNR greater than 40, no noise is added
N = 20; % number of snapshots


[X_gen, A_gen, S_gen] = gendata(M, N, Delta, theta, f, SNR);

% Apply ESPRIT
theta_estimated = sort(esprit(X_gen, d));
% Apply Espritfreq
freq_estimated = sort(espritfreq(X_gen, d));

% Construct steering matrix using estimated DoAs
A_est_theta = zeros(M, d);
for i = 1:d
    A_est_theta(:, i) = exp(1j * 2 * pi * Delta * (0:M-1)' * sind(theta_estimated(i)));
end

W_H_theta = pinv(A_est_theta);
S_theta_rec = W_H_theta * X_gen;


% Construct steering matrix using estimated freqs
A_est_freq = zeros(M, d);
for i = 1:d
    A_est_freq(:, i) = exp(1j * 2 * pi * freq_estimated(i) * (0:M-1)');
end

W_H_freq = pinv(A_est_freq);
S_freq_rec = W_H_freq * X_gen;

% Display rank
S_mtrx_theta = [S_gen; S_theta_rec];
fprintf('Rank of S matrix:       %s\n', mat2str(rank(S_mtrx_theta)));

S_mtrx_freq = [S_gen; S_freq_rec];
fprintf('Rank of S matrix:       %s\n', mat2str(rank(S_mtrx_freq)));

%% Testing purposes
% Plot of recovered first source using DoA
figure;
plot(real(S_gen(1, :)), 'b-', 'DisplayName', 'Re\{S\_gen(1,:)\}');
hold on;
plot(real(S_theta_rec(1, :)), 'r--', 'DisplayName', 'Re\{S\_theta\_rec(1,:)\}');
hold off;
legend;
xlabel('Index');
ylabel('Real Value');
title('Real Part of First Recovered Source using DoA');

% Plot of recovered second source using DoA
figure;
plot(real(S_gen(2, :)), 'b-', 'DisplayName', 'Re\{S\_gen(1,:)\}');
hold on;
plot(real(S_theta_rec(2, :)), 'r--', 'DisplayName', 'Re\{S\_theta\_rec(1,:)\}');
hold off;
legend;
xlabel('Index');
ylabel('Real Value');
title('Real Part of Second Recovered Source using DoA');

% Plot of recovered first source using freq.
figure;
plot(real(S_gen(1, :)), 'b-', 'DisplayName', 'Re\{S\_gen(1,:)\}');
hold on;
plot(real(S_freq_rec(1, :)), 'r--', 'DisplayName', 'Re\{S\_theta\_rec(1,:)\}');
hold off;
legend;
xlabel('Index');
ylabel('Real Value');
title('Real Part of First Recovered Source using Freq.');

% Plot of recovered first source using freq.
figure;
plot(real(S_gen(2, :)), 'b-', 'DisplayName', 'Re\{S\_gen(1,:)\}');
hold on;
plot(real(S_freq_rec(2, :)), 'r--', 'DisplayName', 'Re\{S\_theta\_rec(1,:)\}');
hold off;
legend;
xlabel('Index');
ylabel('Real Value');
title('Real Part of Second Recovered Source using Freq.');

%%  Comparison - 3. Plotting the spatial response

% Change SNR to 10 db (Other parameters are same)
SNR = 10; 

[X, A, S] = gendata(M, N, Delta, theta, f, SNR);

% Apply ESPRIT
theta_estimated = sort(esprit(X, d));

% Apply Espritfreq
freq_estimated = sort(espritfreq(X, d));


theta_scan = -90:1:90; % Scanning angles for beam pattern

% Construct steering vectors using estimated thetas
a_theta1 = exp(1j * 2 * pi * Delta * (0:M-1)' * sind(theta_estimated(1)));
a_theta2 = exp(1j * 2 * pi * Delta * (0:M-1)' * sind(theta_estimated(2)));
% Compute ZF beamformer
w_H_theta = a_theta1' / norm(a_theta1); %ZF beamformer for source 1
w_H_theta2 = a_theta2' / norm(a_theta2); %ZF beamformer for source 2

% Construct steering vectors using estimated freqs
a_freq1 = exp(1j * 2 * pi * freq_estimated(1) * (0:M-1)');
a_freq2 = exp(1j * 2 * pi * freq_estimated(2) * (0:M-1)');
% Compute ZF beamformer -> norm: Euclidean norm
w_H_freq = a_freq1' / norm(a_freq1); %ZF beamformer for source 1
w_H_freq2 = a_freq2' / norm(a_freq2); %ZF beamformer for source 2

% Calculate zero-forcing beamformer output power for DOA estimation
zf_output_theta = zeros(size(theta_scan));
zf_output_freq  = zeros(size(theta_scan));
zf_output_theta2 = zeros(size(theta_scan)); 
zf_output_freq2  = zeros(size(theta_scan)); 
for i = 1:length(theta_scan)
    a_theta = exp(1i * 2 * pi * Delta * (0:M-1)' * sind(theta_scan(i))); %For uniform linear array
    zf_output_theta(i) = abs(w_H_theta * a_theta); %from DoA, Euclidean norm
    zf_output_freq(i)  = abs(w_H_freq  * a_theta);  % from freq
    zf_output_theta2(i) = abs(w_H_theta2 * a_theta); %from DoA
    zf_output_freq2(i)  = abs(w_H_freq2  * a_theta);  % from freq
end

% Plotting the spatial responses
figure;
plot(theta_scan, zf_output_theta, 'b', 'LineWidth', 1.5); hold on;
plot(theta_scan, zf_output_freq,  'r--', 'LineWidth', 1.5);
xlabel('Angle θ (degrees)');
ylabel('|w^H a(θ)|');
title('Comparison of ZF Beamformer Spatial Responses for First Source');
legend('ZF from DoA estimates', 'ZF from frequency estimates');
grid on;

figure;
plot(theta_scan, zf_output_theta2, 'b', 'LineWidth', 1.5); hold on;
plot(theta_scan, zf_output_freq2,  'r--', 'LineWidth', 1.5);
xlabel('Angle θ (degrees)');
ylabel('|w^H a(θ)|');
title('Comparison of ZF Beamformer Spatial Responses for Second Source');
legend('ZF from DoA estimates', 'ZF from frequency estimates');
grid on;

%% Channel equalization - 1. Signal model

% Parameters
N = 2;          % number of QPSK symbols
P = 2;           % oversampling factor
sigma = 0.1;      % noise std deviation

% Generate QPSK symbols 
constellation = [1+1j, 1-1j, -1+1j, -1-1j]/sqrt(2);
s = constellation(randi(4, N, 1));  % random QPSK symbols

% Generate received signal
x = gendata_conv(s, P, N, sigma);


%% Channel equalization - 2. Signal model
N_s = 500;          % number of QPSK symbols
P = 8;           % oversampling factor
sigma = 0;      % noise std deviation

% Generate QPSK symbols 
constellation = [1+1j, 1-1j, -1+1j, -1-1j]/sqrt(2);
s = constellation(randi(4, 1, N_s));  % random QPSK symbols

% Generate received signal
x = gendata_conv(s, P, N_s, sigma);

X = zeros(2*P,N_s-1);

for k = 1:(N_s-1)
    % Each column contains samples: x(k), x(k+1), ..., x(k + 2P - 1)
    X(:, k) = x(P*(k-1)+1 : P*(k+1));  % MATLAB is 1-indexed
end

rank_X = rank(X)
[~, e] = svd(X);

diag(e)

%% Channel equalization - 1&2. Zero-forcing and Wiener Equalizer
N_s = 500;          % number of QPSK symbols
P = 4;           % oversampling factor
sigma = 0.5;      % noise std deviation

% Generate QPSK symbols 
constellation = [1+1j, 1-1j, -1+1j, -1-1j]/sqrt(2);
s = constellation(randi(4, 1, N_s));  % random QPSK symbols

% Generate received signal
x = gendata_conv(s, P, N_s, sigma);

X = zeros(2*P,N_s-1);

for k = 1:(N_s-1)
    % Each column contains samples: x(k), x(k+1), ..., x(k + 2P - 1)
    X(:, k) = x(P*(k-1)+1 : P*(k+1));  % MATLAB is 1-indexed
end

rank_X_noisy = rank(X)



%% Zero-forcing implementation

% assume we have access to h(t)
if P == 4
    h_sample = [1, 1, -1, 1];  
elseif P == 8
    h_sample = [1, 1, 1, -1, -1, 1, 1, -1]; %for P=8
end

L = 2; %for just span of two symbols (L) 


% Construct H for just span of two symbols (L)
% H = [[h_0 0], [0 h_1]], h_0, 0, h1 are size of 4x1


H = zeros(2*P, L);
H(1:length(h_sample), 1) = h_sample;
H(P+1:end, 2) = h_sample;

% does not matter since symmetric
% H = zeros(2*P, L);
% H(1:length(h_sample), 2) = h_sample;
% H(P+1:end, 1) = h_sample;

% We are trying to build shifted versions of s -> the way we did with H in
% class
W_ZF_H = pinv(H);
S_ZF = W_ZF_H * X; % it should have the size (L)*(N-1)
s_ZF = S_ZF(1,1:end);

check_reconstruction_ZF = [s_ZF; s(1:end-1)]; % s is delayed by L-1

rank_ZF = rank(check_reconstruction_ZF)

[~,sv] = svd(check_reconstruction_ZF);

diag(sv)
%% Wiener Filter Implementation

W_Wiener_H = (inv(H*H' + sigma^2*eye(2*P))*H)';

S_Wiener = W_Wiener_H * X; 
s_Wiener = S_Wiener(1,1:end);

check_reconstruction_Wiener = [s_Wiener; s(1:end-1)];
rank_Wiener = rank(check_reconstruction_Wiener)

[~,sv] = svd(check_reconstruction_Wiener);

diag(sv)


%%

% Plot the constellation for ZF
figure;
hold on; grid on; axis equal;

% Plot the original QPSK constellation
plot(real(s), imag(s), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Original');

% Plot the Zero Forcing equalized constellation
plot(real(s_ZF), imag(s_ZF), 'bx', 'DisplayName', 'ZF Equalized');

% Reference lines
plot([-2 2], [0 0], 'k--');  % x-axis
plot([0 0], [-2 2], 'k--');  % y-axis

% Labels and title
xlabel('In-Phase');
ylabel('Quadrature');
title('QPSK Constellation: Original vs. ZF');
xlim([-2 2]);
ylim([-2 2]);
hold off;

%%

% Plot the constellation for Wiener
figure;
hold on; grid on; axis equal;

% Plot the original QPSK constellation
plot(real(s), imag(s), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Original');

% Plot the Zero Forcing equalized constellation
plot(real(s_Wiener), imag(s_Wiener), 'bx', 'DisplayName', 'ZF Equalized');

% Reference lines
plot([-2 2], [0 0], 'k--');  % x-axis
plot([0 0], [-2 2], 'k--');  % y-axis

% Labels and title
xlabel('In-Phase');
ylabel('Quadrature');
title('QPSK Constellation: Original vs. Wiener');
xlim([-2 2]);
ylim([-2 2]);
hold off;

