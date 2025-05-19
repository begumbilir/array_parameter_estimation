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
f = [0.1; 0.3]; % normalized frequencies of the sources w.r.t to the carrier -> for narrowband model
N = 20; % number of snapshots/samples
SNR_dB = 50; % perfect reconstruction, but with order ambiguity
m = 5; % smoothing factor -> maximum the number of antennas
d = size(f,1);

% Array response and signal generation 
[~, A_gen, S_gen] = gendata(M, N, Delta, theta, f, SNR_dB);

% Define the vectors phi and theta
phi_vec = (exp(1i * 2 * pi * f'));
theta_vec = (exp(1i * 2 * pi * Delta * sind(theta)'));

% Construct the Khatri Rao structure
F = [ones(1, d); phi_vec; theta_vec];

% Received signal at the array
K = khatrirao(F,A_gen) * S_gen;

% Add noise
SNR_lin = 10^(SNR_dB/10); % linear scale

% Signal power = num_sources 
noise_power_spectral_density = (d+1)*m / SNR_lin;  % SNR defined per source

% Each noise sample is complex -> add normalization factor to preserve the scaling 
Noise = sqrt(noise_power_spectral_density / 2) * (randn(size(K)) + 1i * randn(size(K)));

K = K + Noise;


% first try -> use K
%Y_gen = A_gen * Phi;
%Z_gen = A_gen * Theta;
%K = [A_gen *S_gen; Y_gen*S_gen; Z_gen*S_gen];


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
n = N - m + 1; %number of samples of each source for freq estimation (smoothing in time by factor M)
num_runs = 1000;


% Preallocate
theta_estimates = zeros(length(SNR_db), d, num_runs);
freq_estimates = zeros(length(SNR_db), d, num_runs);


element_positions = (0:M-1) * Delta; % positions of the array elements
for idx = 1:length(SNR_db)
    SNR = SNR_db(idx);
    fprintf('Running SNR = %d dB\n', SNR);
    for run = 1:num_runs        
        [X, ~, ~] = gendata(M, N, Delta, theta, f_true, SNR);

        % Apply ESPRIT
        theta_hat = esprit(X, d);
        % Apply Espritfreq
        freq_hat = espritfreq(X(:,1:n), d);

        % Sort estimates to match true order
        theta_hat = sort(real(theta_hat));
        freq_hat = sort(freq_hat);
        
        theta_estimates(idx, :, run) = theta_hat;
        freq_estimates(idx, :, run) = freq_hat;
    end
end


% Compute mean and std dev
mean_theta = squeeze(mean(theta_estimates, 3));
std_theta = squeeze(std(theta_estimates, 0, 3));
mean_freq = squeeze(mean(freq_estimates, 3));
std_freq = squeeze(std(freq_estimates, 0, 3));

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




%%  Comparison - 2. Compute two zero-forcing beamformers
M = 5; % number of antennas
Delta = 0.5; % distance between elements (in m)
theta = [-20; 30]; % true directions of arrival in degrees
f = [0.1; 0.3]; %normalized frequencies of the sources
SNR = 100; % for SNR greater than 40, no noise is added
N = 20; % number of snapshots
m = M; % smoothing factor in time
n = N - m + 1; % number of samples of each source for freq estimation (smoothing in time by factor M)

[X_gen, A_gen, S_gen] = gendata(M, N, Delta, theta, f, SNR);

% Apply ESPRIT
theta_estimated = esprit(X_gen, d);
% Apply Espritfreq
X_gen_freq = X_gen(1:m, 1:n);
freq_estimated = espritfreq(X_gen_freq, d);

% Construct steering matrix using estimated DoAs
A_est_theta = zeros(M, d);
for i = 1:d
    A_est_theta(:, i) = exp(1j * 2 * pi * Delta * (0:M-1)' * sind(theta_estimated(i)));
end

W_H_theta = pinv(A_est_theta);
S_theta_rec = W_H_theta * X_gen;


% Construct steering matrix using estimated freqs
A_est_freq = zeros(m, d);
for i = 1:d
    A_est_freq(:, i) = exp(1j * 2 * pi * freq_estimated(i) * (0:m-1)');
end

W_H_freq = pinv(A_est_freq);
S_freq_rec = W_H_freq * X_gen_freq;

% Display rank
S_mtrx_theta = [S_gen; S_theta_rec];
fprintf('Rank of S matrix:       %s\n', mat2str(rank(S_mtrx_theta)));

S_mtrx_freq = [S_gen; S_theta_rec];
fprintf('Rank of S matrix:       %s\n', mat2str(rank(S_mtrx_freq)));

%%  Comparison - 3. Plotting the spatial response

% Change SNR to 10 db (Other parameters are same)
SNR = 10; 

[X, A, S] = gendata(M, N, Delta, theta, f, SNR);

% Apply ESPRIT
theta_estimated = sort(esprit(X, d));

% Apply Espritfreq
X_freq = X(1:m, 1:n);
freq_estimated = sort(espritfreq(X_freq, d));

% Construct steering matrix using estimated DoAs
A_est_theta = zeros(M, d);
for i = 1:d
    A_est_theta(:, i) = exp(1j * 2 * pi * Delta * (0:M-1)' * sind(theta_estimated(i)));
end

W_H_theta = pinv(A_est_theta);

% Construct steering matrix using estimated freqs
A_est_freq = zeros(m, d);
for i = 1:d
    A_est_freq(:, i) = exp(1j * 2 * pi * freq_estimated(i) * (0:m-1)');
end

W_H_freq = pinv(A_est_freq);

theta_scan = -90:1:90; % Scanning angles for beam pattern
freq_scan = -0.5:0.1:0.5; % Scanning freqs for beam pattern (T = 1)
                        % -pi < 2pi*f*T < pi

% Calculate zero-forcing beamformer output power for DOA estimation
zf_output_theta = zeros(size(theta_scan));
for i = 1:length(theta_scan)
    a_theta = exp(1i * 2 * pi * Delta * (0:M-1)' * sind(theta_scan(i))); %For uniform linear array
    zf_output_theta(i) = norm(W_H_theta * a_theta); %Euclidean norm
end

% Calculate zero-forcing beamformer output power for freq estimation
zf_output_freq = zeros(size(freq_scan));
for i = 1:length(freq_scan) 
    a_freq = exp(1i * 2 * pi * (0:m-1)' * freq_scan(i)); 
    zf_output_freq(i) = norm(W_H_freq * a_freq); %Euclidean norm
end

% Plotting
figure;
plot(theta_scan, 10*log10(zf_output_theta));
xlabel('Angle (degrees)');
ylabel('Zero-forcing Output Power (dB)');
title('Zero-forcing Beamformer DOA Estimation');
grid on;

figure;
plot(freq_scan, 10*log10(zf_output_freq));
xlabel('Frequency (Hz)');
ylabel('Zero-forcing Output Power (dB)');
title('Zero-forcing Beamformer Frequency Estimation');
grid on;