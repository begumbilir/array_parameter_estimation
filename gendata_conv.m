function x_oversampled  = gendata_conv(s, P, N, sigma)
%INPUTS: 
% s = QPSK signals (vector) 
% N = length of s (Ns in the book)
% P = oversampling factor
% sigma = standard deviation of Gaussian zero-mean noise 

%OUTPUTS:
% x_oversampled = sampled received signal (vector form) -> x = [x(0) x(1/P) · · · x(N − 1/P)]T


% Define h(t) over [0, 1] with L samples
L = 4;         % number of taps (channel memory)

h_sample = [1, -1, 1, -1];  % corresponds to h(t) as piecewise described
% Convolve using vectorized model: x[n] = sum_k h[k] * s[n-k]
x = conv(s, h_sample);

% Add complex Gaussian noise + normalization
noise = sigma/sqrt(2) * (randn(size(x)) + 1j * randn(size(x)));
x = x + noise;

% Oversample s by P 
x_oversampled = upsample(x, P);

end
