function x  = gendata_conv(s, P, N, sigma)
%INPUTS: 
% s = QPSK signals (vector) 
% N = length of s (Ns in the book)
% P = oversampling factor
% sigma = standard deviation of Gaussian zero-mean noise 

%OUTPUTS:
% x = sampled received signal (vector form) -> x = [x(0) x(1/P) · · · x(N − 1/P)]T


% Define h(t) over [0, 1] with P samples
L = 4;
t = 0: 1/P : 1-1/P;
h_t = zeros(1, P);

h_t(t >= 0   & t < 0.25)  = 1;
h_t(t >= 0.25 & t < 0.5)  = -1;
h_t(t >= 0.5 & t < 0.75)  = 1;
h_t(t >= 0.75 & t <= 1.0) = -1;


% Oversample s by P (ideally multiple of 4)
s_upsampled = zeros(N*P, 1);
%s_upsampled(1:P:end) = s;
s_upsampled = repelem(s, P);



% Convolve with channel h(t)
x = conv(s_upsampled, h_t, 'same'); %% same returns the central part of the conv which is the same size as s_upsampled

% Add complex Gaussian noise + normalization
noise = sigma/sqrt(2) * (randn(size(x)) + 1j * randn(size(x)));
x = x + noise;






end
