function x  = gendata_conv(s, P, N, sigma)
%INPUTS: 
% s = QPSK signals (vector) 
% N = length of s (Ns in the book)
% P = oversampling factor
% sigma = standard deviation of Gaussian zero-mean noise 

%OUTPUTS:
% x = sampled received signal (vector form) -> x = [x(0) x(1/P) · · · x(N − 1/P)]T

h_sample = [1, -1, 1, -1];  % corresponds to h(t) as piecewise described

x = zeros(1,N*P);
% Create a time vector t from 0 to (N - 1/P) with increments of 1/P
% This ensures t has N*P samples, evenly spaced
t = 0 : 1/P : N - 1/P;

% Each row h(k+1, :) will correspond to the kth shifted channel function
h = zeros(N, length(t));  % h will be N rows (k) by length(t) columns

for k = 0:N-1
    tau = t - k;
    
    % Define a piecewise function h_k(tau):
    % Set values to +1 in the intervals [0, 0.25] and (0.5, 0.75]
    % Set values to -1 in the intervals (0.25, 0.5] and (0.75, 1)
    h(k+1, (tau >= 0   & tau <= 0.25) | (tau > 0.5  & tau <= 0.75))  = 1;
    h(k+1, (tau > 0.25 & tau <= 0.5)  | (tau > 0.75 & tau < 1.0))  = -1;
    
    % Form the output signal x by summing scaled versions of h_k(t)
    % Each row h(k+1, :) is weighted by the corresponding s(k+1) 
    x = x + s(k+1) * h(k+1, :); % Multiply the row by s(k+1)
end


% Add complex Gaussian noise + normalization
noise = sigma/sqrt(2) * (randn(size(x)) + 1j * randn(size(x)));
x = x + noise;

end
