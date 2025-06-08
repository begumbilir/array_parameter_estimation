%% Joint estimation of directions and frequencies:  1&2. Khatri Rao product + time smoothing + diagonalization

% INPUTS:
% X: data matrix
% d: number of sources
% m : time smoothing factor 

% OUTPUTS:
% estimated_joint_thetas: the estimated normalized angles of the source (in degrees)
% estimated_joint_freq: the estimated normalized frequencies of the source signals

function [estimated_joint_thetas, estimated_joint_freq] = joint(K, d, m)

[U, sigma, ~] = svd(K);
U_econ = U(:, 1:d);
M = size(K,1) / m;
% Identity matrices
Im = eye(m);
IM = eye(M);
IM_1 = eye(M-1);
Im_1 = eye(m-1);

% Phi selectors
Jx_phi = kron([Im_1, zeros(m-1,1)], IM);  % size: M*(m-1) x M*m
Jy_phi = kron([zeros(m-1,1), Im_1], IM);  % size: M*(m-1) x M*m

% Theta selectors, selection in every m blocks
Jx_theta = kron(Im, [IM_1, zeros(M-1,1)]);  % size: m*(M-1) x M*m
Jy_theta = kron(Im, [zeros(M-1,1), IM_1]);  % size: m*(M-1) x M*m

% Apply selection matrices
Ux_phi = Jx_phi * U_econ;
Uy_phi = Jy_phi * U_econ;
Ux_theta = Jx_theta * U_econ;
Uy_theta = Jy_theta * U_econ;


M_phi = pinv(Ux_phi)*Uy_phi;
M_theta = pinv(Ux_theta)* Uy_theta;


stacked_M = [M_phi M_theta];

 % Peform EVD to get theta and phi
 [~, Diag] = joint_diag(stacked_M, 1.0e-8);
  
 % Separate Diag to get phi and theta
 eigvals_Phi = Diag(:, 1:d);
 eigvals_Theta = Diag(:, d+1:end);

 % Estimate the frequency 
 angles_f = angle(diag(eigvals_Phi));
 estimated_joint_freq = angles_f/(2 * pi); %T=1

 % Estimate the DOA
 Delta = 0.5;
 angles_theta = angle(diag(eigvals_Theta));
 estimated_joint_thetas = asin(angles_theta/(2 * pi * Delta)) * 180/pi;


end