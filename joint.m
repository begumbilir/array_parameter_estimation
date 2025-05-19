%% Joint estimation of directions and frequencies:  1&2. Khatri Rao product + time smoothing + diagonalization

% INPUTS:
% X: data matrix
% d: number of sources
% m : time smoothing factor 

% OUTPUTS:
% estimated_joint_thetas: the estimated normalized angles of the source (in degrees)
% estimated_joint_freq: the estimated normalized frequencies of the source signals

function [estimated_joint_thetas, estimated_joint_freq] = joint(K, d, m)

% assume P = 1 ( no oversampling -> we use normalized frequencies)
% assume B  full ones, so no attenuation
 P = 1; 

 [U, sigma, ~] = svd(K);

 U_x = U( 1:m , 1:d ); % principal d left singular vectors
 U_y = U( m+1:2*m, 1: d);
 U_z = U( 2*m+1:end, 1: d);

 % Assume m >= num sources (which is the case) for psuedo inv of U_x
 M_y = pinv(U_x) * U_y;
 M_z = pinv(U_x) * U_z;

 stacked_M = [M_y M_z];

 % Peform EVD to get theta and phi
 [~, Diag] = joint_diag(stacked_M, 1.0e-8);
  
 % Separate Diag to get phi and theta
 eigvals_Phi = Diag(:, 1:d);
 eigvals_Theta = Diag(:, d+1:end);

 % Estimate the frequency 
 angles_f = angle(diag(eigvals_Phi));
 estimated_joint_freq = angles_f/(2 * pi/P); %T=1

 % Estimate the DOA

 Delta = 0.5;
 angles_theta = angle(diag(eigvals_Theta));
 estimated_joint_thetas = asin(angles_theta/(2 * pi * Delta)) * 180/pi;


end