% Script to generate one channel matrix and estimate its parameters using TSDCE

%% Channel parameters

L  = 3;                        % number of paths 
nt = 16;                       % number of transmit antennas  
nr = 16;                       % mumber of receive antennas   
SNR = 10;                      % SNR value [dB]

%% Initialize Method

% Codebook initialization for TSDCE
P = 32;                        % Number of quantized angles in transmission
Q = 32;                        % Number of quantized angles in reception

% steering vector function
ev = @(nt,angle) (1/sqrt(nt))*exp(-1i*pi*cos(angle)*(0:nt-1)');

% P beamforming vectors covering the range [0,pi];
cosp = -(1/pi)*angle(exp(-1i*(2*pi/P)*(0:P-1))); 
phi_p = acos(cosp);
F = zeros(nt,P);               % Beamforming weight vector
for p = 1:P
    F(:,p) = ev(nt,phi_p(p));
end

% Q beamforming vectors covering the range [0,pi];
cosq = (1/pi)*angle(exp(-1i*(2*pi/Q)*(0:Q-1)));
psi_q = acos(cosq);
W = zeros(nr,Q);               % Combiner weight vector 
for q = 1:Q
    W(:,q) = ev(nr,psi_q(q));   
end

K = 2;                         % Refinement parameter

%% Generate observation

rho = 1;                       % Transmit power
varn = rho*10^(-SNR/10);       % Noise variance 
Sigma_n = (varn)*eye(Q);       % Noise covariance matrix at the receiver

% Alpha coefficient
alpha_l = sqrt(1/L)*sqrt(1/2)*randn(L,1) + 1i*sqrt(1/L)*sqrt(1/2)*randn(L,1);
% AoD (Angle of Departure) for each path
phi_l = pi*rand(L,1);
% AoA (Angle of Arrival) for each path
psi_l = pi*rand(L,1);
% angle_vector
angle_v = [phi_l; psi_l];

% Note:
% Structure of angle_v is:
% [phi_1; phi_2; ...; phi_L; psi_1; psi_2; ...; psi_L]
% Structure of alpha_l is:
% [alpha_1; alpha_2; ...; alpha_L]

% Generate channel 
 Hl = zeros(nr,nt,L);
 for l = 1:L
     Hl(:,:,l) = alpha_l(l)*ev(nr,angle_v(l+L))*ev(nt,angle_v(l))';
 end
 Hl = sqrt(nt*nr)*Hl;
 H  = sum(Hl,3);
 
 % Generate noise
 N = (1/sqrt(2))*(mvnrnd(zeros(Q,P),Sigma_n) + 1i*mvnrnd(zeros(Q,P),Sigma_n));
 
 % Noiseless measurement
 G = sqrt(rho)*W'*H*F;
 
 % Noisy measurement
 Y = G + N;
 figure(1)
 imagesc(0:P-1,0:Q-1,abs(Y))
 xlabel('p')
 ylabel('q')
 title('Observation matrix Y')
  
 %================================
 % Call TSDCE method
 %================================
 
 [angles, alphas, Hl_est] = runtsdce(Y,nt,nr,L,K,rho);

 H_est = sum(Hl_est,3);
  
 % Performance metric: Normalized squared error
 NSE = (norm(H_est-H,'fro')/norm(H,'fro'))^2
                
 
