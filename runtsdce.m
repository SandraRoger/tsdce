function [angle_est, alpha_est, Hl_est] = runtsdce(Y,nt,nr,L,K,rho)

% Implements Algorithm 1 of manuscript S. Roger, M. Cobos, C. Botella-Mascarell 
% and G. Fodor, "Fast Channel Estimation in the Transformed Spatial Domain 
% for Analog Millimeter Wave Systems", IEEE Trans. on Wireless Comm., 2021"

% Allocate output vectors
angle_est  = zeros(2*L,1);
alpha_est = zeros(L,1);
Hl_est = zeros(nt,nr,L);

% steering vector function
ev = @(nt,angle) (1/sqrt(nt))*exp(-1i*pi*cos(angle)*(0:nt-1)');

% Get codebook size from observation matrix
[Q,P] = size(Y);

% Transformed spatial domain observation
D = ifft2(Y);
figure(2)
imagesc(0:P-1,0:Q-1,abs(D))
xlabel('p')
ylabel('q')
title('Transformed spatial domain observation D')

% Get informative part of D
Dc = D(1:nt,1:nr);

% Initialize variables
Ct_est = zeros(nt,nr,L);
Lset = 1:L;

for k = 1:K
    for l = 1:L
        %Reconstruct and remove previously estimated paths (if any)
        Previous_paths = sqrt(rho)*sum(Ct_est(:,:,Lset~=l),3);
        Dc_prime = Dc - Previous_paths;
        
        % SVD Path Separation
        if k == 1 && l < L
            [U,S,V] = svd(Dc_prime);
            Dc_tilde = S(1,1)*U(:,1)*V(:,1)';
        else
            Dc_tilde = Dc_prime;
        end
        
         % 2D Autocorrelation
        kappa = ((nt:-1:1)'*(nr:-1:1));
        R = xcorr2_fft(Dc_tilde);
        R = R(1:nr, 1:nt);
        R = R./kappa;
        
        % Frequency estimation
        c_psi_est = R(:,1);
        c_phi_est = R(1,:).';
        i_v = 0:nr-1;
        wi = ((nr+1)*(nr-i_v))./(i_v +1);
      
        [omega_psi_est, omega_phi_est] = wlsfreqest(c_psi_est,c_phi_est,wi);

        % Complex coefficient estimation
        Kp = (0.25*nr*(nr+1)*nt*(nt+1)-nt*nr)^-1;
        A2_est = Kp*(sum(abs(kappa(2:end).*R(2:end))));
        alpha_mag = sqrt(nt*nr/rho)*sqrt(A2_est);
        Sphase_est = exp(1i*omega_psi_est*(0:1:nr-1)).'*exp(1i*omega_phi_est*(0:1:nt-1));
        alpha_phase = angle(mean(Dc_tilde(:).*conj(Sphase_est(:))));
        
        alpha_est(l) = alpha_mag.*exp(1i*alpha_phase);
        angle_est(l)   = acos(omega_phi_est/pi);
        angle_est(l+L) = acos(-omega_psi_est/pi);
        
        % Recover channel from angles
        Hn = sqrt(nt*nr)*alpha_est(l)*ev(nr,angle_est(l+L))*ev(nt,angle_est(l))';
        Hl_est(:,:,l) = Hn;
                
        Al_est = alpha_est(l)/sqrt(nt*nr);
        Ct_est(:,:,l) = Al_est*Sphase_est;
        
    end
end
            
