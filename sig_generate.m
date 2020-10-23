function [Rx,Rs,Rn,Ri,ps,rec_sig] = sig_generate(inr1,inr2,snr,snr_noise,theta_s,theta1,theta2,M)
% signal's model:
% x(k) = s(k) + i(k) + n(k)
%      = s(k)*a + i(k) + n(k)
% inr1: interference-noise-ratio of interference signal 1
% inr2: interference-noise-ratio of interference signal 2
% snr: signal-noise-ratio of target signal
% snr_noise: noise power
% theta_s: target signal's DOA
% theta1: interference signal1's DOA
% theta2: interference signal2's DOA
% M: number of elements

% default parameters:
% f_s:sample frequency
% f:signal's frequency
% T:duration of signal
% c:sound's velocity
% N:snapshot's number
% d:gap between element
f_s = 5000; % sample frequency
f = 1000; % signal frequency
N = 60; % snapshot
T = .05; t = 1/f_s:1/f_s:T;
c = 340;
lamda = c/f;
d = .5*lamda; % default gap between element

% steering vector formulation
ps = exp(-1j*2*pi*d*sind(theta_s)*(0:M-1)'/lamda);
p1 = exp(-1j*2*pi*d*sind(theta1)*(0:M-1)'/lamda);
p2 = exp(-1j*2*pi*d*sind(theta2)*(0:M-1)'/lamda);

% signal generation
tar_sig = wgn(1,length(t), snr+snr_noise); 
inf1 = wgn(1,length(t),inr1+snr_noise); % interference signal 1
inf2 = wgn(1,length(t),inr2+snr_noise); % interference signal 2
noise = wgn(M,length(t),snr_noise);

rec_sig = ps*tar_sig + p1*inf1 + p2*inf2 + noise; % received signal
interference = p1*inf1 + p2*inf2;
sig = ps * tar_sig;

% snapshot estimation of covariance matrix
Rx = rec_sig(:,1:N)*rec_sig(:,1:N)'/N; % total received signal covariance
Rs = sig(:,1:N)*sig(:,1:N)'/N; % signal covariance
Ri = interference(:,1:N)*interference(:,1:N)'/N; % interference covariance
Rn = noise(:,1:N)*noise(:,1:N)'/N; % noise covariance
end


