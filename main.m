function main(inr1,inr2,snr,snr_noise,theta_s,theta1,theta2,method,M)
% Implements the robust adaptive beamforming algorithm
% signal's model:
%           x(k) = s(k) + i(k) + n(k)
%               = s(k)*a + i(k) + n(k)
%           inr1: interference-noise-ratio of interference signal 1
%           inr2: interference-noise-ratio of interference signal 2
%           snr: signal-noise-ratio of target signal
%           snr_noise: noise power
%           theta_s: target signal's DOA
%           theta1: interference signal1's DOA
%           theta2: interference signal2's DOA
%           M: number of elements
% method - robust beamforming algorithms
%          'WCP' = worst-case performance algorithm[1]
%          'MCBP' = given error steering vector,probability constrained
%          algorithm[2]
%          'QCQP' = covariance matrix reconstrunction based on Capon
%          estimation[3]
%          'Shrinkage' = covariance matrix reconstrution based on MMSE
%          estimation of the received covariance matrix[4]
% Example call: rab_algo(30,30,10,0,10,40,-30,'WCP',10)

if nargin<9
    fprintf('Usage: rab_algo(inr1,inr2,snr,snr_noise,theta_s,theta1,theta2,method) \n');
    fprintf('where''method'' is of the following:'\n);
    fprintf('WCP = worst-case performance'\n);
    fprintf('MCBP = probability constrianed based on given mean and variance'\n);
    fprintf('QCQP = Yu_method which based on Capon estimation and cvx program'\n);
    fprintf('Shrinkage = MMSE based covariance matrix reconstruction'\n);
    fprintf('LSMI = loading sample matrix inversion method'\n);
    fprintf('Subspace = subspace based steering vector estimation'\n);
    
    return;
end

% signal generation
[Rx,Rs,Rn,Ri,ps,rec_sig] = sig_generate(inr1,inr2,snr,snr_noise,theta_s,theta1,theta2,M); 

% optimal weight vector calculation
parameters = initialise_parameters(theta_s,theta1,theta2,Rx,method); % initialization
w_opt = cal_weight(rec_sig,ps,Rx,method,parameters); % calculate the optimal weight


% plotation
myplot(Rx,w_opt);
snr_output = snr_calculation(Rs,Ri,Rn,w_opt);
% caution
%！！yus_method-work
%！！wcp-work
%！！lsmi-work

end






