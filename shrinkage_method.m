function optimal_weight = shrinkage_method(rec_sig,ps,Rx.parameters)
    M = size(Rx,2); % size of matrix
    v = trace(Rx)/M;
    N = parameters.snap;
    K = parameters.principle; % principle eigenvector
    delta = parameters.Delta;
    
    xx = rec_sig(:,1:N);
    for iter = 1:N
        x_norm(iter) = norm(xx(:,iter),2)^4;
    end
    x_sum = sum(x_norm);
    p = 1/N^2 * x_sum - 1/N * norm(Rx,2)^2;
    arfa_0 = min(v*p/norm(Rx-v*eye(size(Rx,1)),2)^2, v);
    beta_0 = 1 - arfa_0/v;
    R_recon = beta_0*Rx + arfa_0*eye(size(Rx,1));
    
    % reconstruct covariance matrix
    % step1,angle region definition
    angle_sector = parameters.rec_sig_angle;
    outside_scan_angle1 = parameters.outside_scan_angle1;
    outside_scan_angle2 = parameters.outside_scan_angle2;
    angle_sector_disc = linspace(min(angle_sector),max(angle_sector),length(angle_sector));
    outside_angle1_disc = linspace(min(outside_scan_angle1),max(outside_scan_angle1),length(outside_scan_angle1));
    outside_angle2_disc = linspace(min(outside_scan_angle2),max(outside_scan_angle2),length(outside_scan_angle2));
    % matrix reconstruction
    C_recon_sig = exp(-1j*2*pi*.5*(0:M-1)'*sind(angle_sector_disc))*exp(-1j*2*pi*.5*(0:M-1)'*sind(...
        angle_sector_disc)); % default distance/lamda = .5
    C_recon_out1 = exp(-1j*2*pi*.5*(0:M-1)'*sind(outside_angle1_disc))*exp(-1j*2*pi*,5*(0:M-1)'*sind(...
        outside_angle1_disc));
    C_recon_out2 = exp(-1j*2*pi*.5*(0:M-1)'*sind(outside_angle2_disc))*exp(-1j*2*pi*,5*(0:M-1)'*sind(...
        outside_angle2_disc));
    C_recon_out_total = C_recon_out1 + C_recon_out2;
    
%     % eigen-decomposition
%     [U,V] = eig(C_recon_in);
%     U_in = U(:,1:K); % principle eigenvector
%     P = eye(M) - U_in*U_in'; % projection matrix
    
    R_inv = inv(R_recon); % inversion of shrinkage matrix
    cvx_begin 
    variable e(M) complex
    minimize norm(sqrtm(R_inv)*(ps+e))
    subject to 
        norm((ps+e),2) < sqrt(M) + delta; % slack constraint
        ps'*e == 0 ;
        (ps+e)'*P == 0;
        norm(sqrtm(C_recon_out_total)*(ps+e)) <= norm(sqrtm(C_recon_out_total)*ps);
    cvx_end
    ps_recon = ps + e;
    optimal_weight = R_inv*ps_recon*(ps_recon'*R_inv*ps_recon);
end
    
        