function optimal_weight = yus_method(ps,Rx,parameters)
    M = size(Rx,2); % matrix size
    theta_range1 = parameters.range1;
    theta_range2 = parameters.range2;
    
    C1 = zeros(M,M);
    for iter1 = 1:length(theta_range1)
        a1 = exp(-1j*2*pi*.5*sind(theta_range1(iter1))*(0:M-1)'); % default distance/lamda ratio is .5
        C1 = C1 + inv(a1'*inv(Rx)*a1)*a1*a1';
    end
    
    C2 = zeros(M,M);
    for iter2 = 1:length(theta_range1)
        a2 = exp(-1j*2*pi*.5*sind(theta_range2(iter2))*(0:M-1)'); % default distance/lamda ratio is .5
        C2 = C2 + inv(a2'*inv(Rx)*a2)*a2*a2';
    end
    
    C_recon = C1 + C2; % Capon estimated covariance matrix
    R_inv = inv(C_recon);
    
    cvx_begin 
    variable e(M) complex
    minimize norm(sqrtm(R_inv)*(ps+e))
    subject to
        ps' * e == 0;
        norm(sqrtm(C_recon)*(ps+e)) <= norm(sqrtm(C_recon)*ps);
    cvx_end
    ps_recon = ps + e;
    optimal_weight = R_inv*ps_recon*inv(ps_recon'*R_inv*ps_recon);
end

    
    
        