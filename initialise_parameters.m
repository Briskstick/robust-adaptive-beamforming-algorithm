function parameters = initialise_parameters(theta_s,theta1,theta2,Rx,method) % number of elements,covariance matrix
mat_size = size(Rx,1); % covariance matrix size
switch lower(method)
    case 'wcp'
        Sigma = .3;
        wcp_param = Sigma*mat_size;
        parameters = struct('sigma',wcp_param); % default wcp parameter
%     case 'pro_constrained'
%         pro = .95; % default probability constrained parameter
%         parameters = struct('probability',pro); 

    case 'mean_cov_based_pro_method'
        ste_vec_mean = 0 * eye(mat_size);
        ste_vec_var = sqrt(.3*(mat_size)/mat_size*eye(mat_size));
        pro = .95;
        parameters = struct('probability',pro,'mean',ste_vec_mean,'var',ste_vec_var);
        
    case 'yus_method'
        theta_1_disc_range = theta1-5:.1:theta1+5; % .1 stands for sample gap
        theta_2_disc_range = theta2-5:.1:theta2+5; 
        parameters = struct('range1',theta_1_disc_range,'range2',theta_2_disc_range);
        
    case 'shrinkage_method'
        N = 60; % numbers of snapshot
        K = 6; % principle eigenvector
        Delta = .1;
        rec_sig_angle = [theta_s - 5, theta_s + 5];
        % scan angle outside of received signal region
        outside_scan_angle1 = [-90, min(rec_sig_angle)];
        outside_scan_angle2 = [max(rec_sig_angle), 90];
        parameters = struct('snap',N,'principle',K,'delta',Delta,'sig_range',rec_sig_angle,'out_range1',outside_scan_angle1,...
            'out_range2',outside_scan_angle2);
        
    case 'lsmi_method'
        LNR = 10; % default loading noise level-10dB
        parameters = struct('lnr',LNR);
        
    case 'subspace_method'
        theta_s_range = theta_s-8:.5:theta_s+8;
        theta_1_range = theta1-8:.5:theta1+8;
        theta_2_range = theta2-8:.5:theta2+8;
        L = 3; % principle eigenvalue
        thr = .9; % threhold
        parameters = struct('source_range',theta_s_range,'interference1_range',theta_1_range,...
            'interference2_range',theta_2_range,'prin',L,'threshold',thr);
    otherwise
        error('Method not implemented. Check spelling.');
end
