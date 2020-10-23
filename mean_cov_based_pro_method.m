function optimal_weight = mean_cov_based_pro(ps,Rx,parameters)

    pro = parameters.pro;
    mean_val = parameters.ste_vec_mean; % its defalut value equal to zero
    var_val = parameters.ste_vec_var;
    M = size(Rx,2);
    Cov = Rx; 
    
    cvx_begin
    variable w(M) complex
    minimize norm(sqrtm(Cov)*w)
    subject to
        norm(var_val*w,2) <= 1/sqrt(-log2(1-pro)) * (real(w'*ps)-1);
    cvx_end
    optimal_weight = w;
end