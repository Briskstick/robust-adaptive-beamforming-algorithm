function optimal_weight = cal_weight(rec_sig,ps,Rx,method,parameters) % optimal weight's calculation
switch lower(method)
    case 'wcp'
        optimal_weight = wcp_method(ps,Rx,parameters);
    case 'mean_cov_based_pro_method'
        optimal_weight = pro_constrained(ps,Rx,parameters);
    case 'yus_method'
        optimal_weight = yus_method(ps,Rx,parameters);
    case 'shrinkage_method'
        optimal_weight = shrinkage_method(rec_sig,ps,Rx,parameters);
    case 'lsmi_method'
        optimal_weight = lsmi(ps,Rx,parameters);
    case 'subspace_method'
        optimal_weight = subspace_method(Rx,parameters);
end
return;