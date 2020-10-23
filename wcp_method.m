function w = wcp_method(ps,Rx,parameters)
    sigma = parameters.sigma; % default algo parameter

    Cov = Rx; % covariance matrix of received signal
    % algo formulation
    M = size(Rx,2);
    cvx_begin
        variable w(M) complex;
        minimize norm(sqrtm(Cov)*w)
        subject to
            sigma * norm(w,2) <= real(w'*ps) - 1;
            imag(w'*ps) == 0;
    cvx_end
    optimal_weight = w;
end
