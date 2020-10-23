function optimal_weight = lsmi(ps,Rx,parameters)
    M = size(Rx,2);
    LNR = parameters.lnr; 
    Loading = LNR * eye(M);
    optimal_weight = inv(Rx+Loading)*ps*inv(ps'*(Rx+Loading)*ps);
end
