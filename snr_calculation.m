function snr_output = snr_calculation(Rs,Ri,Rn,w_opt)
    R_in = Rn + Ri;
    snr_output = 10*log10(abs(w_opt'*Rs*w_opt*inv(w_opt'*R_in*w_opt)));
    fprintf('snr_output = %f\n', snr_output);
end