function myplot(Rx,optimal_weight)
    M = size(Rx,2); % matrix size
    theta = -90:1:90; % scan angle
    p = exp(-1j*2*pi*.5*(0:M-1)'*sind(theta));
    y = optimal_weight'*p;
    figure;
    plot(theta,20*log10(abs(y)/max(abs(y))));
end