function optimal_weight = subspace_method(Rx,parameters)
% 代码可能存在问题
    theta_s_range = parameters.source_range;
    theta1_range = parameters.interference1_range;
    theta2_range = parameters.interference2_range;
    L = parameters.prin; % principle component
    thr = parameters.threshold; % default threshold
    M = size(Rx,1);
    
    % signal's subspace
    [U1,V1] = eig(Rx);
    [V1_sort,index] = sort(diag(V1),'descend');
    U1_sort = U1(:,index);
    U_sig_noise = U1_sort(:,1:L); % principle component
    
    % noise's variance calculation
    eig_noise = V1_sort(L+1:M); % smallest eigenvalue
    eig_noise_sum = sum(eig_noise);
    delta_noise = 1/(M-L) * eig_noise_sum;
    delta_noise = delta_noise * eye(M);
    
    % interference matrix recosntruction
    C_int1 = zeros(M,M);
    for num1 = 1:length(theta1_range)
        a_int1 = exp(-1j*2*pi*.5*sind(theta1_range(num1))*(0:M-1)'); % default distance/lamda ratio is .5
        C_int1 = C_int1 + a_int1*a_int1'*inv(a_int1'*inv(Rx)*a_int1);
    end
    % select M largest eigenvalue
    [U_int1, V_int1] = eig(C_int1);
    [V_int1_sort,index1] = sort(diag(V_int1),'descend');
    U_int1_sort = U_int1(:,index1);
    V_int1_sum = sum(V_int1_sort);
    V_int1_thr = 0; int1_index = 0;
    for thr1 = 1:M
        V_int1_thr = V_int1_thr + V_int1_sort(thr1);
        V_int1_bound = V_int1_thr/V_int1_sum;
        if(V_int1_bound >= thr)
            if(thr1 == 1)
                int1_index = thr1;
                break;
            else
                int1_index = thr1 - 1;
                break;
            end
        else
            int1_index = thr1;
        end
    end
    U_int1_actual = U_int1_sort(:,1:int1_index);
    
    C_int2 = zeros(M,M); % interference 2 covariance matrix
    for num2 = 1:length(theta2_range)
        a_int2 = exp(-1i*2*pi*.5*sind(theta2_range(num2))*(0:M-1)');
        C_int2 = C_int2 + a_int2*a_int2'*inv(a_int2'*inv(Rx)*a_int2);
    end
    [U_int2, V_int2] = eig(C_int2);
    [V_int2_sort, index2] = sort(diag(V_int2),'descend');
    U_int2_sort = U_int2(:,index2);
    V_int2_sum = sum(V_int2_sort);
    V_int2_thr = 0; int2_index = 0;
    for thr2 = 1:M
        V_int2_thr = V_int2_thr + V_int2_sort(thr2);
        V_int2_bound = V_int2_thr / V_int2_sum;
        if(V_int2_bound >= thr)
            if(thr2 == 1)
                int2_index = 1;
                break;
            else
                int2_index = thr2 - 1;
                break;
            end
        else
            int2_index = thr2;
        end
    end
    U_int2_actual = U_int2_sort(:,1:int2_index);
    
    P_E = U_sig_noise * U_sig_noise';
    P_B1 = U_int1_actual * U_int1_actual';
    P_int1 = P_B1 * P_E;
    [U_P1, V_P1] = eig(P_int1);
    [V_P1_sort, index_P1] = sort(diag(V_P1), 'descend');
    U_P1_sort = U_P1(:,index_P1);
    U_pri1 = U_P1_sort(:,1); % principle eigenvector correspond to the largest eigenvalue
    a_int1 = sqrt(M) * U_pri1;
    
    P_B2 = U_int2_actual * U_int2_actual';
    P_int2 = P_B2 * P_E;
    [U_P2, V_P2] = eig(P_int2);
    [V_P2_sort, index_P2] = sort(diag(V_P2), 'descend');
    U_P2_sort = U_P2(:,index_P2);
    U_pri2 = U_P2_sort(:,1); % principle eigenvector correspond to the largest eigenvalue
    a_int2 = sqrt(M) * U_pri2;
    
    delta_int1 = 1*inv(a_int1'*inv(Rx)*a_int1);
    delta_int2 = 1*inv(a_int2'*inv(Rx)*a_int2);
    R_int = delta_int1*a_int1*a_int1' + delta_int2*a_int2*a_int2'; % reconstructed interference covariance matrix
    R_in = R_int + delta_noise; % reconstructed interference plus noise matrix
    
    C_s = zeros(M,M);
    for num3 = 1:length(theta_s_range)
        a_s = exp(-1i*2*pi*.5*sind(theta_s_range(num3))*(0:M-1)');
        C_s = C_s + a_s*a_s'*inv(a_s'*inv(Rx)*a_s);
    end
    [U_s, V_s] = eig(C_s);
    [V_s_sort, index_s] = sort(diag(V_s), 'descend');
    U_s_sort = U_s(:,index_s);
    U_sig = U_s_sort(:,1); % eigenvector correspond to the largest eigenvalue
    a_s = sqrt(M) * U_sig; % reconstructed signal steering vector
    
    optimal_weight = inv(R_in)*a_s*inv(a_s'*R_in*a_s); % optimal weighted vector
end