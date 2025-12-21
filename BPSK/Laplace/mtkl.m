clear;
clc;
num_samples = 1e6;
SNR = 1:2:7;
data = zeros(length(SNR),4);
for j = 1:length(SNR)
    S = SNR(j);
    P = 10^(S/10);
    clear points;
    UU = rand(1,num_samples) - 0.5;       % U ~ Uniform(-0.5,0.5)
    sigma = 1; 
    b = sigma/sqrt(2);             
    Z = -b * sign(UU) .* log(1 - 2*abs(UU));
    points = zeros(num_samples,3);
    cdf = cumsum([0.5,0.5]);
    U = rand(1, num_samples);
    idx = arrayfun(@(u) find(u < cdf, 1), U);
    angle = [+1, -1];
    points(:,1) = angle(idx).*sqrt(P)+Z;
    points(:,2) = sqrt(P)+Z;
    points(:,3) = -sqrt(P)+Z;
    % file_path = fullfile(pwd, sprintf('points_BPSK_SNR%d.mat',S));
    % save(file_path, 'points');
    Y = points(1:end, 1);
    num_samples = length(Y);
    q_0 = 1/sqrt(2)*exp(-sqrt(2).*abs(Y-sqrt(P)));
    q_1 = 1/sqrt(2)*exp(-sqrt(2).*abs(Y+sqrt(P)));
    Psi_values = abs(log(q_0./q_1));
    sum11 = 0;
    sum12 = 0;
    C11 = 0;
    C12 = 0;
    Y1 = points(1:end, 2);
    Y2 = points(1:end, 3);
    for i = 1:num_samples
        y1 = Y1(i);
        q_0=1/(sqrt(2))*exp(-sqrt(2)*abs(y1-sqrt(P)));
        q_1=1/(sqrt(2))*exp(-sqrt(2)*abs(y1+sqrt(P)));
        if q_0 < q_1
            t = abs(log(q_0/q_1));
            val = sum(Psi_values <= t) / length(Psi_values);
            sum11 = sum11 + val;
        end
        C11 = C11 + log(1+q_1/q_0);
        y2 = Y2(i);
        q_0=1/(sqrt(2))*exp(-sqrt(2)*abs(y2-sqrt(P)));
        q_1=1/(sqrt(2))*exp(-sqrt(2)*abs(y2+sqrt(P)));
        if q_0 > q_1
            t = abs(log(q_0/q_1));
            val = sum(Psi_values <= t) / length(Psi_values);
            sum12 = sum12 + val;
        end
        C12 = C12 + log(1+q_0/q_1);
        tot_jifen = 1/2 * (sum11/i + sum12/i);
        C1 = 1/2 * (C11/i + C12/i);
        if mod(i, 1e4)==1
            disp('SNR:');
            disp(S);
            disp('Points：');
            disp(i);
            disp('tot_jifen');
            disp(tot_jifen);
            disp('C1');
            disp(C1);
        end
    end
    max_re = -inf;       
    best_theta = NaN;   
    for theta = -2500:0.1:-0.1   
        f = @(t) log(1 + exp(theta*t)); 
        re_int = integral(f, 0, 1); 
        re = (log(2)-re_int+theta*tot_jifen)/log(2);
        if re > max_re
            max_re = re;          
            best_theta = theta;   
        end
    end
    Capability = (log(2)-C1)/log(2);
    data(j,1) = S;
    data(j,2) = max_re;
    data(j,3) = best_theta;
    data(j,4) = Capability;
end
file_path = fullfile(pwd, sprintf('data_Laplacian_BPSK.mat'));
save(file_path, 'data');