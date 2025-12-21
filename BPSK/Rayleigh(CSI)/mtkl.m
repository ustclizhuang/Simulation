clear;
clc;
SNR = 0:1:10;
data = zeros(length(SNR),4);
num_samples1 = 5e3;
num_samples2 = 5e3;
for j = 1:length(SNR)
    clear points;
    points = zeros(num_samples1*num_samples2,4);
    S = SNR(j);
    P = 10^(S/10);
    H = raylrnd(0.71, 1, num_samples1);
    index = 1;
    for i = 1:num_samples1
        h = H(i);
        group = zeros(num_samples2,4);
        cdf = cumsum([0.5,0.5]);
        U = rand(1, num_samples2);
        idx = arrayfun(@(u) find(u < cdf, 1), U);
        angle = [+1, -1];
        group(:,1) = h;
        group(:,2) = angle(idx).*h.*sqrt(P)+randn(1, num_samples2);
        group(:,3) = h.*sqrt(P)+randn(1, num_samples2);
        group(:,4) = -h.*sqrt(P)+randn(1, num_samples2);
        points(index:index+num_samples2-1,:) = group;
        index = index + num_samples2;
    end
    % file_path = fullfile(pwd, sprintf('points_BPSK_SNR%d.mat',S));
    % save(file_path, 'points');
    H = points(1:end, 1);
    Y = points(1:end, 2);
    det = 20;
    idx  = (0:num_samples2:(num_samples1-1)*num_samples2)' + (1:det);
    Y_idx = Y(idx(:));
    H_idx = H(idx(:));
    Psi_values = abs(2.*Y_idx.*H_idx.*sqrt(P));
    num_samples = length(Y);
    sum11 = 0;
    sum12 = 0;
    C11 = 0;
    C12 = 0;
    Y1 = points(1:end, 3);
    Y2 = points(1:end, 4);
    for i = 1:num_samples
        h = H(i);
        y1 = Y1(i);
        q_0=1/(sqrt(2*pi))*exp(-(y1-h*sqrt(P))^2/2);
        q_1=1/(sqrt(2*pi))*exp(-(y1+h*sqrt(P))^2/2);
        if q_0 < q_1
            t = abs(log(q_0/q_1));
            val = sum(Psi_values <= t) / length(Psi_values);
            sum11 = sum11 + val;
        end
        C11 = C11 + log(1+q_1/q_0);
        y2 = Y2(i);
        q_0=1/(sqrt(2*pi))*exp(-(y2-h*sqrt(P))^2/2);
        q_1=1/(sqrt(2*pi))*exp(-(y2+h*sqrt(P))^2/2);
        if q_0 > q_1
            t = abs(log(q_0/q_1));
            val = sum(Psi_values <= t) / length(Psi_values);
            sum12 = sum12 + val;
        end
        C12 = C12 + log(1+q_0/q_1);
        tot_jifen = 1/2 * (sum11/i + sum12/i);
        C1 = 1/2 * (C11/i + C12/i);
        if mod(i, 1e4)==1
            disp('SNR：');
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
file_path = fullfile(pwd, sprintf('data_Rayleigh_BPSK.mat'));
save(file_path, 'data');