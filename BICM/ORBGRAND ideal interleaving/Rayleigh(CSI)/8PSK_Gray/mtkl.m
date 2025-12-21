clear;
clc;
num_samples1 = 1e3;
num_samples2 = 1e4;
SNR=-10:2:30;
data = zeros(length(SNR),12);
for j = 1:length(SNR)
    S = SNR(j);
    a = sqrt(10^(S/10));
    sigma = sqrt(0.5);
    clear points;
    points = complex(double(zeros(num_samples1*num_samples2,9)));
    H = raylrnd(0.71, 1, num_samples1);
    Phi = -pi + 2*pi*rand(1, num_samples1);
    index = 1;
    for k = 1:num_samples1
        h = H(k);
        phi = Phi(k);
        angle_set = [0, 1, 2, 3, 4, 5, 6, -1]*pi/4 + phi;
        group = zeros(num_samples2, 9);  % [h, phi, Z, Z1, Z2, Z3, Z4]
        y = generate_symbols(angle_set, [0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125], h, a, sigma, num_samples2);
        % q+(y)
        y1 = generate_symbols(angle_set([3,4,5,6]), [0.25,0.25,0.25,0.25], h, a, sigma, num_samples2);
        % q-(y)
        y2 = generate_symbols(angle_set([2,1,8,7]), [0.25,0.25,0.25,0.25], h, a, sigma, num_samples2);
        % q+(y)
        y3 = generate_symbols(angle_set([5,6,7,8]), [0.25,0.25,0.25,0.25], h, a, sigma, num_samples2);
        % q-(y)
        y4 = generate_symbols(angle_set([1,2,3,4]), [0.25,0.25,0.25,0.25], h, a, sigma, num_samples2);
        % q+(y)
        y5 = generate_symbols(angle_set([5,4,1,8]), [0.25,0.25,0.25,0.25], h, a, sigma, num_samples2);
        % q-(y)
        y6 = generate_symbols(angle_set([3,2,7,6]), [0.25,0.25,0.25,0.25], h, a, sigma, num_samples2);
        group(:,1) = h;
        group(:,2) = phi;
        group(:,3) = y;
        group(:,4) = y1;
        group(:,5) = y2;
        group(:,6) = y3;
        group(:,7) = y4;
        group(:,8) = y5;
        group(:,9) = y6;
        points(index:index+num_samples2-1,:) = group;
        index = index + num_samples2;
    end
    H = points(1:end, 1);
    Phi = points(1:end, 2);
    Y = points(1:end, 3);
    num_samples = length(Y);
    sum11 = 0;
    sum12 = 0;
    sum21 = 0;
    sum22 = 0;
    sum31 = 0;
    sum32 = 0;
    C11 = 0;
    C12 = 0;
    C21 = 0;
    C22 = 0;
    C31 = 0;
    C32 = 0;
    Y1 = points(1:end, 4);
    Y2 = points(1:end, 5);
    Y3 = points(1:end, 6);
    Y4 = points(1:end, 7);
    Y5 = points(1:end, 8);
    Y6 = points(1:end, 9);
    
    det = 20;
    idx  = (0:num_samples2:(num_samples1-1)*num_samples2)' + (1:det);
    Y_idx = Y(idx(:));
    H_idx = H(idx(:));
    Phi_idx = Phi(idx(:));
    Psi1_values = abs(log((exp(-(abs(Y_idx-(H_idx*a.*(cos(2*pi/4+Phi_idx)+1i*sin(2*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(3*pi/4+Phi_idx)+1i*sin(3*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(4*pi/4+Phi_idx)+1i*sin(4*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(5*pi/4+Phi_idx)+1i*sin(5*pi/4+Phi_idx))))).^2))./(exp(-(abs(Y_idx-(H_idx*a.*(cos(1*pi/4+Phi_idx)+1i*sin(1*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(0*pi/4+Phi_idx)+1i*sin(0*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(-1*pi/4+Phi_idx)+1i*sin(-1*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(6*pi/4+Phi_idx)+1i*sin(6*pi/4+Phi_idx))))).^2))));
    Psi2_values = abs(log((exp(-(abs(Y_idx-(H_idx*a.*(cos(4*pi/4+Phi_idx)+1i*sin(4*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(5*pi/4+Phi_idx)+1i*sin(5*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(6*pi/4+Phi_idx)+1i*sin(6*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(-1*pi/4+Phi_idx)+1i*sin(-1*pi/4+Phi_idx))))).^2))./(exp(-(abs(Y_idx-(H_idx*a.*(cos(3*pi/4+Phi_idx)+1i*sin(3*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(2*pi/4+Phi_idx)+1i*sin(2*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(1*pi/4+Phi_idx)+1i*sin(1*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(0*pi/4+Phi_idx)+1i*sin(0*pi/4+Phi_idx))))).^2))));
    Psi3_values = abs(log((exp(-(abs(Y_idx-(H_idx*a.*(cos(4*pi/4+Phi_idx)+1i*sin(4*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(3*pi/4+Phi_idx)+1i*sin(3*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(0*pi/4+Phi_idx)+1i*sin(0*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(-1*pi/4+Phi_idx)+1i*sin(-1*pi/4+Phi_idx))))).^2))./(exp(-(abs(Y_idx-(H_idx*a.*(cos(2*pi/4+Phi_idx)+1i*sin(2*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(1*pi/4+Phi_idx)+1i*sin(1*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(6*pi/4+Phi_idx)+1i*sin(6*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(H_idx*a.*(cos(5*pi/4+Phi_idx)+1i*sin(5*pi/4+Phi_idx))))).^2))));
    
    
    for i =1:num_samples
        h = H(i);
        phi = Phi(i);
        y1 = Y1(i);
        q_0=1/(4*pi)*(exp(-(abs(y1-(h*a.*(cos(2*pi/4+phi)+1i*sin(2*pi/4+phi))))).^2)+exp(-(abs(y1-(h*a.*(cos(3*pi/4+phi)+1i*sin(3*pi/4+phi))))).^2)+exp(-(abs(y1-(h*a.*(cos(4*pi/4+phi)+1i*sin(4*pi/4+phi))))).^2)+exp(-(abs(y1-(h*a.*(cos(5*pi/4+phi)+1i*sin(5*pi/4+phi))))).^2));
        q_1=1/(4*pi)*(exp(-(abs(y1-(h*a.*(cos(1*pi/4+phi)+1i*sin(1*pi/4+phi))))).^2)+exp(-(abs(y1-(h*a.*(cos(0*pi/4+phi)+1i*sin(0*pi/4+phi))))).^2)+exp(-(abs(y1-(h*a.*(cos(-1*pi/4+phi)+1i*sin(-1*pi/4+phi))))).^2)+exp(-(abs(y1-(h*a.*(cos(6*pi/4+phi)+1i*sin(6*pi/4+phi))))).^2));
        if q_0 < q_1
            t = abs(log(q_0/q_1));
            val = sum(Psi1_values <= t) / length(Psi1_values);
            sum11 = sum11 + val;
        end
        C11 = C11 + log(1+q_1/q_0);
        y2 = Y2(i);
        q_0=1/(4*pi)*(exp(-(abs(y2-(h*a.*(cos(2*pi/4+phi)+1i*sin(2*pi/4+phi))))).^2)+exp(-(abs(y2-(h*a.*(cos(3*pi/4+phi)+1i*sin(3*pi/4+phi))))).^2)+exp(-(abs(y2-(h*a.*(cos(4*pi/4+phi)+1i*sin(4*pi/4+phi))))).^2)+exp(-(abs(y2-(h*a.*(cos(5*pi/4+phi)+1i*sin(5*pi/4+phi))))).^2));
        q_1=1/(4*pi)*(exp(-(abs(y2-(h*a.*(cos(1*pi/4+phi)+1i*sin(1*pi/4+phi))))).^2)+exp(-(abs(y2-(h*a.*(cos(0*pi/4+phi)+1i*sin(0*pi/4+phi))))).^2)+exp(-(abs(y2-(h*a.*(cos(-1*pi/4+phi)+1i*sin(-1*pi/4+phi))))).^2)+exp(-(abs(y2-(h*a.*(cos(6*pi/4+phi)+1i*sin(6*pi/4+phi))))).^2));
        if q_0 > q_1
            t = abs(log(q_0/q_1));
            val = sum(Psi1_values <= t) / length(Psi1_values);
            sum12 = sum12 + val;
        end
        C12 = C12 + log(1+q_0/q_1);
        y3 = Y3(i);
        q_0=1/(4*pi)*(exp(-(abs(y3-(h*a.*(cos(4*pi/4+phi)+1i*sin(4*pi/4+phi))))).^2)+exp(-(abs(y3-(h*a.*(cos(5*pi/4+phi)+1i*sin(5*pi/4+phi))))).^2)+exp(-(abs(y3-(h*a.*(cos(6*pi/4+phi)+1i*sin(6*pi/4+phi))))).^2)+exp(-(abs(y3-(h*a.*(cos(-1*pi/4+phi)+1i*sin(-1*pi/4+phi))))).^2));
        q_1=1/(4*pi)*(exp(-(abs(y3-(h*a.*(cos(3*pi/4+phi)+1i*sin(3*pi/4+phi))))).^2)+exp(-(abs(y3-(h*a.*(cos(2*pi/4+phi)+1i*sin(2*pi/4+phi))))).^2)+exp(-(abs(y3-(h*a.*(cos(1*pi/4+phi)+1i*sin(1*pi/4+phi))))).^2)+exp(-(abs(y3-(h*a.*(cos(0*pi/4+phi)+1i*sin(0*pi/4+phi))))).^2));
        if q_0 < q_1
            t = abs(log(q_0/q_1));
            val = sum(Psi2_values <= t) / length(Psi2_values);
            sum21 = sum21 + val;
        end
        C21 = C21 + log(1+q_1/q_0);
        y4 = Y4(i);
        q_0=1/(4*pi)*(exp(-(abs(y4-(h*a.*(cos(4*pi/4+phi)+1i*sin(4*pi/4+phi))))).^2)+exp(-(abs(y4-(h*a.*(cos(5*pi/4+phi)+1i*sin(5*pi/4+phi))))).^2)+exp(-(abs(y4-(h*a.*(cos(6*pi/4+phi)+1i*sin(6*pi/4+phi))))).^2)+exp(-(abs(y4-(h*a.*(cos(-1*pi/4+phi)+1i*sin(-1*pi/4+phi))))).^2));
        q_1=1/(4*pi)*(exp(-(abs(y4-(h*a.*(cos(3*pi/4+phi)+1i*sin(3*pi/4+phi))))).^2)+exp(-(abs(y4-(h*a.*(cos(2*pi/4+phi)+1i*sin(2*pi/4+phi))))).^2)+exp(-(abs(y4-(h*a.*(cos(1*pi/4+phi)+1i*sin(1*pi/4+phi))))).^2)+exp(-(abs(y4-(h*a.*(cos(0*pi/4+phi)+1i*sin(0*pi/4+phi))))).^2));
        if q_0 > q_1
            t = abs(log(q_0/q_1));
            val = sum(Psi2_values <= t) / length(Psi2_values);
            sum22 = sum22 + val;
        end
        C22 = C22 + log(1+q_0/q_1);
        y5 = Y5(i);
        q_0=1/(4*pi)*(exp(-(abs(y5-(h*a.*(cos(4*pi/4+phi)+1i*sin(4*pi/4+phi))))).^2)+exp(-(abs(y5-(h*a.*(cos(3*pi/4+phi)+1i*sin(3*pi/4+phi))))).^2)+exp(-(abs(y5-(h*a.*(cos(0*pi/4+phi)+1i*sin(0*pi/4+phi))))).^2)+exp(-(abs(y5-(h*a.*(cos(-1*pi/4+phi)+1i*sin(-1*pi/4+phi))))).^2));
        q_1=1/(4*pi)*(exp(-(abs(y5-(h*a.*(cos(2*pi/4+phi)+1i*sin(2*pi/4+phi))))).^2)+exp(-(abs(y5-(h*a.*(cos(1*pi/4+phi)+1i*sin(1*pi/4+phi))))).^2)+exp(-(abs(y5-(h*a.*(cos(6*pi/4+phi)+1i*sin(6*pi/4+phi))))).^2)+exp(-(abs(y5-(h*a.*(cos(5*pi/4+phi)+1i*sin(5*pi/4+phi))))).^2));
        if q_0 < q_1
            t = abs(log(q_0/q_1));
            val = sum(Psi3_values <= t) / length(Psi3_values);
            sum31 = sum31 + val;
        end
        C31 = C31 + log(1+q_1/q_0);
        y6 = Y6(i);
        q_0=1/(4*pi)*(exp(-(abs(y6-(h*a.*(cos(4*pi/4+phi)+1i*sin(4*pi/4+phi))))).^2)+exp(-(abs(y6-(h*a.*(cos(3*pi/4+phi)+1i*sin(3*pi/4+phi))))).^2)+exp(-(abs(y6-(h*a.*(cos(0*pi/4+phi)+1i*sin(0*pi/4+phi))))).^2)+exp(-(abs(y6-(h*a.*(cos(-1*pi/4+phi)+1i*sin(-1*pi/4+phi))))).^2));
        q_1=1/(4*pi)*(exp(-(abs(y6-(h*a.*(cos(2*pi/4+phi)+1i*sin(2*pi/4+phi))))).^2)+exp(-(abs(y6-(h*a.*(cos(1*pi/4+phi)+1i*sin(1*pi/4+phi))))).^2)+exp(-(abs(y6-(h*a.*(cos(6*pi/4+phi)+1i*sin(6*pi/4+phi))))).^2)+exp(-(abs(y6-(h*a.*(cos(5*pi/4+phi)+1i*sin(5*pi/4+phi))))).^2));
        if q_0 > q_1
            t = abs(log(q_0/q_1));
            val = sum(Psi3_values <= t) / length(Psi3_values);
            sum32 = sum32 + val;
        end
        C32 = C32 + log(1+q_0/q_1);
        tot_jifen1 = 1/2 * (sum11/i + sum12/i);
        C1 = 1/2 * (C11/i + C12/i);
        tot_jifen2 = 1/2 * (sum21/i + sum22/i);
        C2 = 1/2 * (C21/i + C22/i);
        tot_jifen3 = 1/2 * (sum31/i + sum32/i);
        C3 = 1/2 * (C31/i + C32/i);
        if mod(i, 1e4)==1
            disp('SNR:');
            disp(S);
            disp('Points：');
            disp(i);
            disp('tot_jifen1');
            disp(tot_jifen1);
            disp('C1');
            disp(C1);
            disp('tot_jifen2');
            disp(tot_jifen2);
            disp('C2');
            disp(C2);
            disp('tot_jifen3');
            disp(tot_jifen3);
            disp('C3');
            disp(C3);
        end
    end
    max_re1 = -inf;
    best_theta1 = NaN;
    for theta = -2500:0.1:-0.1
        f = @(t) log(1 + exp(theta*t));
        re_int = integral(f, 0, 1);
        re1 = (log(2)-re_int+theta*tot_jifen1)/log(2);
        if re1 > max_re1
            max_re1 = re1;
            best_theta1 = theta;
        end
    end
    Capability1 = (log(2)-C1)/log(2);
    max_re2 = -inf;
    best_theta2 = NaN;
    for theta = -2500:0.1:-0.1
        f = @(t) log(1 + exp(theta*t));
        re_int = integral(f, 0, 1);
        re2 = (log(2)-re_int+theta*tot_jifen2)/log(2);
        if re2 > max_re2
            max_re2 = re2;
            best_theta2 = theta;
        end
    end
    Capability2 = (log(2)-C2)/log(2);
    max_re3 = -inf;
    best_theta3 = NaN;
    for theta = -2500:0.1:-0.1
        f = @(t) log(1 + exp(theta*t));
        re_int = integral(f, 0, 1);
        re3 = (log(2)-re_int+theta*tot_jifen3)/log(2);
        if re3 > max_re3
            max_re3 = re3;
            best_theta3 = theta;
        end
    end
    Capability3 = (log(2)-C3)/log(2);
    data(j,1) = S;
    data(j,2) = max_re1;
    data(j,3) = best_theta1;
    data(j,4) = Capability1;
    data(j,5) = max_re2;
    data(j,6) = best_theta2;
    data(j,7) = Capability2;
    data(j,8) = max_re3;
    data(j,9) = best_theta3;
    data(j,10) = Capability3;
    data(j,11) = max_re1+max_re2+max_re3;
    data(j,12) = Capability1+Capability2+Capability3;
end
file_path = fullfile(pwd, sprintf('data_Rayleigh_8PSK_Gray.mat'));
save(file_path, 'data');