clear;
clc;
num_samples1 = 1e3;
num_samples2 = 1e4;
SNR=-10:2:30;
data = zeros(length(SNR),9);
for j = 1:length(SNR)
    S = SNR(j);
    a = sqrt(10^(S/10)/2);
    sigma = sqrt(0.5);
    clear points;
    points = complex(double(zeros(num_samples1*num_samples2,7)));
    H = raylrnd(0.71, 1, num_samples1);
    Phi = -pi + 2*pi*rand(1, num_samples1);
    index = 1;
    for k = 1:num_samples1
        h = H(k);
        phi = Phi(k);
        angle_set = [5, 3, 1, -1]*pi/4 + phi;
        group = zeros(num_samples2, 7);  % [h, phi, Z, Z1, Z2, Z3, Z4]
        % q(y|h)
        y = generate_symbols(angle_set, [0.25, 0.25, 0.25, 0.25], h, a, sigma, num_samples2);
        % q+(y)
        y1 = generate_symbols(angle_set([1,4]), [0.5,0.5], h, a, sigma, num_samples2);
        % q-(y)
        y2 = generate_symbols(angle_set([3,2]), [0.5,0.5], h, a, sigma, num_samples2);
        % q+(y)
        y3 = generate_symbols(angle_set([1,2]), [0.5,0.5], h, a, sigma, num_samples2);
        % q-(y)
        y4 = generate_symbols(angle_set([3,4]), [0.5,0.5], h, a, sigma, num_samples2);
        group(:,1) = h;
        group(:,2) = phi;
        group(:,3) = y;
        group(:,4) = y1;
        group(:,5) = y2;
        group(:,6) = y3;
        group(:,7) = y4;
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
    C11 = 0;
    C12 = 0;
    C21 = 0;
    C22 = 0;
    Y1 = points(1:end, 4);
    Y2 = points(1:end, 5);
    Y3 = points(1:end, 6);
    Y4 = points(1:end, 7);
    det = 20;
    idx  = (0:num_samples2:(num_samples1-1)*num_samples2)' + (1:det);
    Y_idx = Y(idx(:));
    H_idx = H(idx(:));
    Phi_idx = Phi(idx(:));

    Psi1_values = abs(log((exp(-(abs(Y_idx-(sqrt(2)*H_idx*a.*(cos(5*pi/4+Phi_idx)+1i*sin(5*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(sqrt(2)*H_idx*a.*(cos(-1*pi/4+Phi_idx)+1i*sin(-1*pi/4+Phi_idx))))).^2))./(exp(-(abs(Y_idx-(sqrt(2)*H_idx*a.*(cos(1*pi/4+Phi_idx)+1i*sin(1*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(sqrt(2)*H_idx*a.*(cos(3*pi/4+Phi_idx)+1i*sin(3*pi/4+Phi_idx))))).^2))));
    Psi2_values = abs(log((exp(-(abs(Y_idx-(sqrt(2)*H_idx*a.*(cos(3*pi/4+Phi_idx)+1i*sin(3*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(sqrt(2)*H_idx*a.*(cos(5*pi/4+Phi_idx)+1i*sin(5*pi/4+Phi_idx))))).^2))./(exp(-(abs(Y_idx-(sqrt(2)*H_idx*a.*(cos(1*pi/4+Phi_idx)+1i*sin(1*pi/4+Phi_idx))))).^2)+exp(-(abs(Y_idx-(sqrt(2)*H_idx*a.*(cos(-1*pi/4+Phi_idx)+1i*sin(-1*pi/4+Phi_idx))))).^2))));

    
    
    for i =1:num_samples
        h = H(i);
        phi = Phi(i);
        y1 = Y1(i);
        q_0=1/(2*pi)*(exp(-(abs(y1-(sqrt(2)*h*a.*(cos(5*pi/4+phi)+1i*sin(5*pi/4+phi))))).^2)+exp(-(abs(y1-(sqrt(2)*h*a.*(cos(-1*pi/4+phi)+1i*sin(-1*pi/4+phi))))).^2));
        q_1=1/(2*pi)*(exp(-(abs(y1-(sqrt(2)*h*a.*(cos(1*pi/4+phi)+1i*sin(1*pi/4+phi))))).^2)+exp(-(abs(y1-(sqrt(2)*h*a.*(cos(3*pi/4+phi)+1i*sin(3*pi/4+phi))))).^2));
        if q_0 < q_1
            t = abs(log(q_0/q_1));
            val = sum(Psi1_values <= t) / length(Psi1_values);
            sum11 = sum11 + val;
        end
        C11 = C11 + log(1+q_1/q_0);
        y2 = Y2(i);
        q_0=1/(2*pi)*(exp(-(abs(y2-(sqrt(2)*h*a.*(cos(5*pi/4+phi)+1i*sin(5*pi/4+phi))))).^2)+exp(-(abs(y2-(sqrt(2)*h*a.*(cos(-1*pi/4+phi)+1i*sin(-1*pi/4+phi))))).^2));
        q_1=1/(2*pi)*(exp(-(abs(y2-(sqrt(2)*h*a.*(cos(1*pi/4+phi)+1i*sin(1*pi/4+phi))))).^2)+exp(-(abs(y2-(sqrt(2)*h*a.*(cos(3*pi/4+phi)+1i*sin(3*pi/4+phi))))).^2));
        if q_0 > q_1
            t = abs(log(q_0/q_1));
            val = sum(Psi1_values <= t) / length(Psi1_values);
            sum12 = sum12 + val;
        end
        C12 = C12 + log(1+q_0/q_1);
        y3 = Y3(i);
        q_0=1/(2*pi)*(exp(-(abs(y3-(sqrt(2)*h*a.*(cos(3*pi/4+phi)+1i*sin(3*pi/4+phi))))).^2)+exp(-(abs(y3-(sqrt(2)*h*a.*(cos(5*pi/4+phi)+1i*sin(5*pi/4+phi))))).^2));
        q_1=1/(2*pi)*(exp(-(abs(y3-(sqrt(2)*h*a.*(cos(1*pi/4+phi)+1i*sin(1*pi/4+phi))))).^2)+exp(-(abs(y3-(sqrt(2)*h*a.*(cos(-1*pi/4+phi)+1i*sin(-1*pi/4+phi))))).^2));
        if q_0 < q_1
            t = abs(log(q_0/q_1));
            val = sum(Psi2_values <= t) / length(Psi2_values);
            sum21 = sum21 + val;
        end
        C21 = C21 + log(1+q_1/q_0);
        y4 = Y4(i);
        q_0=1/(2*pi)*(exp(-(abs(y4-(sqrt(2)*h*a.*(cos(3*pi/4+phi)+1i*sin(3*pi/4+phi))))).^2)+exp(-(abs(y4-(sqrt(2)*h*a.*(cos(5*pi/4+phi)+1i*sin(5*pi/4+phi))))).^2));
        q_1=1/(2*pi)*(exp(-(abs(y4-(sqrt(2)*h*a.*(cos(1*pi/4+phi)+1i*sin(1*pi/4+phi))))).^2)+exp(-(abs(y4-(sqrt(2)*h*a.*(cos(-1*pi/4+phi)+1i*sin(-1*pi/4+phi))))).^2));
        if q_0 > q_1
            t = abs(log(q_0/q_1));
            val = sum(Psi2_values <= t) / length(Psi2_values);
            sum22 = sum22 + val;
        end
        C22 = C22 + log(1+q_0/q_1);
        tot_jifen1 = 1/2 * (sum11/i + sum12/i);
        C1 = 1/2 * (C11/i + C12/i);
        tot_jifen2 = 1/2 * (sum21/i + sum22/i);
        C2 = 1/2 * (C21/i + C22/i);
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
    data(j,1) = S;
    data(j,2) = max_re1;
    data(j,3) = best_theta1;
    data(j,4) = Capability1;
    data(j,5) = max_re2;
    data(j,6) = best_theta2;
    data(j,7) = Capability2;
    data(j,8) = max_re1+max_re2;
    data(j,9) = Capability1+Capability2;
end
file_path = fullfile(pwd, sprintf('data_Rayleigh_4PSK_Gray.mat'));
save(file_path, 'data');