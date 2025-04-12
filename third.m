N = 2048; 
M = 600;  
SNR_dB = 10;

% QPSK
QPSK_signal = [1+1j, 1-1j, -1+1j, -1-1j];
S_f = repmat(QPSK_signal, 1, M/4);

%C
X = [0, S_f(301:600), zeros(1,1447), S_f(1:300)];
s_m = ifft(X, N);

% 加入 AWGN 噪聲
r_m = awgn(s_m, SNR_dB, 'measured');

% 前60個s[m] 和 r[m]
figure;
subplot(2,1,1);
plot(1:60, real(s_m(1:60)), 'bo', 'DisplayName', 'Real(s[m])'); hold on;
plot(1:60, real(r_m(1:60)), 'r*', 'DisplayName', 'Real(r[m])');
%plot(1:60, real(n_m(1:60)), 'g*', 'DisplayName', 'Real(n[m])');
xlabel('Sample Index');
ylabel('Amplitude');
title('Real Part of s[m] and r[m]');
legend;
grid on;

subplot(2,1,2);
plot(1:60, imag(s_m(1:60)), 'bo', 'DisplayName', 'Imag(s[m])'); hold on;
plot(1:60, imag(r_m(1:60)), 'r*', 'DisplayName', 'Imag(r[m])');
%plot(1:60, imag(n_m(1:60)), 'g*', 'DisplayName', 'Imag(n[m])');
xlabel('Sample Index');
ylabel('Amplitude');
title('Imaginary Part of s[m] and r[m]');
legend;
grid on;


%s[m] 的平均功率
power_s = mean(abs(s_m).^2);

%10次平均功率
num_trials = 10;
power_s_trials = zeros(1, num_trials);
power_n_trials = zeros(1, num_trials);

for i = 1:num_trials
    r_m = awgn(s_m, SNR_dB, 'measured'); % 重新加入噪聲
    n_m = r_m - s_m; %噪聲信號
    
    power_s_trials(i) = mean(abs(s_m).^2); %傳輸信號功率
    power_n_trials(i) = mean(abs(n_m).^2); %噪聲功率
end

avg_power_s = mean(power_s_trials); %10次平均傳輸信號功率
avg_power_n = mean(power_n_trials); %10次平均噪聲功率

%理論 SNR
SNR_calculated = 10 * log10(avg_power_s / avg_power_n);

disp(['平均傳輸信號功率: ', num2str(avg_power_s)]);
disp(['平均噪聲功率: ', num2str(avg_power_n)]);
disp(['計算出的 SNR (dB): ', num2str(SNR_calculated)]);
disp(['設定的 SNR (dB): ', num2str(SNR_dB)]);