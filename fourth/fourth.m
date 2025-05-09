clc;clear;
N = 2048;
M = 600;     %600 QPSK signal
SNR_dB = 10; %SNR in 3th
num_trials = 10;

% QPSK Mapping
QPSK_signal = [1+1j, 1-1j, -1+1j, -1-1j];  
S_f = repmat(QPSK_signal, 1, M/4);

X = [0, S_f(301:600), zeros(1,1447), S_f(1:300)];
s_m = ifft(X, N);

r_m = awgn(s_m, SNR_dB, 'measured');
R_f = fft(r_m, N);
n_m = r_m - s_m;
N_f = R_f - X;

% ---------- 題目 1 ----------
Re_N = real(N_f);
Im_N = imag(N_f);

mean_Re = mean(Re_N);
var_Re = var(Re_N);
mean_Im = mean(Im_N);
var_Im = var(Im_N);

fprintf('Re(N(f)) -> Mean: %.4f, Variance: %.4f\n', mean_Re, var_Re);
fprintf('Im(N(f)) -> Mean: %.4f, Variance: %.4f\n', mean_Im, var_Im);

% ---------- 題目 2 ----------
figure;
histogram(Re_N, 100, 'Normalization', 'pdf'); hold on;   
x = linspace(min(Re_N), max(Re_N), 1000);
plot(x, normpdf(x, mean_Re, sqrt(var_Re)), 'r', 'LineWidth', 2); 
title('Re(N(f)) vs Gaussian');
legend('Re(N(f)) Histogram', 'Gaussian Fit');
grid on;

% ---------- 題目 3 ----------
figure;
histogram(Im_N, 100, 'Normalization', 'pdf'); hold on;
x = linspace(min(Im_N), max(Im_N), 1000);
plot(x, normpdf(x, mean_Im, sqrt(var_Im)), 'r', 'LineWidth', 2);
title('Im(N(f)) vs Gaussian');
legend('Im(N(f)) Histogram', 'Gaussian Fit');
grid on;

% ---------- 題目 4 SNR ----------
Es_N0_freq = 2/ (var_Re+var_Im);
SNR_est_freq = 10 * log10(Es_N0_freq);
fprintf('實際計算的 SNR (Es/N0): %.4f dB\n', SNR_est_freq);
fprintf('理論設定的 SNR: %.4f dB\n', SNR_dB);
fprintf('差異: %.4f dB\n', abs(SNR_est_freq - SNR_dB));


% ---------- 題目 5  ----------
N = 2048;               
M = 600;                
SNR_dB_set = 2:2:14;    
num_bits = M * 2;       
num_trials = 1e3;        %試驗次數 (如果想跑快一點改1000 大概四個點)
QPSK_signal = [1+1j, 1-1j, -1+1j, -1-1j]; % QPSK星座圖

% ===== 預生成正確的1200資料 =====
rng(0);
tx_bits = randi([0 1], 1, num_bits);       
tx_symbols = reshape(tx_bits, 2, [])';     
symbol_idx = bi2de(tx_symbols, 'left-msb')';
S_f = QPSK_signal(symbol_idx + 1);         % 固定頻域信號


BER_sim = zeros(1, length(SNR_dB_set));
EbN0_measured = zeros(1, length(SNR_dB_set));


for k = 1:length(SNR_dB_set)
    SNR_dB = SNR_dB_set(k);
    total_var_Re = 0;   
    total_var_Im = 0;   
    error_count = 0;   
    
    for trial = 1:num_trials
        
        X = [0, S_f(301:600), zeros(1, 1447), S_f(1:300)];
        s_m = ifft(X, N); 
        
        % ===== 加噪訊 (改跟前面算法一樣) =====
        r_m = awgn(s_m, SNR_dB, 'measured'); 
        R_f = fft(r_m, N);                   
        
        
        N_f = R_f - X; 
        valid_indices = [2:301, 1749:2048];   % 有效子載波
        N_f_valid = N_f(valid_indices);       % 有效雜訊
        
        total_var_Re = total_var_Re + var(real(N_f_valid));
        total_var_Im = total_var_Im + var(imag(N_f_valid));
        
       
        received_symbols = R_f(valid_indices); 
        
        % 接收符號 = [S_f(301:600), S_f(1:300)]  
        % 原始符號 = [S_f(1:300), S_f(301:600)]

        received_symbols = [received_symbols(301:600), received_symbols(1:300)]; % 前後交換
        
        % 最小歐式距離來判斷
        rx_symbol_idx = zeros(1, M);
        for i = 1:M
            [~, idx] = min(abs(received_symbols(i) - QPSK_signal));
            rx_symbol_idx(i) = idx - 1;
        end
        
        
        rx_bits = reshape(de2bi(rx_symbol_idx, 2, 'left-msb')', 1, []);
        error_count = error_count + sum(rx_bits ~= tx_bits);
    end
    
    % ===== 計算Eb/N0與BER =====
    var_Re = total_var_Re / num_trials;
    var_Im = total_var_Im / num_trials;
    %N0 = var_Re + var_Im;           
    %Es = mean(abs(S_f).^2);          
    EbN0 = (1/var_Re) / log2(4);        
    EbN0_measured(k) = 10 * log10(EbN0);
    BER_sim(k) = error_count / (num_bits * num_trials);
    % ===== BER曲線 =====
EbN0_lin = 10.^(EbN0_measured/10);
BER_theory = 0.5 * erfc(sqrt(EbN0_lin)); % QPSK理論BER
    fprintf('SNR = %2d dB  var_Re =%.4f  Eb/N0 = %.2f dB, BER = %.4e  理論BER = %.4e\n', SNR_dB,var_Re, EbN0_measured(k), BER_sim(k),BER_theory(k));

end



figure;
semilogy(EbN0_measured, BER_sim, 'bo-', 'LineWidth', 2); hold on;
semilogy(EbN0_measured, BER_theory, 'r--', 'LineWidth', 2);
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('模擬結果', '理論曲線', 'Location', 'southwest');
title('QPSK BER性能 (修正符號順序後)');
grid on;