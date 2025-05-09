N = 2048;           
M = 600;            % 有效子載波
cp_len = 160;       % CP長度
SNR_dB = 15;        

% ------ Doopler參數------
fc = 1.8e9;         % 載波頻率（1.8 GHz）
v_km_hour = 3;      % 移動速度（3 km/h）
c = 3e8;            % 光速
v_m_sec = v_km_hour*1000/3600;    % 速度轉換為公尺/秒
fd = v_m_sec*fc/c;  % 計算Doopler Freq

% ================== 生成QPSK訊號 ==================
QPSK_table = (1/sqrt(2)) * [1+1j, 1-1j, -1+1j, -1-1j]; % 單位功率歸一化
S_f = QPSK_table(randi(4, 1, M));  % 隨機生成600個QPSK符號

% ================== 頻域配置 ==================
X = zeros(N, 1);                    % 頻域訊號（列向量）
start_idx = N/2 - M/2 + 1;          % 有效子載波起始
X(start_idx:start_idx + M - 1) = S_f.'; % 將QPSK符號居中放置

% ================== IFFT生成時域訊號 ==================
s_m = ifft(X, N);       % 生成時域訊號s[m]

% ================== 添加循環前綴 ==================
s_cp = [s_m(1889:2048,1); s_m]; % 添加CP

% ================== 多徑通道配置==================
seednum = randi([1 2^31-7],1,1);     % 隨機種子
chcfg.Seed = seednum;                % 設定通道種子
chcfg.NRxAnts = 1;                   % 接收端單天線數量
chcfg.NormalizeTxAnts = 'Off';       % 關閉天線歸一化
chcfg.DelayProfile = 'EVA';          % 使用EVA多徑模型
chcfg.MIMOCorrelation = 'High';      % 高天線相關性
chcfg.DopplerFreq = fd;              % 帶入fd
chcfg.SamplingRate = 30.72e6;        % 取樣率30.72 MHz
chcfg.InitTime = 0;                  % 初始時間，每次只傳一個frame 所以=0
chcfg.NTerms = 16;                   % 預設值
chcfg.ModelType = 'GMEDS';           % 通道模型類型
chcfg.NormalizePathGains = 'On';     % 路徑增益歸一化
chcfg.InitPhase = 'Random';          % 初始相位隨機

% ================== 通過多徑通道 ==================
% 使用LTE的lteFadingChannel函數
[u, ~] = lteFadingChannel(chcfg, s_cp); % 輸出為列向量

% ==================================================================
%                       問題1：無雜訊通道估計
% ==================================================================
% --- 去除CP並執行FFT ---
r1 = u(cp_len + 1:end);           % 去除CP
R_f1 = fft(r1, N);                % 頻域訊號

% --- 計算通道響應H(f) ---
valid_idx = start_idx:start_idx + M - 1;  % 有效子載波索引
H_f = zeros(N, 1);
H_f(valid_idx) = R_f1(valid_idx) ./ X(valid_idx);

% --- 繪製H(f)的幅度與相位 ---
figure('Position', [100 100 800 600]);
subplot(2,1,1);
plot(abs(H_f(valid_idx)), 'LineWidth', 1.5);
title('|H(f)| (無雜訊)'); xlabel('子載波索引'); grid on;
subplot(2,1,2);
plot(angle(H_f(valid_idx)), 'LineWidth', 1.5);
title('∠H(f) (無雜訊)'); xlabel('子載波索引'); grid on;

% ==================================================================
%          問題2：有雜訊下的星座圖比較（有fading vs 無fading）
% ==================================================================
% --- 計算訊號功率 ---
s_avg_power = mean(abs(s_cp).^2);
s_avg_power_dB = 10*log10(s_avg_power);

% --- 添加AWGN雜訊（有fading）---
r_noise_fading = awgn(u, SNR_dB, s_avg_power_dB);
r2 = r_noise_fading(cp_len + 1:end);
R_f2 = fft(r2, N);

% --- 繪製有fading的星座圖 ---
figure;
scatter(real(R_f2(valid_idx)), imag(R_f2(valid_idx)), 10, 'filled');
axis equal; grid on;
title(['R(f) 星座圖 (有衰落, SNR=', num2str(SNR_dB), 'dB)']);
xlabel('實部'); ylabel('虛部');

% --- 無fading對照組：直接添加雜訊 ---
r_no_fading = awgn(s_cp, SNR_dB, s_avg_power_dB);
R_f_no_fading = fft(r_no_fading(cp_len + 1:end), N);

% --- 繪製無fading的星座圖 ---
figure;
scatter(real(R_f_no_fading(valid_idx)), imag(R_f_no_fading(valid_idx)), 10, 'filled');
axis equal; grid on;
title(['R(f) 星座圖 (無衰落, SNR=', num2str(SNR_dB), 'dB)']);
xlabel('實部'); ylabel('虛部');

% ==================================================================
%                 問題3：通道補償後的星座圖
% ==================================================================
% --- R'(f) = R(f)/H(f)---
R_compensated = zeros(N, 1);
R_compensated(valid_idx) = R_f2(valid_idx) ./ H_f(valid_idx);

% --- 繪製補償後的星座圖 ---
figure;
scatter(real(R_compensated(valid_idx)), imag(R_compensated(valid_idx)), 10, 'filled');
hold on;
plot(QPSK_table, 'ro', 'MarkerSize', 10, 'LineWidth', 2); % 理想QPSK位置
axis equal; grid on;
title('補償後 R''(f) 星座圖');
xlabel('實部'); ylabel('虛部');
legend('接收符號', '理想位置');