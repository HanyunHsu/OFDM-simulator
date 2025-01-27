fs = 1000;
t = 0:1/fs:1; 
f = [16, 18, 22, 24]; 
a = [1, -1, 1, -1]; 
b = [1, 1, -1, 1];
m_t = zeros(size(t)); 
fc = 20;
for n = 1:4 
    m_t = m_t + a(n) * cos(2*pi*f(n)*t) - b(n) * sin(2*pi*f(n)*t);
end

%第一題的圖 m(t)
figure; 
plot(t, m_t, 'LineWidth', 1.5); 
title('Band-pass Signal m(t)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on; %顯示網格線

m_e_t = zeros(size(t));%基頻等效


for n = 1:length(f)
    m_e_t = m_e_t + (a(n)+1j*b(n))*exp(1j*2*pi*(f(n)-20)*t);
end
m_I = real(m_e_t);
m_Q = imag(m_e_t);
r = sqrt(m_I.^2 + m_Q.^2);
%直接取RE 跟IM畫圖
figure;
plot(t, m_I, 'r', 'LineWidth', 1.5); hold on;
plot(t, m_Q, 'b', 'LineWidth', 1.5);
plot(t, r, 'g', 'LineWidth', 1.5);
plot(t, -r, 'g--', 'LineWidth', 1.5); % 負包絡
legend('m_I(t)', 'm_Q(t)', 'r(t)');
xlabel('Time (s)');
ylabel('Amplitude');
title('In-Phase Component, Quadrature Component, and Envelope');
grid on;

m_It = zeros(size(t)); 
m_Qt = zeros(size(t)); 
for n = 1:length(f)
    m_It= m_It + (a(n)*cos(2*pi*f(n)*t) - b(n)*sin(2*pi*f(n)*t));
    m_Qt= m_Qt+(a(n)*sin(2*pi*f(n)*t) + b(n)*cos(2*pi*f(n)*t));
end

%第七題的圖 r(t) m_I(t) m_Q(t)
figure;
plot(t, m_It, 'r', 'LineWidth', 1.5); hold on;
plot(t, m_Qt, 'b', 'LineWidth', 1.5);
plot(t, r, 'g', 'LineWidth', 1.5);
plot(t, -r, 'g--', 'LineWidth', 1.5); % 負包絡
legend('m_I(t)', 'm_Q(t)', 'r(t)');
xlabel('Time (s)');
ylabel('Amplitude');
title('In-Phase Component, Quadrature Component, and Envelope');
grid on;

%第八題的圖 m(t) r(t)
figure;
plot(t, m_t, 'r', 'LineWidth', 1.5); hold on;
plot(t, r, 'g', 'LineWidth', 1.5);
plot(t, -r, 'g--', 'LineWidth', 1.5); 
xlabel('Time (s)');
ylabel('Amplitude');
title('Bandpass signal and envelope');
grid on;

f = 20; 
a_t = m_I .* cos(2*pi*fc*t);
b_t = m_Q .* sin(2*pi*fc*t);

c_t = a_t - b_t;
%先畫他們三個的圖做檢查
figure;
plot(t,a_t,'r','LineWidth',1.5); hold on;
plot(t,b_t,'b','LineWidth',1.5); 
plot(t,c_t,'g','LineWidth',1.5); 

xlabel('Time (s)');
ylabel('Amplitude');
title('a(t) and b(t) and c(t)');
grid on;
%第九題的圖 c(t) m(t)
figure;
plot(t,c_t,'g','LineWidth',1.5); hold on;
plot(t,m_t,'r--','LineWidth',1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('c(t) and m(t)');
grid on;




