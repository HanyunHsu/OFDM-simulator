fs = 1000;
t = 0:1/fs:1; 
f = [16, 18, 22, 24]; 
a = [1, -1, 1, -1]; 
b = [1, 1, -1, 1];
m_t = zeros(size(t)); 
fc = 20;
fn = [-4,-2,2,4]; %bass
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

m_It = zeros(size(t)); 
m_Qt = zeros(size(t)); 
for n = 1:length(fn)
    m_It= m_It + (a(n)*cos(2*pi*fn(n)*t) - b(n)*sin(2*pi*fn(n)*t));
    m_Qt= m_Qt+(a(n)*sin(2*pi*fn(n)*t) + b(n)*cos(2*pi*fn(n)*t));
end

r = sqrt(m_It.^2 +m_Qt.^2);

%第七題的圖 r(t) m_I(t) m_Q(t)
figure;
plot(t, m_It, 'r', 'LineWidth', 1.5); hold on;
plot(t, m_Qt, 'b', 'LineWidth', 1.5);
plot(t, r, 'g', 'LineWidth', 1.5);
plot(t, -r, 'g--', 'LineWidth', 1.5); % 負包絡
legend('m_I(t)', 'm_Q(t)', 'r(t)');
xlabel('Time (s)');
ylabel('Amplitude');
title('m_I(t), m_Q(t) and r(t)');
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



a_t = m_It.* cos(2*pi*fc*t);
b_t = m_Qt.* sin(2*pi*fc*t);
c_t = a_t - b_t;
%第九題的圖 c(t) m(t)
figure;
plot(t,c_t,'g','LineWidth',1.5); hold on;
plot(t,m_t,'r--','LineWidth',1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('c(t) and m(t)');
grid on;




