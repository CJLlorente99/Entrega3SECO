clear all; clc; close all;

%% Initiation Parameters
K = 1500;
p = 35;

x = p.*(0:0.0005/p:0.05);

aux = -4:1:1;
auxTi = -1.5:0.5:1;
auxTd1 = -4:1:1;
auxTd2 = -4:1:1;

Kp = 10.^aux;
Ti = 10.^auxTi;
Td1 = 10.^auxTd1;
Td2 = 10.^auxTd2;

Kpfija = 1;
Td1fija = 0.001;
Td2fija = 0.001;
Tifija = 1;

%% Rlocus
% Td = Td1fija + Td2fija;

% for i = 1:length(auxTd2)
%     figure (i)
%     
%     num = [K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
%     den = [1 p+K*Kpfija*Td1fija K*Kpfija K*Kpfija/Tifija];
%     
%     rlocus(num,den)
% end

%% Td1, Td2, Ti fijo
Td = Td1fija + Td2fija;

% Step
figure(1)
grid on;

Legend = cell(length(aux),1);

u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

for i = 1:length(Kp)
    num = [K*Kp(i)*Td K*Kp(i) K*Kp(i)/Tifija];
    den = [1 p+K*Kp(i)*Td1fija K*Kp(i) K*Kp(i)/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{p} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_DIPIDEscalonKp.png')

% Ramp
figure(2)
grid on;

Legend = cell(length(aux),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Kp)
    num = [K*Kp(i)*Td K*Kp(i) K*Kp(i)/Tifija];
    den = [1 p+K*Kp(i)*Td1fija K*Kp(i) K*Kp(i)/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{p} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a una rampa (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_DIPIDRampaKp.png')

% Parable
figure(3)
grid on;

Legend = cell(length(aux),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Kp)
    num = [K*Kp(i)*Td K*Kp(i) K*Kp(i)/Tifija];
    den = [1 p+K*Kp(i)*Td1fija K*Kp(i) K*Kp(i)/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{p} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a una parabola (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_DIPIDParabolaKp.png')

%% Ti, Td2, Kp fijo
Td = Td1 + Td2fija;

% Step
figure(4)
grid on;

Legend = cell(length(auxTd1),1);

u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

for i = 1:length(Td1)
    num = [K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
    den = [1 p+K*Kpfija*Td1(i) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d1} = ', num2str(Td1(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Td1)')
legend(Legend);
hold off;
saveas(gcf,'img/01_DIPIDEscalonTd1.png')

% Ramp
figure(5)
grid on;

Legend = cell(length(auxTd1),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Td1)
    num = [K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
    den = [1 p+K*Kpfija*Td1(i) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d1} = ', num2str(Td1(i)));
    hold on;
end

title('Comportamiento frente a una rampa (Td1)')
legend(Legend);
hold off;
saveas(gcf,'img/01_DIPIDRampaTd1.png')

% Parable
figure(6)
grid on;

Legend = cell(length(auxTd1),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Td1)
    num = [K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
    den = [1 p+K*Kpfija*Td1(i) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d1} = ', num2str(Td1(i)));
    hold on;
end

title('Comportamiento frente a una parabola (Td1)')
legend(Legend);
hold off;
saveas(gcf,'img/01_DIPIDParabolaTd1.png')

%% Ti, Td1, Kp fijo
Td = Td1fija + Td2;

% Step
figure(7)
grid on;

Legend = cell(length(auxTd2),1);

u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

for i = 1:length(Td2)
    num = [K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
    den = [1 p+K*Kpfija*Td1fija K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d2} = ', num2str(Td2(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Td2)')
legend(Legend);
hold off;
saveas(gcf,'img/01_DIPIDEscalonTd2.png')

% Ramp
figure(8)
grid on;

Legend = cell(length(auxTd2),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Td2)
    num = [K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
    den = [1 p+K*Kpfija*Td1fija K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d2} = ', num2str(Td2(i)));
    hold on;
end

title('Comportamiento frente a una rampa (Td2)')
legend(Legend);
hold off;
saveas(gcf,'img/01_DIPIDRampaTd2.png')

% Parable
figure(9)
grid on;

Legend = cell(length(auxTd2),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Td2)
    num = [K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
    den = [1 p+K*Kpfija*Td1fija K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d2} = ', num2str(Td2(i)));
    hold on;
end

title('Comportamiento frente a una parabola (Td2)')
legend(Legend);
hold off;
saveas(gcf,'img/01_DIPIDParabolaTd2.png')


%% Kp, Td1, Td2 fijo
Td = Td1fija + Td2fija;

% Step
figure(10)
grid on;

Legend = cell(length(auxTi),1);

u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

for i = 1:length(Ti)
    num = [K*Kpfija*Td K*Kpfija K*Kpfija/Ti(i)];
    den = [1 p+K*Kpfija*Td1fija K*Kpfija K*Kpfija/Ti(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{i} = ', num2str(Ti(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Ti)')
legend(Legend);
hold off;
saveas(gcf,'img/01_DIPIDEscalonTi.png')

% Ramp
figure(11)
grid on;

Legend = cell(length(auxTi),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Ti)
    num = [K*Kpfija*Td K*Kpfija K*Kpfija/Ti(i)];
    den = [1 p+K*Kpfija*Td1fija K*Kpfija K*Kpfija/Ti(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{i} = ', num2str(Ti(i)));
    hold on;
end

title('Comportamiento frente a una rampa (Ti)')
legend(Legend);
hold off;
saveas(gcf,'img/01_DIPIDRampaTi.png')

% Parable
figure(12)
grid on;

Legend = cell(length(auxTi),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Ti)
    num = [K*Kpfija*Td K*Kpfija K*Kpfija/Ti(i)];
    den = [1 p+K*Kpfija*Td1fija K*Kpfija K*Kpfija/Ti(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{i} = ', num2str(Ti(i)));
    hold on;
end

title('Comportamiento frente a una parabola (Ti)')
legend(Legend);
hold off;
saveas(gcf,'img/01_DIPIDParabolaTi.png')

%% Derivada en 0

figure(10)

Kpfija = 0.15;
Td1fija = -0.001;
Td2fija = -0.05;
Tifija = 1;
Td = Td1fija + Td2fija;

u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

num = [K*Kpfija*Td K*Kpfija K*Kpfija/Tifija];
den = [1 p+K*Kpfija*Td1fija K*Kpfija K*Kpfija/Tifija];
sys = tf(num,den);
lsim(sys, u, t)

title('Comportamiento frente a un escalón para Td1 = -0.05, Td2 = 0.001, Ti = 1 y Kp = 0.15');

saveas(gcf,'img/01_DIPIDDerivada.png')

figure(11)

num = K.*[Td 1 1/Tifija];
den = [1 (p-K*Td2fija) 0 0];
rlocus(num,den)

saveas(gcf,'img/01_DIPIDRlocusDerivada.png')