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
    num = [K*Kp(i)*Td1fija K*Kp(i) K*Kp(i)/Tifija];
    den = [1 p+K*Kp(i)*Td K*Kp(i) K*Kp(i)/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{p} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PID-DEscalonKp.png')

% Ramp
figure(2)
grid on;

Legend = cell(length(aux),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Kp)
    num = [K*Kp(i)*Td1fija K*Kp(i) K*Kp(i)/Tifija];
    den = [1 p+K*Kp(i)*Td K*Kp(i) K*Kp(i)/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{p} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a una rampa (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PID-DRampaKp.png')

% Parable
figure(3)
grid on;

Legend = cell(length(aux),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Kp)
    num = [K*Kp(i)*Td1fija K*Kp(i) K*Kp(i)/Tifija];
    den = [1 p+K*Kp(i)*Td K*Kp(i) K*Kp(i)/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{p} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a una parabola (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PID-DParabolaKp.png')

%% Ti, Td2, Kp fijo
Td = Td1 + Td2fija;

% Step
figure(4)
grid on;

Legend = cell(length(auxTd1),1);

u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

%     num = [K*Kpfija*Td1fija K*Kpfija K*Kpfija/Tifija];
%     den = [1 p+K*Kpfija*Td K*Kpfija K*Kpfija/Tifija];

for i = 1:length(Td1)
    num = [K*Kpfija*Td1(i) K*Kpfija K*Kpfija/Tifija];
    den = [1 p+K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d1} = ', num2str(Td1(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Td1)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PID-DEscalonTd1.png')

% Ramp
figure(5)
grid on;

Legend = cell(length(auxTd1),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Td1)
    num = [K*Kpfija*Td1(i) K*Kpfija K*Kpfija/Tifija];
    den = [1 p+K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d1} = ', num2str(Td1(i)));
    hold on;
end

title('Comportamiento frente a una rampa (Td1)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PID-DRampaTd1.png')

% Parable
figure(6)
grid on;

Legend = cell(length(auxTd1),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Td1)
    num = [K*Kpfija*Td1(i) K*Kpfija K*Kpfija/Tifija];
    den = [1 p+K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d1} = ', num2str(Td1(i)));
    hold on;
end

title('Comportamiento frente a una parabola (Td1)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PID-DParabolaTd1.png')

%% Ti, Td1, Kp fijo
Td = Td1fija + Td2;

% Step
figure(7)
grid on;

Legend = cell(length(auxTd2),1);

u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

%     num = [K*Kpfija*Td1fija K*Kpfija K*Kpfija/Tifija];
%     den = [1 p+K*Kpfija*Td K*Kpfija K*Kpfija/Tifija];

for i = 1:length(Td2)
    num = [K*Kpfija*Td1fija K*Kpfija K*Kpfija/Tifija];
    den = [1 p+K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d2} = ', num2str(Td2(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Td2)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PID-DEscalonTd2.png')

% Ramp
figure(8)
grid on;

Legend = cell(length(auxTd2),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Td2)
    num = [K*Kpfija*Td1fija K*Kpfija K*Kpfija/Tifija];
    den = [1 p+K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d2} = ', num2str(Td2(i)));
    hold on;
end

title('Comportamiento frente a una rampa (Td2)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PID-DRampaTd2.png')

% Parable
figure(9)
grid on;

Legend = cell(length(auxTd2),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Td2)
    num = [K*Kpfija*Td1fija K*Kpfija K*Kpfija/Tifija];
    den = [1 p+K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d2} = ', num2str(Td2(i)));
    hold on;
end

title('Comportamiento frente a una parabola (Td2)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PID-DParabolaTd2.png')


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
    num = [K*Kpfija*Td1fija K*Kpfija K*Kpfija/Ti(i)];
    den = [1 p+K*Kpfija*Td K*Kpfija K*Kpfija/Ti(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{i} = ', num2str(Ti(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Ti)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PID-DEscalonTi.png')

% Ramp
figure(11)
grid on;

Legend = cell(length(auxTi),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Ti)
    num = [K*Kpfija*Td1fija K*Kpfija K*Kpfija/Ti(i)];
    den = [1 p+K*Kpfija*Td K*Kpfija K*Kpfija/Ti(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{i} = ', num2str(Ti(i)));
    hold on;
end

title('Comportamiento frente a una rampa (Ti)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PID-DRampaTi.png')

% Parable
figure(12)
grid on;

Legend = cell(length(auxTi),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Ti)
    num = [K*Kpfija*Td1fija K*Kpfija K*Kpfija/Ti(i)];
    den = [1 p+K*Kpfija*Td K*Kpfija K*Kpfija/Ti(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{i} = ', num2str(Ti(i)));
    hold on;
end

title('Comportamiento frente a una parabola (Ti)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PID-DParabolaTi.png')