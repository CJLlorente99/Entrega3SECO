clear all; clc; close all;

%% Initiation Parameters
K = 1500;
p = 35;

x = p.*(0:0.0005/p:0.05);

aux = -4:1:1;
auxTi = -1.5:0.5:1;
auxTd = -4:1:1;

Kp = 10.^aux;
Ti = 10.^auxTi;
Td = 10.^auxTd;

Kpfija = 1;
Tdfija = 0.001;
Tifija = 1;

%% Td, Ti fijo

% Step
figure(1)
grid on;

Legend = cell(length(aux),1);

u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

for i = 1:length(Kp)
    num = [K*Kp(i) K*Kp(i)/Tifija];
    den = [1 p+K*Kp(i)*Tdfija K*Kp(i) K*Kp(i)/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PI-DEscalonKp.png')

% Ramp
figure(2)
grid on;

Legend = cell(length(aux),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Kp)
    num = [K*Kp(i) K*Kp(i)/Tifija];
    den = [1 p+K*Kp(i)*Tdfija K*Kp(i) K*Kp(i)/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a una rampa (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PI-DRampaKp.png')

% Parable
figure(3)
grid on;

Legend = cell(length(aux),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Kp)
    num = [K*Kp(i) K*Kp(i)/Tifija];
    den = [1 p+K*Kp(i)*Tdfija K*Kp(i) K*Kp(i)/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a una parabola (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PI-DParabolaKp.png')

%% Ti, Kp fijo

% Step
figure(4)
grid on;

Legend = cell(length(auxTd),1);

u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

for i = 1:length(Td)
    num = [K*Kpfija K*Kpfija/Tifija];
    den = [1 p+K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d} = ', num2str(Td(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Td)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PI-DEscalonTd.png')

% Ramp
figure(5)
grid on;

Legend = cell(length(auxTd),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Td)
    num = [K*Kpfija K*Kpfija/Tifija];
    den = [1 p+K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d} = ', num2str(Td(i)));
    hold on;
end

title('Comportamiento frente a una rampa (Td)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PI-DRampaTd.png')

% Parable
figure(6)
grid on;

Legend = cell(length(auxTd),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Td)
    num = [K*Kpfija K*Kpfija/Tifija];
    den = [1 p+K*Kpfija*Td(i) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d} = ', num2str(Td(i)));
    hold on;
end

title('Comportamiento frente a una parabola (Td)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PI-DParabolaTd.png')

%% Kp, Td fijo

% Step
figure(7)
grid on;

Legend = cell(length(auxTi),1);

u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

for i = 1:length(Ti)
    num = [K*Kpfija K*Kpfija/Ti(i)];
    den = [1 p+K*Kpfija*Tdfija K*Kpfija K*Kpfija/Ti(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{i} = ', num2str(Ti(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Ti)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PI-DEscalonTi.png')

% Ramp
figure(8)
grid on;

Legend = cell(length(auxTi),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Ti)
    num = [K*Kpfija K*Kpfija/Ti(i)];
    den = [1 p+K*Kpfija*Tdfija K*Kpfija K*Kpfija/Ti(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{i} = ', num2str(Ti(i)));
    hold on;
end

title('Comportamiento frente a una rampa (Ti)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PI-DRampaTi.png')

% Parable
figure(9)
grid on;

Legend = cell(length(auxTi),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Ti)
    num = [K*Kpfija K*Kpfija/Ti(i)];
    den = [1 p+K*Kpfija*Tdfija K*Kpfija K*Kpfija/Ti(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{i} = ', num2str(Td(i)));
    hold on;
end

title('Comportamiento frente a una parabola (Ti)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PI-DParabolaTi.png')
