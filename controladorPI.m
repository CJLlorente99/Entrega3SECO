clear all; clc; close all;

%% Initiation Parameters
K = 1500;
p = 35;

x = p.*(0:0.0005/p:0.1);

aux = -4:1:1;
auxTi = -1.5:0.5:1;

Kp = 10.^aux;
Ti = 10.^auxTi;

Kpfija = 10;
Tifija = 1;

%% Ti fijo

% Step
figure(1)
grid on;

Legend = cell(length(aux),1);

u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

for i = 1:length(Kp)
    num = [K*Kp(i) K*Kp(i)/Tifija];
    den = [1 p K*Kp(i) K*Kp(i)/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PIEscalonKp.png')

% Ramp
figure(2)
grid on;

Legend = cell(length(aux),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Kp)
    num = [K*Kp(i) K*Kp(i)/Tifija];
    den = [1 p K*Kp(i) K*Kp(i)/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un rampa (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PIRampaKp.png')

% Parable
figure(3)
grid on;

Legend = cell(length(aux),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Kp)
    num = [K*Kp(i) K*Kp(i)/Tifija];
    den = [1 p K*Kp(i) K*Kp(i)/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un parabola (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PIParabolaKp.png')

%% Kp fijo

% Step
figure(4);
grid on;

Legend = cell(length(auxTi),1);

u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

for i = 1:length(auxTi)
    num = [K*Kpfija K*Kpfija/Ti(i)];
    den = [1 p K*Kpfija K*Kpfija/Ti(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('\tau_{i} = ', num2str(Ti(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Ti)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PIEscalonTi.png')

% Ramp
figure(5)
grid on;

Legend = cell(length(auxTi),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(auxTi)
    num = [K*Kpfija K*Kpfija/Ti(i)];
    den = [1 p K*Kpfija K*Kpfija/Ti(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('\tau_{i} = ', num2str(Ti(i)));
    hold on;
end

title('Comportamiento frente a un rampa (Ti)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PIRampaTi.png')

% Parable
figure(6)
grid on;

Legend = cell(length(auxTi),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(auxTi)
    num = [K*Kpfija K*Kpfija/Ti(i)];
    den = [1 p K*Kpfija K*Kpfija/Ti(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('\tau_{i} = ', num2str(Ti(i)));
    hold on;
end

title('Comportamiento frente a un parabola (Ti)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PIParabolaTi.png')

%% Rlocus

figure(7)

num = [K*Tifija K];
den = [Tifija Tifija*p 0 0];
rlocus(num,den)
saveas(gcf,'img/01_PIRLocusKp.png')

% fig = 8;
% 
% for i = 1:length(auxTi)
%     figure(fig)
%     num = [K*Ti(i) K];
%     den = [Ti(i) Ti(i)*p 0 0];
%     rlocus(num,den)
%     fig = fig + 1;
% end