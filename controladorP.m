clear all; clc; close all;

K = 1500;
p = 35;

x = p.*(0:0.005/p:0.05);

aux = -3:2:3;

Kp = 2.^aux;

figure(1)
grid on;

Legend = cell(length(aux),1);

u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

for i = 1:length(Kp)
    num = Kp(i)*K;
    den = [1 p K*Kp(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PEscalon.png')

figure(2)
grid on;

Legend = cell(length(aux),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Kp)
    num = Kp(i)*K;
    den = [1 p K*Kp(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un rampa (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PRampa.png')

figure(3)
grid on;

Legend = cell(length(aux),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Kp)
    num = Kp(i)*K;
    den = [1 p K*Kp(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un parabola (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PParabola.png')

%% Rlocus

figure(4)

num = [K];
den = [1 p 0];

rlocus(num,den)

saveas(gcf,'img/01_PRlocus.png')