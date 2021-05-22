clear all; close all;

%% Initiation Parameters
K = 1500;
p = 35;

x = p.*(0:0.0001/p:0.05);

aux = -4:1:1;
auxTd = -4:1:1;

Kp = 10.^aux;
Td = 10.^aux;

Kpfija = 1;
Tdfija = 0.001;

%% Td fijo

% Step
figure(1)
grid on;

Legend = cell(length(aux),1);

u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

for i = 1:length(Kp)
    num = [Kp(i)*K*Tdfija Kp(i)*K];
    den = [1 p+Kp(i)*K*Tdfija K*Kp(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PDEscalonKp.png')

% Ramp
figure(2)
grid on;

Legend = cell(length(aux),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Kp)
    num = [Kp(i)*K*Tdfija Kp(i)*K];
    den = [1 p+Kp(i)*K*Tdfija K*Kp(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un rampa (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PDRampaKp.png')

% Parable
figure(3)
grid on;

Legend = cell(length(aux),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Kp)
    num = [Kp(i)*K*Tdfija Kp(i)*K];
    den = [1 p+Kp(i)*K*Tdfija K*Kp(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un parabola (Kp)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PDParabolaKp.png')

%% Kp fijo

% Step
figure(4)
grid on;

Legend = cell(length(auxTd),1);

u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

for i = 1:length(Td)
    num = [Kpfija*K*Td(i) Kpfija*K];
    den = [1 p+Kpfija*K*Td(i) K*Kpfija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d} = ', num2str(Td(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Td)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PDEscalonTd.png')

% Ramp
figure(5)
grid on;

Legend = cell(length(auxTd),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Td)
    num = [Kpfija*K*Td(i) Kpfija*K];
    den = [1 p+Kpfija*K*Td(i) K*Kpfija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d} = ', num2str(Td(i)));
    hold on;
end

title('Comportamiento frente a un rampa (Td)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PDRampaTd.png')

% Parable
figure(6)
grid on;

Legend = cell(length(auxTd),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:ceil(length(u)/6)) = 0;

for i = 1:length(Td)
    num = [Kpfija*K*Td(i) Kpfija*K];
    den = [1 p+Kpfija*K*Td(i) K*Kpfija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('T_{d} = ', num2str(Td(i)));
    hold on;
end

title('Comportamiento frente a un parabola (Td)')
legend(Legend);
hold off;
saveas(gcf,'img/01_PDParabolaTd.png')

%% Rlocus

figure(7)

num = [K*Kpfija*Tdfija K*Kpfija];
den = [1 p 0];
rlocus(num,den)

saveas(gcf,'img/01_PDRlocusKp.png')