clear all; clc; close all;

K = 1500;
p = 35;

x = p.*(0:0.005/p:0.05);

auxKp = -4:1:1;
auxTi = -1:1:0;
auxTd1 = -4:1:1;
auxTd2 = -4:1:1;

Kp = 10.^auxKp;
Ti = 10.^auxTi;
Td1 = 10.^auxTd1;
Td2 = 10.^auxTd2;

Kpfija = 1;
Tifija = 1;
Td1fija = 1;
Td2fija = 1;

%% Kp libre

figure(1)
grid on;

Legend = cell(length(auxKp),1);

u = ones(1,length(x));
u(1:length(u)/6) = 0;
t = x;

for i = 1:length(Kp)
    num = K*Kp(i)*(Td1fija+Td2fija).*[1 1/(Td1fija+Td2fija) 1/((Td1fija+Td2fija)*Tifija)];
    den = [1 (p+K*Kp(i)*Td1fija) K*Kp(i) K*Kp(i)/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Kp)')
legend(Legend);
hold off;

figure(2)
grid on;

Legend = cell(length(auxKp),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:length(u)/6) = 0;

for i = 1:length(Kp)
    num = K*Kp(i)*(Td1fija+Td2fija).*[1 1/(Td1fija+Td2fija) 1/((Td1fija+Td2fija)*Tifija)];
    den = [1 (p+K*Kp(i)*Td1fija) K*Kp(i) K*Kp(i)/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un rampa (Kp)')
legend(Legend);
hold off;

figure(3)
grid on;

Legend = cell(length(auxKp),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:length(u)/6) = 0;

for i = 1:length(Kp)
    num = K*Kp(i)*(Td1fija+Td2fija).*[1 1/(Td1fija+Td2fija) 1/((Td1fija+Td2fija)*Tifija)];
    den = [1 (p+K*Kp(i)*Td1fija) K*Kp(i) K*Kp(i)/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un parabola (Kp)')
legend(Legend);
hold off;

%% Ti libre

figure(4)
grid on;

Legend = cell(length(auxTi),1);

u = ones(1,length(x));
u(1:length(u)/6) = 0;
t = x;

for i = 1:length(Ti)
    num = K*Kpfija*(Td1fija+Td2fija).*[1 1/(Td1fija+Td2fija) 1/((Td1fija+Td2fija)*Ti(i))];
    den = [1 (p+K*Kpfija*Td1fija) K*Kpfija K*Kpfija/Ti(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Ti)')
legend(Legend);
hold off;

figure(5)
grid on;

Legend = cell(length(auxTi),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:length(u)/6) = 0;

for i = 1:length(Ti)
    num = K*Kpfija*(Td1fija+Td2fija).*[1 1/(Td1fija+Td2fija) 1/((Td1fija+Td2fija)*Ti(i))];
    den = [1 (p+K*Kpfija*Td1fija) K*Kpfija K*Kpfija/Ti(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un rampa (Ti)')
legend(Legend);
hold off;

figure(6)
grid on;

Legend = cell(length(auxTi),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:length(u)/6) = 0;

for i = 1:length(Ti)
    num = K*Kpfija*(Td1fija+Td2fija).*[1 1/(Td1fija+Td2fija) 1/((Td1fija+Td2fija)*Ti(i))];
    den = [1 (p+K*Kpfija*Td1fija) K*Kpfija K*Kpfija/Ti(i)];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un parabola (Ti)')
legend(Legend);
hold off;

%% Td1 libre

figure(7)
grid on;

Legend = cell(length(auxTd1),1);

u = ones(1,length(x));
u(1:length(u)/6) = 0;
t = x;

for i = 1:length(Td1)
    num = K*Kpfija*(Td1(i)+Td2fija).*[1 1/(Td1(i)+Td2fija) 1/((Td1(i)+Td2fija)*Tifija)];
    den = [1 (p+K*Kpfija*Td1(i)) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Td1)')
legend(Legend);
hold off;

figure(8)
grid on;

Legend = cell(length(auxTd1),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:length(u)/6) = 0;

for i = 1:length(Td1)
    num = K*Kpfija*(Td1(i)+Td2fija).*[1 1/(Td1(i)+Td2fija) 1/((Td1(i)+Td2fija)*Tifija)];
    den = [1 (p+K*Kpfija*Td1(i)) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un rampa (Td1)')
legend(Legend);
hold off;

figure(9)
grid on;

Legend = cell(length(auxTd1),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:length(u)/6) = 0;

for i = 1:length(Td1)
    num = K*Kpfija*(Td1(i)+Td2fija).*[1 1/(Td1(i)+Td2fija) 1/((Td1(i)+Td2fija)*Tifija)];
    den = [1 (p+K*Kpfija*Td1(i)) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un parabola (Td1)')
legend(Legend);
hold off;

%% Td2 libre

figure(10)
grid on;

Legend = cell(length(auxTd2),1);

u = ones(1,length(x));
u(1:length(u)/6) = 0;
t = x;

for i = 1:length(Td2)
    num = K*Kpfija*(Td1fija+Td2(i)).*[1 1/(Td1fija+Td2(i)) 1/((Td1fija+Td2(i))*Tifija)];
    den = [1 (p+K*Kpfija*Td1fija) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un escalon (Td2)')
legend(Legend);
hold off;

figure(11)
grid on;

Legend = cell(length(auxTd2),1);

t = x;
u = t-t(ceil(length(t)/6));
u(1:length(u)/6) = 0;

for i = 1:length(Td2)
    num = K*Kpfija*(Td1fija+Td2(i)).*[1 1/(Td1fija+Td2(i)) 1/((Td1fija+Td2(i))*Tifija)];
    den = [1 (p+K*Kpfija*Td1fija) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un rampa (Td2)')
legend(Legend);
hold off;

figure(12)
grid on;

Legend = cell(length(auxTd2),1);

t = x;
u = (t-t(ceil(length(t)/6))).^2;
u(1:length(u)/6) = 0;

for i = 1:length(Td2)
    num = K*Kpfija*(Td1fija+Td2(i)).*[1 1/(Td1fija+Td2(i)) 1/((Td1fija+Td2(i))*Tifija)];
    den = [1 (p+K*Kpfija*Td1fija) K*Kpfija K*Kpfija/Tifija];
    sys = tf(num,den);
    lsim(sys, u, t)
    Legend{i} = strcat('K_{P} = ', num2str(Kp(i)));
    hold on;
end

title('Comportamiento frente a un parabola (Td2)')
legend(Legend);
hold off;

%% Rlocus
% No idea de si se puede despejar la ecuación de lazo abierto