clear all; clc; close all;

K = 1500;
p = 35;

syms s;

f_step = 1/s;
f_ramp = 1/s^2;
f_parab = 1/s^3;

aux = -4:1:1;
auxTi = -1:1:0;
auxTd = -4:1:1;

Kp = 10.^aux;
Ti = 10.^auxTi;
Td = 10.^auxTd;

%% Td, Ti fijo

figure(1)

plot(0:10,ones(1,length(0:10)), '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Kp)
    F = f_step*(K*Kp(i)*(s+1))/(s^2*(s+p) + K*Kp(i)*(s^2+s+1));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('Kp = ', num2str(Kp(i))))
    hold on;
end

legend
hold off;

figure(2)

clear graph;
plot(0:10,0:10, '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Kp)
    F = f_ramp*(K*Kp(i)*(s+1))/(s^2*(s+p) + K*Kp(i)*(s^2+s+1));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('Kp = ', num2str(Kp(i))))
    hold on;
end

legend
hold off;

figure(3)

clear graph;
plot(0:0.05:10,(0:0.05:10).^2, '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Kp)
    F = f_parab*(K*Kp(i)*(s+1))/(s^2*(s+p) + K*Kp(i)*(s^2+s+1));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('Kp = ', num2str(Kp(i))))
    hold on;
end

legend
hold off;

%% Ti, Kp fijo

figure(4)

plot(0:10,ones(1,length(0:10)), '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Td)
    F = f_step*(K*(s + 1))/(s^2*(s+p)+K*Td(i)*(s^2+s/Td(i) + 1/Td(i)));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('\tau_{d} = ', num2str(Td(i))))
    hold on;
end

legend
hold off;

figure(5)

clear graph;
plot(0:10,0:10, '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Td)
    F = f_ramp*(K*(s+ 1))/(s^2*(s+p)+K*Td(i)*(s^2+s/Td(i) + 1/Td(i)));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('\tau_{d} = ', num2str(Td(i))))
    hold on;
end

legend
hold off;

figure(6)

clear graph;
plot(0:0.05:10,(0:0.05:10).^2, '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Td)
    F = f_parab*(K*(s + 1))/(s^2*(s+p)+K*Td(i)*(s^2+s/Td(i) + 1/Td(i)));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('\tau_{d} = ', num2str(Td(i))))
    hold on;
end

legend
hold off;

%% Kp. Td fijo

figure(7)

plot(0:10,ones(1,length(0:10)), '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Ti)
    F = f_step*(K*(s+1/Ti(i)))/(s^2*(s+p) + (K*(s^2+s+1/Ti(i))));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('\tau_{i} = ', num2str(Ti(i))))
    hold on;
end

legend
hold off;

figure(8)

clear graph;
plot(0:10,0:10, '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Ti)
    F = f_ramp*(K*(s+1/Ti(i)))/(s^2*(s+p) + (K*(s^2+s+1/Ti(i))));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('\tau_{i} = ', num2str(Ti(i))))
    hold on;
end

legend
hold off;

figure(9)

clear graph;
plot(0:0.05:10,(0:0.05:10).^2, '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Ti)
    F = f_parab*(K*(s+1/Ti(i)))/(s^2*(s+p) + (K*(s^2+s+1/Ti(i))));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('\tau_{i} = ', num2str(Ti(i))))
    hold on;
end

legend
hold off;