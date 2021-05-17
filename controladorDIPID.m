clear all; clc; close all;

K = 1500;
p = 35;

syms s;

f_step = 1/s;
f_ramp = 1/s^2;
f_parab = 1/s^3;

aux = -4:1:1;
auxTi = -1:1:0;
auxTd1 = -4:1:1;
auxTd2 = -4:1:1;

Kp = 10.^aux;
Ti = 10.^auxTi;
Td1 = 10.^auxTd1;
Td2 = 10.^auxTd2;

%% Td1, Td2, Ti fijo

figure(1)

plot(0:10,ones(1,length(0:10)), '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Kp)
    F = f_step*(K*Kp(i)*2*(s^2+s/2+1/2))/(s^2*(s+p) + K*Kp(i)*(s^2+s+1));
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
    F = f_ramp*(K*Kp(i)*2*(s^2+s/2+1/2))/(s^2*(s+p) + K*Kp(i)*(s^2+s+1));
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
    F = f_parab*(K*Kp(i)*2*(s^2+s/2+1/2))/(s^2*(s+p) + K*Kp(i)*(s^2+s+1));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('Kp = ', num2str(Kp(i))))
    hold on;
end

legend
hold off;

%% Ti, Td2, Kp fijo

figure(4)

plot(0:10,ones(1,length(0:10)), '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Td)
    F = f_step*(K*(1+Td1(i))*(s^2+s/(1+Td1(i))+1/(1+Td1(i))))/(s^2*(s+p)+K*Td1(i)*(s^2+s/Td1(i) + 1/Td(i)));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('\tau_{d1} = ', num2str(Td1(i))))
    hold on;
end

legend
hold off;

figure(5)

clear graph;
plot(0:10,0:10, '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Td)
    F = f_ramp*(K*(1+Td1(i))*(s^2+s/(1+Td1(i))+1/(1+Td1(i))))/(s^2*(s+p)+K*Td1(i)*(s^2+s/Td1(i) + 1/Td(i)));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('\tau_{d1} = ', num2str(Td1(i))))
    hold on;
end

legend
hold off;

figure(6)

clear graph;
plot(0:0.05:10,(0:0.05:10).^2, '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Td)
    F = f_parab*(K*(1+Td1(i))*(s^2+s/(1+Td1(i))+1/(1+Td1(i))))/(s^2*(s+p)+K*Td1(i)*(s^2+s/Td1(i) + 1/Td(i)));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('\tau_{d1} = ', num2str(Td1(i))))
    hold on;
end

legend
hold off;

%% Ti, Td1, Kp fijo

figure(7)

plot(0:10,ones(1,length(0:10)), '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Td)
    F = f_step*(K*(1+Td2(i))*(s^2 + s/(1+Td2(i)) + 1/(1+Td2(i))))/(s^2*(s+p)+K*(s^2+s + 1));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('\tau_{d2} = ', num2str(Td2(i))))
    hold on;
end

legend
hold off;

figure(8)

clear graph;
plot(0:10,0:10, '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Td)
    F = f_ramp*(K*(1+Td2(i))*(s^2 + s/(1+Td2(i)) + 1/(1+Td2(i))))/(s^2*(s+p)+K*(s^2+s + 1));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('\tau_{d2} = ', num2str(Td2(i))))
    hold on;
end

legend
hold off;

figure(9)

clear graph;
plot(0:0.05:10,(0:0.05:10).^2, '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Td)
    F = f_parab*(K*(1+Td2(i))*(s^2 + s/(1+Td2(i)) + 1/(1+Td2(i))))/(s^2*(s+p)+K*(s^2+s + 1));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('\tau_{d2} = ', num2str(Td2(i))))
    hold on;
end

legend
hold off;


%% Kp, Td1, Td2 fijo

figure(10)

plot(0:10,ones(1,length(0:10)), '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Ti)
    F = f_step*(K*(s^2+s/2+1/(2*Ti(i))))/(s^2*(s+p) + (K*(s^2+s+1/Ti(i))));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('\tau_{i} = ', num2str(Ti(i))))
    hold on;
end

legend
hold off;

figure(11)

clear graph;
plot(0:10,0:10, '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Ti)
    F = f_ramp*(K*(s^2+s/2+1/(2*Ti(i))))/(s^2*(s+p) + (K*(s^2+s+1/Ti(i))));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('\tau_{i} = ', num2str(Ti(i))))
    hold on;
end

legend
hold off;

figure(12)

clear graph;
plot(0:0.05:10,(0:0.05:10).^2, '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Ti)
    F = f_parab*(K*(s^2+s/2+1/(2*Ti(i))))/(s^2*(s+p) + (K*(s^2+s+1/Ti(i))));
    f = ilaplace(F);
    graph = fplot(f, [0 10]);
    set(graph, 'DisplayName', strcat('\tau_{i} = ', num2str(Ti(i))))
    hold on;
end

legend
hold off;