clear all; clc;

K = 1500;
p = 35;

syms s;

f_step = 1/s;
f_ramp = 1/s^2;
f_parab = 1/s^3;

aux = -3:2:3;

Kp = 2.^aux;

figure(1)

plot(0:2.5,ones(1,length(0:2.5)), '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Kp)
    F = f_step*(Kp(i)*K)/(s^2 + p*s+ K*Kp(i));
    f = ilaplace(F);
    graph = fplot(f, [0 2.5]);
    set(graph, 'DisplayName', strcat('Kp = ', num2str(Kp(i))))
    hold on;
end

legend
hold off;

figure(2)

clear graph;
plot(0:2.5,0:2.5, '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Kp)
    F = f_ramp*(Kp(i)*K)/(s^2 + p*s+ K*Kp(i));
    f = ilaplace(F);
    graph = fplot(f, [0 2.5]);
    set(graph, 'DisplayName', strcat('Kp = ', num2str(Kp(i))))
    hold on;
end

legend
hold off;

figure(3)

clear graph;
plot(0:0.05:2.5,(0:0.05:2.5).^2, '--', 'DisplayName', 'Reference')
hold on;

for i = 1:length(Kp)
    F = f_parab*(Kp(i)*K)/(s^2 + p*s+ K*Kp(i));
    f = ilaplace(F);
    graph = fplot(f, [0 2.5]);
    set(graph, 'DisplayName', strcat('Kp = ', num2str(Kp(i))))
    hold on;
end

legend
hold off;