%% Constantes

reductora = 23;
p = 64.986;
K = 2652.28*reductora;

%% Creacion de señales de referencia

x = p.*(0:0.0001/p:1/6000);
t = x;
nInit = ceil(length(x)/6);

% Escalon

uEscalon = ones(1,length(x));
uEscalon(1:nInit) = 0;

% Rampa

uRampa = t-t(nInit);
uRampa(1:nInit) = 0;

% Parabola

uParabola = (t-t(nInit)).^2;
uParabola(1:nInit) = 0;

%% Experimento 1

Kp = 50.9597; Kd1 = 0.063379; Kd2 = -0.0010653; Ki = 0.036504;
Td1 = Kd1/Kp; Td2 = Kd2/Kp; Ti = Kp/Ki; Td = Td1 + Td2;
num = [K*Kp*Td1 K*Kp K*Kp/Ti];den = [1 p+K*Kp*Td K*Kp K*Kp/Ti];
sys = tf(num,den);

% Escalon

figure(1)

plot(t(nInit).*[1 1], [-50 50],'c');
lsim(sys, uEscalon, t);
saveas(gcf,'img/respuestaIdealEscalonKp50.png')

% Rampa

figure(2)

plot(t(nInit).*[1 1], [-50 50],'c');
lsim(sys, uRampa, t);
saveas(gcf,'img/respuestaIdealRampaKp50.png')

% Parabola

figure(3)

plot(t(nInit).*[1 1], [-50 50],'c');
lsim(sys, uParabola, t);
saveas(gcf,'img/respuestaIdealParabolaKp50.png')

%% Creacion de señales de referencia

x = p.*(0:0.0001/p:1/1500);
t = x;
nInit = ceil(length(x)/6);

% Escalon

uEscalon = ones(1,length(x));
uEscalon(1:nInit) = 0;

% Rampa

uRampa = t-t(nInit);
uRampa(1:nInit) = 0;

% Parabola

uParabola = (t-t(nInit)).^2;
uParabola(1:nInit) = 0;

%% Experimento 2

Kp = 4.2558; Kd1 = 0.027703; Kd2 = -0.0010653; Ki = 750.4625;
Td1 = Kd1/Kp; Td2 = Kd2/Kp; Ti = Kp/Ki; Td = Td1 + Td2;
num = [K*Kp*Td1 K*Kp K*Kp/Ti];den = [1 p+K*Kp*Td K*Kp K*Kp/Ti];
sys = tf(num,den);

% Escalon

figure(4)

plot(t(nInit).*[1 1], [-50 50],'c');
lsim(sys, uEscalon, t);
saveas(gcf,'img/respuestaIdealEscalonKp4.png')

% Rampa

figure(5)

plot(t(nInit).*[1 1], [-50 50],'c');
lsim(sys, uRampa, t);
saveas(gcf,'img/respuestaIdealRampaKp4.png')

% Parabola

figure(6)

plot(t(nInit).*[1 1], [-50 50],'c');
lsim(sys, uParabola, t);
saveas(gcf,'img/respuestaIdealParabolaKp4.png')