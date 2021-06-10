close all; clear all;

%% Inicializar valores
reductora = 23;
p = 64.986;
K = 2652.28*reductora;

minResol = 0.02;

experimentosKp = 25;
Kpmin = -0.75;
% Kpmin = minResol;
Kpmax = 0;
% Kpmax = 100;

experimentosTd1 = 25;
Td1min = -4.6;
% Td1min = -2;
Td1max = -2;

experimentosTi = 25;
Timin = -8;
% Timin = 0;
Timax = -4.6;

% [kps, Td1s, Tis] = ndgrid(  linspace(Kpmin, Kpmax, experimentosKp),...
%                             linspace(Td1min, Td1max, experimentosTd1),...
%                             linspace(Timin, Timax, experimentosTi));
                        
[kps, Td1s, Tis] = ndgrid(  exp(linspace(Kpmin, Kpmax, experimentosKp)),...
                            exp(linspace(Td1min, Td1max, experimentosTd1)),...
                            exp(linspace(Timin, Timax, experimentosTi)));
                        
n = 1;

% Valores correctos

kpCorrecto = [];
tiCorrecto = [];
tdCorrecto = [];
    
% Init somethings

x = p.*(0:0.0001/p:0.05);
u = ones(1,length(x));
t = x;

% figure(1)

Mpmin = 0.08;
% Mpmin = 0;
% Mpmax = 0.15;
Mpmax = 0.15;

tsmax = 0.5; %segundos
tolerancia = 0.02;

trmax = 0.3; %segundos
trmin = 0; %segundos
            
for kpi = 1: experimentosKp
    for td1i = 1: experimentosTd1
        for tii = 1: experimentosTi
    
%           disp(['Kp ' num2str(Kp) ' Td1 ' num2str(Td1) ' Td2 ' num2str(Td2) ' Ti ' num2str(Ti)]);

            Kp = kps(kpi,td1i,tii);
            Td1 = Td1s(kpi,td1i,tii);
            Ti = Tis(kpi,td1i,tii);

            disp([num2str(n) '/' num2str(experimentosKp*experimentosTd1*experimentosTi)]);
%           %disp(['Kp ' num2str(Kp) ' Td1 ' num2str(Td1) ' Td2 ' num2str(Td2) ' Ti ' num2str(Ti)]);
            n = n + 1;
        %% Calculo de parametros

            Td2 = -p/(K*Kp);
            Td = Td1+Td2;

            num = [K*Kp*Td1 K*Kp K*Kp/Ti];
            den = [1 p+K*Kp*Td K*Kp K*Kp/Ti];
            sys = tf(num,den);

            %% Condiciones de estabilidad
            epsilon = p+K*Kp*Td;

            if Kp < 0
%                 disp(['Kp < 0 ' 'Kp ' num2str(Kp) ' Td1 ' num2str(Td1) ' Td2 ' num2str(Td2) ' Ti ' num2str(Ti)]);
                continue
                % no estable
            end

            if Ti < 0
%                 disp(['Ti < 0 ' 'Kp ' num2str(Kp) ' Td1 ' num2str(Td1) ' Td2 ' num2str(Td2) ' Ti ' num2str(Ti)]);
                continue
                % no estable
            end

            if Kp*K*(1-1/(Ti*epsilon)) < 0
%                 disp(['Kp*K*(1-1/(Ti*epsilon)) < 0 ' 'Kp ' num2str(Kp) ' Td1 ' num2str(Td1) ' Td2 ' num2str(Td2) ' Ti ' num2str(Ti)]);
                continue
                % no estable
            end

            if epsilon <= 0
%                 disp(['epsilon <= 0 ' 'Kp ' num2str(Kp) ' Td1 ' num2str(Td1) ' Td2 ' num2str(Td2) ' Ti ' num2str(Ti) ' epsilon ' num2str(epsilon)]);
                continue
                % no estable
            end

            %% Respuesta ante el escalon

            y = lsim(sys, u, t);

            %% Condicion de sobreelongacion maxima y minima

            if sum(y > (1 + Mpmax)) >= 1
%                 disp(['Cond. sobreelongacion max ' 'Kp ' num2str(Kp) ' Td1 ' num2str(Td1) ' Td2 ' num2str(Td2) ' Ti ' num2str(Ti)]);
                continue
            end

            if sum(y > (1 + Mpmin)) == 0
%                 disp(['Cond. sobreelongacion min ' 'Kp ' num2str(Kp) ' Td1 ' num2str(Td1) ' Td2 ' num2str(Td2) ' Ti ' num2str(Ti)]);
                continue
            end

            % Condicion de tiempo de establecimiento

            if sum(y(find(t >= tsmax,1):end) > (1 + tolerancia)) >= 1
%                 disp(['Cond. Ts max ' 'Kp ' num2str(Kp) ' Td1 ' num2str(Td1) ' Td2 ' num2str(Td2) ' Ti ' num2str(Ti)]);
                continue
            end

            if sum(y(find(t >= tsmax,1):end) < (1 - tolerancia)) >= 1
%                 disp(['Cond. Ts min ' 'Kp ' num2str(Kp) ' Td1 ' num2str(Td1) ' Td2 ' num2str(Td2) ' Ti ' num2str(Ti)]);
                continue
            end

            %% Condicion de tiempo de subida

            if sum(t(find(y >= 1,1)) >= trmax) >= 1
%                 disp(['Cond. tiempo subida max ' 'Kp ' num2str(Kp) ' Td1 ' num2str(Td1) ' Td2 ' num2str(Td2) ' Ti ' num2str(Ti)]);
                continue
            end

            if sum(t(find(y >= 1,1)) <= trmin) >= 1
%                 disp(['Cond. tiempo subida min ' 'Kp ' num2str(Kp) ' Td1 ' num2str(Td1) ' Td2 ' num2str(Td2) ' Ti ' num2str(Ti)]);
                continue
            end
            
            %% Condiciones resolucion
            
            if abs(Kp) < minResol
                continue
            end
            
            if abs(Td1*Kp) < minResol
                continue
            end
            
%             if abs(Td2*Kp) < minResol
%                 continue
%             end
            
            if abs(Kp/Ti) < minResol
                continue
            end

            if abs(Kp) > 999
                continue
            end
            
            if abs(Td1*Kp) > 999
                continue
            end
            
            if abs(Td2*Kp) > 999
                continue
            end
            
            if abs(Kp/Ti) > 999
                continue
            end

            %% Plot
            info = ['Kp ' num2str(Kp) ' Kd1 ' num2str(Kp*Td1) ' Kd2 ' num2str(Kp*Td2) ' Ki ' num2str(Kp/Ti)];
            h = plot(t,y, 'DisplayName', info);
            h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Trace',repmat({h.DisplayName},size(h.XData)));

            hold on;
            
            kpCorrecto(end+1) = Kp;
            tiCorrecto(end+1) = Ti;
            tdCorrecto(end+1) = Td1;
        end
    end
end

% Plotear escalon
plot(t, u)

Tolerancia en tiempo de establecimiento
plot(t, ones(1,length(t))+tolerancia,'--k');
plot(t, ones(1,length(t))-tolerancia,'--k');
plot(tsmax.*[1 1], [-50 50],'k');

Tiempo de subida
plot(trmax.*[1 1], [-50 50],'c');
plot(trmin.*[1 1], [-50 50],'c');

Sobreelongacion maxima
plot(t, ones(1,length(t))+Mpmax,'--g');
plot(t, ones(1,length(t))+Mpmin,'--g');

axis([0 t(end) -0.25 1.2]);
hold off;

figure(2)

s = scatter3(kps(:),Td1s(:),Tis(:), 'red');
hold on;
scatter3(kpCorrecto, tdCorrecto, tiCorrecto, 'green','filled');
alpha(s, 0.2);
axis([exp(Kpmin) exp(Kpmax) exp(Td1min) exp(Td1max) exp(Timin) exp(Timax)]);
xlabel('Kp');ylabel('Td1');zlabel('Ti');title('Combinaciones de Kp, Td1 y Ti admitidos')
set(gca,'xscale','log');set(gca,'yscale','log');set(gca,'zscale','log')
