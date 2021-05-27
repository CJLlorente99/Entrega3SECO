close all; clear all;

%% Inicializar valore
p = 64.986;
K = 2652.28;

experimentosbeta2 = 20;
beta2min = 10;
beta2max = 15;

experimentosbeta = 20;
betamin = 0;
betamax = 5;

experimentoszeta = 20;
zetamin = 2;
zetamax = 10;

experimentosomega = 1;
omegamin = 25;
omegamax = 25;

[betas2, betas, zetas, omegas] = ndgrid(   linspace(beta2min, beta2max, experimentosbeta2),...
                                        linspace(betamin, betamax, experimentosbeta),...
                                        linspace(zetamin, zetamax, experimentoszeta),...
                                        linspace(omegamin, omegamax, experimentosomega));
                                    
n = 1;
    
% Init somethings

x = p.*(0:0.0005/p:0.05);
u = ones(1,length(x));
u(1:ceil(length(u)/6)) = 0;
t = x;

Mpmin = 0.08;
Mpmax = 0.15;

tsmax = 0.5 + t(ceil(length(u)/6)); %segundos
tolerancia = 0.02;

trmax = 0.3 + t(ceil(length(u)/6)); %segundos
trmin = 0.09 + t(ceil(length(u)/6)); %segundos
            
for b2 = 1: experimentosbeta2
    for b = 1: experimentosbeta
        for z = 1: experimentoszeta
            for o = 1: experimentosomega
    
%                 disp(['beta2 ' num2str(b2) ' beta ' num2str(b) ' zeta ' num2str(z) ' omega ' num2str(o)]);

                beta2 = betas2(b2,b,z,o);
                beta = betas(b2,b,z,o);
                zeta = zetas(b2,b,z,o);
                omega = omegas(b2,b,z,o);
                
                disp(['n ' num2str(n) ' de ' num2str(experimentosomega*experimentoszeta*experimentosbeta*experimentosbeta2)]);
%                 disp(['beta2 ' num2str(beta2) ' beta ' num2str(beta) ' zeta ' num2str(zeta) ' omega ' num2str(omega)]);
                n = n + 1;
            %% Calculo de parametros

                Kp = (p^2*(2*beta+1/zeta^2)) / (beta2^2*K);
                Td1 = (beta2*(beta+2)) / (p*(2*beta+1/zeta^2));
                Td2 = -p/(K*Kp);
                Ti = (beta2*zeta^2*(2*beta+1/zeta^2)) / (beta*p);
                Td = Td1+Td2;

                num = [K*Kp*Td1 K*Kp K*Kp/Ti];
                den = [1 p+K*Kp*Td K*Kp K*Kp/Ti];
                sys = tf(num,den);

                %% Condiciones de estabilidad
                if Kp < 0
                    disp(['Kp < 0 ' 'beta2 ' num2str(b2) ' beta ' num2str(b) ' zeta ' num2str(z) ' omega ' num2str(o)]);
                    continue
                    % no estable
                end

                if Ti < 0
                    disp(['Ti < 0 ' 'beta2 ' num2str(b2) ' beta ' num2str(b) ' zeta ' num2str(z) ' omega ' num2str(o)]);
                    continue
                    % no estable
                end

                if Kp*Td <= -p/K
                    disp(['Kp*Td <= -p/K ' 'beta2 ' num2str(b2) ' beta ' num2str(b) ' zeta ' num2str(z) ' omega ' num2str(o)]);
                    continue
                    % no estable
                end

                epsilon = p+K*Kp*Td;

                if epsilon <= 0
                    disp(['epsilon <= 0 ' 'beta2 ' num2str(b2) ' beta ' num2str(b) ' zeta ' num2str(z) ' omega ' num2str(o)]);
                    continue
                    % no estable
                end

                %% Respuesta ante el escalon

                y = lsim(sys, u, t);
                
                %% Condicion de sobreelongacion maxima y minima
                
                if sum(y > (1 + Mpmax)) >= 1
                    disp(['Cond. sobreelongacion max ' 'beta2 ' num2str(b2) ' beta ' num2str(b) ' zeta ' num2str(z) ' omega ' num2str(o)]);
                    continue
                end
                
                if sum(y > (1 + Mpmin)) == 0
                    disp(['Cond. sobreelongacion min ' 'beta2 ' num2str(b2) ' beta ' num2str(b) ' zeta ' num2str(z) ' omega ' num2str(o)]);
                    continue
                end
                    
                %% Condicion de tiempo de establecimiento
                                
                if sum(y(find(t >= tsmax,1):end) > (1 + tolerancia)) >= 1
                    disp(['Cond. Ts max ' 'beta2 ' num2str(b2) ' beta ' num2str(b) ' zeta ' num2str(z) ' omega ' num2str(o)]);
                    continue
                end
                
                if sum(y(find(t >= tsmax,1):end) < (1 - tolerancia)) >= 1
                    disp(['Cond. Ts min ' 'beta2 ' num2str(b2) ' beta ' num2str(b) ' zeta ' num2str(z) ' omega ' num2str(o)]);
                    continue
                end

                %% Condicion de tiempo de subida
                                
                if sum(t(find(y >= 1,1)) >= trmax) >= 1
                    disp(['Cond. tiempo subida max ' 'beta2 ' num2str(b2) ' beta ' num2str(b) ' zeta ' num2str(z) ' omega ' num2str(o)]);
                    continue
                end
                                               
                if sum(t(find(y >= 1,1)) <= trmin) >= 1
                    disp(['Cond. tiempo subida min ' 'beta2 ' num2str(b2) ' beta ' num2str(b) ' zeta ' num2str(z) ' omega ' num2str(o)]);
                    continue
                end
                
                % Plot
                info = ['beta2 ' num2str(beta2) ' beta ' num2str(beta) ' zeta ' num2str(zeta) ' omega ' num2str(omega)];
                h = plot(t,y, 'DisplayName', info);
                h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Trace',repmat({h.DisplayName},size(h.XData)));

                hold on;
            end
        end
    end
end

% Plotear escalon
plot(t, u)

% Tolerancia en tiempo de establecimiento
plot(t, ones(1,length(t))+tolerancia,'--k');
plot(t, ones(1,length(t))-tolerancia,'--k');
plot(tsmax.*[1 1], [-50 50],'k');

% Tiempo de subida
plot(trmax.*[1 1], [-50 50],'c');
plot(trmin.*[1 1], [-50 50],'c');

% Sobreelongacion maxima
plot(t, ones(1,length(t))+Mpmax,'--g');
plot(t, ones(1,length(t))+Mpmin,'--g');

axis([0 t(end) -0.25 1.2]);

