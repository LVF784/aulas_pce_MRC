clear all
close all
clc

%% Conversor Forward
% Parâmetros
Vin = 175;
n = 0.343;
C = 2*10^(-6);
L = 2.88*10^(-3);
H = 0.1;
Kpwm = 0.1;
Vo = 25;
R = Vo^2/100;

%% Modelo em tempo contínuo
s = tf('s');
G = (Vin*n/(L*C*s^2 + (L/R)*s + 1))
G = H*Kpwm*G;                % incluindo ganho do sensor
% Ordem relativa de G = 2 (polos) - 0 (zeros) = 2

% Settling time (ts) em malha aberta (MA)
stepG = stepinfo(G);
settling_time_G = stepG.SettlingTime;
fprintf('O ts em malha aberta é de %f \n',settling_time_G)

%% Modelo de referência
% Objetivos: 
% (1) Erro nulo em regime permanente
% (2) Settling time 20% rápido que em malha aberta
% (3) Sobressinal nulo para a resposta ao degrau

% Note que é necessário incluir 2 ou mais polos em Td para que a ordem
% relativa seja igual ou maior à de G

% Garantindo (2)
settling_time_Td = settling_time_G*0.8
tal = settling_time_Td/4;
Td = 1/(tal*s+1);
Td = zpk(Td);
p1 = cell2mat(Td.P);        %polo dominante
p2 = p1*10;                 %polo não dominante
Td = Td/(s-p2)

%% Checando (2)
stepTd = stepinfo(Td);
settling_time_Td = stepTd.SettlingTime;
times = settling_time_G/settling_time_Td;
fprintf('Checando condição (2) \n')
fprintf('O ts de Td é ~%0.1f vezes mais rápido que em MA \n',times)

%% Checando (1) - Td(s=0) = 1
Td_0 = evalfr(Td,0);
fprintf('\nChecando condição (1) \n')
fprintf('Td(s=0) = %0.1f \n',Td_0)

% Não deu certo, portanto, conserta-se da seguinte forma
Td = Td/(evalfr(Td,0));
fprintf('\nO Td novo é:')
Td = zpk(Td)

% Checando novamente
Td_0 = evalfr(Td,0);
fprintf('\nChecando condição (1) novamente \n')
fprintf('Td(s=0) = %0.1f \n',Td_0)

% A condição (3) pode ser checada analiticamente. Sabendo que só temos
% polos reais em Td, o sobressinal será nulo. Podemos
% comprovar olhando para o degrau aplicado a Td
figure(1)
step(Td)
title('Td - resposta ao degrau')

%% Controlador ideal (MRC)
% Utilizando a expressão 8 do módulo MRC, obtemos o seguinte controlador:
Cd = minreal(Td/(G*(1-Td)),1e-3)

% Obs.: minreal(.) é uma função que elimina polos e zeros que estejam
% próximos suficiente ao ponto de se cancelarem. O valor de precisão
% escolhido é de 1e-3.

%% Fechando a malha com o controlador projetado
T = Cd*G/(1+Cd*G);

% Comparando graficamente com a malha aberta
figure(2)
step(G)
hold on
step(T)
title('Resposta ao degrau em malha aberta e malha fechada')
legend('malha aberta', 'malha fechada')

% Settling time
stepTd = stepinfo(T);
settling_time_T = stepTd.SettlingTime;
times_T = settling_time_G/settling_time_T;
fprintf('O ts em malha aberta é de %f\n',settling_time_G)
fprintf('O ts em malha fechada é de %f\n',settling_time_T)
fprintf('Sendo %.1f vezes mais rápido para a malha fechada \n',times_T)

% Margem de ganho e margem de fase
figure(3)
margin(T)