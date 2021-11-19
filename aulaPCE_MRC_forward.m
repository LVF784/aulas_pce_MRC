clear all
close all
clc

%% Conversor Forward
% Par�metros
Vin = 175;
n = 0.343;
C = 2*10^(-6);
L = 2.88*10^(-3);
H = 0.1;
Kpwm = 0.1;
Vo = 25;
R = Vo^2/100;

%% Modelo em tempo cont�nuo
s = tf('s');
G = (Vin*n/(L*C*s^2 + (L/R)*s + 1))
G = H*Kpwm*G;                % incluindo ganho do sensor
% Ordem relativa de G = 2 (polos) - 0 (zeros) = 2

% Settling time (ts) em malha aberta (MA)
stepG = stepinfo(G);
settling_time_G = stepG.SettlingTime;
fprintf('O ts em malha aberta � de %f \n',settling_time_G)

%% Modelo de refer�ncia
% Objetivos: 
% (1) Erro nulo em regime permanente
% (2) Settling time 20% r�pido que em malha aberta
% (3) Sobressinal nulo para a resposta ao degrau

% Note que � necess�rio incluir 2 ou mais polos em Td para que a ordem
% relativa seja igual ou maior � de G

% Garantindo (2)
settling_time_Td = settling_time_G*0.8
tal = settling_time_Td/4;
Td = 1/(tal*s+1);
Td = zpk(Td);
p1 = cell2mat(Td.P);        %polo dominante
p2 = p1*10;                 %polo n�o dominante
Td = Td/(s-p2)

%% Checando (2)
stepTd = stepinfo(Td);
settling_time_Td = stepTd.SettlingTime;
times = settling_time_G/settling_time_Td;
fprintf('Checando condi��o (2) \n')
fprintf('O ts de Td � ~%0.1f vezes mais r�pido que em MA \n',times)

%% Checando (1) - Td(s=0) = 1
Td_0 = evalfr(Td,0);
fprintf('\nChecando condi��o (1) \n')
fprintf('Td(s=0) = %0.1f \n',Td_0)

% N�o deu certo, portanto, conserta-se da seguinte forma
Td = Td/(evalfr(Td,0));
fprintf('\nO Td novo �:')
Td = zpk(Td)

% Checando novamente
Td_0 = evalfr(Td,0);
fprintf('\nChecando condi��o (1) novamente \n')
fprintf('Td(s=0) = %0.1f \n',Td_0)

% A condi��o (3) pode ser checada analiticamente. Sabendo que s� temos
% polos reais em Td, o sobressinal ser� nulo. Podemos
% comprovar olhando para o degrau aplicado a Td
figure(1)
step(Td)
title('Td - resposta ao degrau')

%% Controlador ideal (MRC)
% Utilizando a express�o 8 do m�dulo MRC, obtemos o seguinte controlador:
Cd = minreal(Td/(G*(1-Td)),1e-3)

% Obs.: minreal(.) � uma fun��o que elimina polos e zeros que estejam
% pr�ximos suficiente ao ponto de se cancelarem. O valor de precis�o
% escolhido � de 1e-3.

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
fprintf('O ts em malha aberta � de %f\n',settling_time_G)
fprintf('O ts em malha fechada � de %f\n',settling_time_T)
fprintf('Sendo %.1f vezes mais r�pido para a malha fechada \n',times_T)

% Margem de ganho e margem de fase
figure(3)
margin(T)