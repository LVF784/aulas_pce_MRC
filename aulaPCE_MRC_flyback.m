clear all
close all
clc

%% Conversor Flyback
% Parâmetros
Vin = 400;
n = 1/12;
C = 47*10^(-6);
Rse = 10e-3;
L = 96.02*10^(-3);
H = 0.1;
Vo = 50;
R = Vo^2/300;
fs = 25e3;

%% Modelo em tempo contínuo
s = tf('s');
G = (Vin/sqrt((2*L*fs)/R))*(1+s*Rse*C)/(1+s*R*C)
% Ordem relativa de G = 1 (polos) - 0 (zeros) = 1

% Settling time (ts) em malha aberta (MA)
stepG = stepinfo(G);
settling_time_G = stepG.SettlingTime;
fprintf('O ts em malha aberta é de %f \n',settling_time_G)

%% Modelo de referência
% Objetivos: 
% (1) Erro nulo em regime permanente
% (2) Settling time 5x mais rápido que em malha aberta
% (3) Sobressinal nulo

% Note que é necessário incluir 1 ou mais polos em Td para que a ordem
% relativa seja igual à de G

% Garantindo (2)
settling_time_Td = settling_time_G/5;
tal = settling_time_Td/4;
Td = 1/(tal*s+1);
Td = zpk(Td)

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

% A condição (3) pode ser checada analiticamente. Sabendo que só temos
% polos reais em Td sem nenhum zero, o sobressinal será nulo. Podemos
% comprovar olhando para o degrau aplicado a Td
figure(1)
step(Td)

%% Controlador ideal (MRC)
% Utilizando a expressão 8 do módulo MRC, obtemos o seguinte controlador:
Cd = minreal(Td/(G*(1-Td*H)),1e-3)

% Obs.: minreal(.) é uma função que elimina polos e zeros que estejam
% próximos suficiente ao ponto de se cancelarem. O valor de precisão
% escolhido é de 1e-3.

%% Fechando a malha com o controlador projetado
T = Cd*G/(1+Cd*G*H);

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