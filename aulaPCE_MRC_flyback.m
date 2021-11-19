clear all
close all
clc

%% Conversor Flyback
% Par�metros
Vin = 100;
Nt = 1/2;
C = 47*10^(-6);
Lm = 96*10^(-3);
H = 0.1;
Kpwm = 0.1;
Vo = 30;
Po = 100;
Ro = Vo^2/Po;
D = Nt*Vo/(Vin + Nt*Vo);

%% Modelo completo em tempo cont�nuo
s = tf('s');
Gd0 = Vo/(D*(1-D));
w0 = Nt*(1-D)/sqrt(Lm*C);
wz = Nt^2*(1-D)^2*Ro/(D*Lm);
Q0 = Nt*(1-D)*Ro*sqrt(C/Lm);
G = Gd0*(1-s/wz)/(s^2/w0^2 + s/(Q0*w0) + 1);
zpk(G)
% Ordem relativa de G = 2 (polos) - 1 (zeros) = 1

G = Kpwm*H*G;          %incluindo ganho do sensor e pwm

% Settling time (ts) em malha aberta (MA)
stepG = stepinfo(G);
settling_time_G = stepG.SettlingTime;
fprintf('O ts em malha aberta � de %f \n',settling_time_G)

%% Modelo de refer�ncia
% Objetivos: 
% (1) Erro nulo em regime permanente
% (2) Settling time 10% mais r�pido que em malha aberta
% (3) Sobressinal nulo para a resposta ao degrau

% Note que � necess�rio incluir 1 ou mais polos em Td para que a ordem
% relativa seja igual � de G

% Tamb�m note a presen�a de um zero de Fase N�o M�nima (FNM) em G, sendo
% necess�rio o incluir no modelo de refer�ncia Td

% Garantindo (2)
xp = 10;                        %10% mais r�pido
settling_time_Td = settling_time_G*((100-xp)/100);
p1 = -4/(settling_time_Td);

% Zero de FNM
z1 = roots(G.Num{1});

% Polo n�o dominante
p2 = p1*10;

% Montando a Td
Td = 1*(s-z1)/((s-p1)*(s-p2));
zpk(Td)

%% Checando (2)
stepTd = stepinfo(Td);
settling_time_Td = stepTd.SettlingTime;
times = settling_time_G/settling_time_Td;
percent = 100-100/times;
fprintf('Checando condi��o (2) \n')
fprintf('O ts de Td � ~%0.1f%% mais r�pido que em MA \n',percent)

%% Checando (1) - Td(s=0) = 1
Td_0 = evalfr(Td,0);
fprintf('\nChecando condi��o (1) \n')
fprintf('Td(s=0) = %0.1f \n',Td_0)

% N�o deu certo, portanto, conserta-se da seguinte forma
Td = Td/(evalfr(Td,0));
fprintf('\nO Td novo �:')
Td = zpk(Td)

% A condi��o (3) pode ser checada analiticamente. Sabendo que s� temos
% polos reais em Td, o sobressinal ser� nulo. Podemos
% comprovar olhando para o degrau aplicado a Td
figure(1)
step(Td)
% Note que a resposta ao degrau inicia na dire��o oposta. Esta
% caracter�stica ocorre pela presen�a de zeros de fase n�o m�nima.

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
percent_T = 100-100/times_T;
fprintf('O ts em malha aberta � de %f\n',settling_time_G)
fprintf('O ts em malha fechada � de %f\n',settling_time_T)
fprintf('Sendo %.1f%% mais r�pido para a malha fechada \n',percent_T)

% Margem de ganho e margem de fase
figure(3)
margin(T)

% Note que diferentemente do caso do conversor Forward e do modelo
% simplificado do conversor Flyback, a margem de ganho � menor (14.8 dB).
% Essa limita��o de robustez � causada pela presen�a de zeros de fase n�o
% m�nima.