clear all
close all
clc

%% Conversor Flyback
% Parâmetros
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

%% Modelo completo em tempo contínuo
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
fprintf('O ts em malha aberta é de %f \n',settling_time_G)

%% Modelo de referência
% Objetivos: 
% (1) Erro nulo em regime permanente
% (2) Settling time 10% mais rápido que em malha aberta
% (3) Sobressinal nulo para a resposta ao degrau

% Note que é necessário incluir 1 ou mais polos em Td para que a ordem
% relativa seja igual à de G

% Também note a presença de um zero de Fase Não Mínima (FNM) em G, sendo
% necessário o incluir no modelo de referência Td

% Garantindo (2)
xp = 10;                        %10% mais rápido
settling_time_Td = settling_time_G*((100-xp)/100);
p1 = -4/(settling_time_Td);

% Zero de FNM
z1 = roots(G.Num{1});

% Polo não dominante
p2 = p1*10;

% Montando a Td
Td = 1*(s-z1)/((s-p1)*(s-p2));
zpk(Td)

%% Checando (2)
stepTd = stepinfo(Td);
settling_time_Td = stepTd.SettlingTime;
times = settling_time_G/settling_time_Td;
percent = 100-100/times;
fprintf('Checando condição (2) \n')
fprintf('O ts de Td é ~%0.1f%% mais rápido que em MA \n',percent)

%% Checando (1) - Td(s=0) = 1
Td_0 = evalfr(Td,0);
fprintf('\nChecando condição (1) \n')
fprintf('Td(s=0) = %0.1f \n',Td_0)

% Não deu certo, portanto, conserta-se da seguinte forma
Td = Td/(evalfr(Td,0));
fprintf('\nO Td novo é:')
Td = zpk(Td)

% A condição (3) pode ser checada analiticamente. Sabendo que só temos
% polos reais em Td, o sobressinal será nulo. Podemos
% comprovar olhando para o degrau aplicado a Td
figure(1)
step(Td)
% Note que a resposta ao degrau inicia na direção oposta. Esta
% característica ocorre pela presença de zeros de fase não mínima.

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
percent_T = 100-100/times_T;
fprintf('O ts em malha aberta é de %f\n',settling_time_G)
fprintf('O ts em malha fechada é de %f\n',settling_time_T)
fprintf('Sendo %.1f%% mais rápido para a malha fechada \n',percent_T)

% Margem de ganho e margem de fase
figure(3)
margin(T)

% Note que diferentemente do caso do conversor Forward e do modelo
% simplificado do conversor Flyback, a margem de ganho é menor (14.8 dB).
% Essa limitação de robustez é causada pela presença de zeros de fase não
% mínima.

%%
%% print secreto
% --------------------------------------------------------------
clc; close all;
 
% -------------------------- PARAM FIG ------------------------
ImageDPI = 1200; salvar = 1;
FontName = 'Times New Roman';
Size = 10.0;
ImageSizeX = 10; %tamanho em cm
ImageSizeY = 0.65*ImageSizeX; %tamanho em cm
FileLabel = ['step_flbk'];

fig1 = figure('Units','centimeters','Position',[4 5 ImageSizeX ImageSizeY],'PaperSize',[ImageSizeX ImageSizeY]);

% -------------------------- PARAM FIG ------------------------
% Dimensões do plot
width = 7.0/ImageSizeX;
height = 0.65*width*ImageSizeX/ImageSizeY;
x = 1.6/ImageSizeX;
y = 1.1/ImageSizeY;

% --------------- EIXO 1 -----------------------
ax1 = axes('position',[x y width height],'FontSize',Size,'FontName',FontName); 
t = 0:1e-4:0.5;
hold on
[plt1,t1] = step(G,t);
plot(t1,plt1,'linewidth',1.3)
hold on
[plt2,t2] = step(T,t);
plot(t2,plt2,'linewidth',1.3)
hold on

% grid on
set(ax1, ... 
  'Box'         , 'on'   , ...
  'TickDir'     , 'in'          , ...
  'XMinorTick'  , 'off'         , ...
  'YMinorTick'  , 'off'         , ...
  'XGrid'       , 'on'          , ...
  'YGrid'       , 'on'          , ... 
  'LineWidth'   , 0.25           , ...
  'FontSize'    , Size        ,...
  'FontName'    , FontName  ,...
  'TickLabelInterpreter', 'Latex');

lg = legend('Malha aberta','Malha fechada');
set(lg,'Interpreter','Latex','Orientation','Vertical',...
    'FontName',FontName,'FontSize',9,'Box', 'on','Location','southeast');

% ------------- Labels -------------------------
xlabel(ax1,'Tempo (s)', 'Interpreter', 'Latex','FontName',FontName,'FontSize',Size);
ylabel(ax1,'Tensão de saída (V)', 'Interpreter', 'Latex','FontName',FontName,'FontSize',Size);

pos = get(gca, 'Position');
pos(1) = 0.12;
pos(2) = 0.17;
pos(3) = 0.82;
pos(4) = 0.8;
set(gca, 'Position', pos)

if salvar == 1
    tic; set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 ImageSizeX ImageSizeY],'PaperSize',[ImageSizeX ImageSizeY]);
    location = '';
    print('-dpdf',strcat([location FileLabel],'.pdf'),strcat('-r',num2str(ImageDPI)))
    fprintf('\nOkk! (%.2f seconds to save)\n',toc);
else
    fprintf('\nOk! Save option disable!!!\n');
end