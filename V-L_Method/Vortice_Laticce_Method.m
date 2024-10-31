%----------------------------------------------------------------------------------------------------------------------------%
clear all;
close all;
clc;

%% Definición del perfil NACA 1216

f_max = 0.01;                   % Curvatura máxima ()
xf = 0.2;                       % Posición curvatura máxima ()
t = 0.16;                       % Espesor ()

%% Definición de parámetros geométricos del ala

lambda = 0.3;                   % Estrechamiento()
Al = 5;                         % Alargamiento ()
kappa = 30;                     % Flecha (grados)
eps = -3;                       % Torsion (grados)
alpha = 5;                      % Angulo de ataque (grados)

% Conversion a radianes de los angulos

kappa = kappa*pi/180;
eps = eps*pi/180;

%% Panelizacion
N_paneles_x = 1;                             % Numero de paneles por cuerda
N_paneles_y = 1000;                           % Numero de paneles por semiala

N_nodos_x = N_paneles_x + 1;                 % Numero de nodos por cuerda
N_nodos_y = N_paneles_y + 1;                 % Numero de nodos por semiala

%----------------------------------------------------------------------------------------------------------------------------%

%% Cálculo de la geometría del ala

semib = sqrt(Al*S)/2;                                                       %% Envergadura (m)
C_gm = semib/Al;                                                            %% Cuerda geométrica media (m)
C_r = 2*C_gm/(1 + lambda);                                                  %% Cuerda en la raíz (m)
C_t = C_r*lambda;                                                           %% Cuerda en la punta (m)
C_am = 2*C_r*(1 + lambda + lambda^2)/(3*(1 + lambda));                      %% Cuerda aerodinámica media (m)


% Coordenadas de los puntos BA, BS, 1/4 y 3/4 de C_r

x_BA_cr = 0;
x_1_4_cr = C_r/4;
x_3_4_cr = C_r*3/4;
x_BS_cr = C_r;

% Coordenadas de los puntos BA, BS, 1/4 y 3/4 de C_t

x_1_4_ct = x_1_4_cr + tan(kappa)*semib/2;
x_BA_ct = x_1_4_ct - C_t/4;
x_3_4_ct = x_BA_ct + C_t*3/4;
x_BS_ct = x_BA_ct + C_t;

% Flecha para BA, 1/4, 3/4 y BS

kappa_BA = double(atan(x_BA_ct/(semib/2)));
kappa_1_4 = kappa;
kappa_3_4 = atan((x_3_4_ct - x_3_4_cr)/(semib/2));
kappa_BS = atan((x_BS_ct - x_BS_cr)/(semib/2));

%% Perfil Aerodinamico NACA 

% Derivada para Nx = 1 en 3/4
    if xf >= 3/4
        dzdx = (2*xf - (3/4)*2)*f_max/(xf^2);
    elseif xf < 3/4        
        dzdx = (2*xf - (3/4)*2)*f_max/((1 - xf)^2);
    end
    
%% Panelizacion del ala
% Calculo Coordenadas BA y BS por panel

y_d = zeros(N_nodos_y, N_paneles_x);
y_i = zeros(N_nodos_y, N_paneles_x);

x_BA_d = zeros(N_nodos_y, N_paneles_x);
x_BA_i = zeros(N_nodos_y, N_paneles_x);

x_BS_d = zeros(N_nodos_y, N_paneles_x);
x_BS_i = zeros(N_nodos_y, N_paneles_x);

x_1_4_d = zeros(N_nodos_y, N_paneles_x);
x_1_4_i = zeros(N_nodos_y, N_paneles_x);

x_3_4_d = zeros(N_nodos_y, N_paneles_x);
x_3_4_i = zeros(N_nodos_y, N_paneles_x);

y_d = linspace(0,semib/2, N_nodos_y);

for i = 1: N_nodos_y
    
    y_i(i) = -y_d(i);

    x_BA_d(i) = y_d(i)*tan(kappa_BA);
    x_BS_d(i) = C_r + y_d(i)*tan(kappa_BS);

    x_BA_i(i) = x_BA_d(i);
    x_BS_i(i) = x_BS_d(i);

    x_1_4_d(i) = x_1_4_cr + y_d(i)*tan(kappa_1_4);
    x_3_4_d(i) = x_3_4_cr + y_d(i)*tan(kappa_3_4);

    x_1_4_i(i) = x_1_4_d(i);
    x_3_4_i(i) = x_3_4_d(i);

end

% Calculo Coordenadas Puntos de Control por Panel
    x_pc_d = zeros(N_paneles_y, N_paneles_x);
    y_pc_d = zeros(N_paneles_y, N_paneles_x);

    x_pc_i = zeros(N_paneles_y, N_paneles_x);
    y_pc_i = zeros(N_paneles_y, N_paneles_x);
    
for i = 1:N_paneles_y

    x_pc_d(i) =  (x_3_4_d(i) + x_3_4_d(i + 1))/2;
    y_pc_d(i) =  (y_d(i) + y_d(i + 1))/2;

    x_pc_i(i) = x_pc_d(i);
    y_pc_i(i) = -y_pc_d(i);
end

% Calculo Coordenadas Esquinas torbellinos
    x_tor1_d = zeros(N_paneles_y, N_paneles_x);
    y_tor1_d = zeros(N_paneles_y, N_paneles_x);

    x_tor2_d = zeros(N_paneles_y, N_paneles_x);
    y_tor2_d = zeros(N_paneles_y, N_paneles_x);

    x_tor1_i = zeros(N_paneles_y, N_paneles_x);
    y_tor1_i = zeros(N_paneles_y, N_paneles_x);

    x_tor2_i = zeros(N_paneles_y, N_paneles_x);
    y_tor2_i = zeros(N_paneles_y, N_paneles_x);
    
for i = 1:N_paneles_y

    x_tor1_d(i) =  x_1_4_d(i);
    y_tor1_d(i) =  y_d(i);

    x_tor2_d(i) =  x_1_4_d(i + 1);
    y_tor2_d(i) =  y_d(i + 1);

    x_tor1_i(i) =  x_tor2_d(i);
    y_tor1_i(i) =  - y_tor2_d(i);

    x_tor2_i(i) =  x_tor1_d(i);
    y_tor2_i(i) =  - y_tor1_d(i);
end

% Envergadura y cuerda media geometrica de cada panel

cij_d = zeros(N_paneles_y, N_paneles_x);
bij_d = zeros(N_paneles_y, N_paneles_x);

for i = 1:N_paneles_y

    bij_d(i) = y_d(i + 1) - y_d(i);
    cij_d(i) = C_r + (C_r - C_t)*y_pc_d(i)/semib;

end


%% Calculo velocidad inducida por torbellinos y velocidad aguas arriba

V_d = zeros(N_paneles_y, N_paneles_y);
V_i = zeros(N_paneles_y, N_paneles_y);

B = zeros(N_paneles_y, N_paneles_x);

for i = 1:N_paneles_y
    for j = 1:N_paneles_y
        a = x_pc_d(i) - x_tor1_d(j);
        b = y_pc_d(i) - y_tor1_d(j);
        c = x_pc_d(i) - x_tor2_d(j);
        d = y_pc_d(i) - y_tor2_d(j);
        e = sqrt(a^2 + b^2);
        f = sqrt(c^2 + d^2);
        g = x_tor2_d(j) - x_tor1_d(j);
        h = y_tor2_d(j) - y_tor1_d(j);

        k = ((g*a + h*b)/e) - ((g*c + h*d)/f);
        l = (-(1 + a/e)/b) + (1 + c/f)/d;


        V_d(i,j) = 1/(4*pi)*(k/(a*d - c*b) + l);
    end

        for j = 1:N_paneles_y
        a = x_pc_i(i) - x_tor1_i(j);
        b = y_pc_i(i) - y_tor1_i(j);
        c = x_pc_i(i) - x_tor2_i(j);
        d = y_pc_i(i) - y_tor2_i(j);
        e = sqrt(a^2 + b^2);
        f = sqrt(c^2 + d^2);
        g = x_tor2_i(j) - x_tor1_i(j);
        h = y_tor2_i(j) - y_tor1_i(j);

        k = ((g*a + h*b)/e) - ((g*c + h*d)/f);
        l = (-(1 + a/e)/b) + (1 + c/f)/d;


        V_i(i,j) = 1/(4*pi)*(k/(a*d - c*b) + l);
        end

        B(i, 1) = alpha + eps*y_pc_d(i)/semib - dzdx;

end

%% Resolucion sistema
A = V_d + V_i;

Gamma_d = -(A^(-1))*B;


%% Calculo variables aerodinamicas semiala derecha

Cl_d = zeros(N_paneles_y, N_paneles_x);
Cm_oy_d = zeros(N_paneles_y, N_paneles_x);
Cd_d = zeros(N_paneles_y, N_paneles_x);

CL = 0;
CM_oy = 0;
CD = 0;
for i = 1: N_paneles_y
    alphai = 0;
    Cl_d(i) = 2*Gamma_d(i)/cij_d(i);

    CL = CL + 2*Cl_d(i)*bij_d(i)*cij_d(i);
    
    Cm_oy_d(i) = -Cl_d(i)*x_1_4_d(i)/cij_d(i);

    CM_oy = CM_oy + 2*Cm_oy_d(i)*cij_d(i)*bij_d(i);

    for j = 1:N_paneles_y
        
        alphai = alphai + Gamma_d(j)*A(i,j);
    
    end
    Cd_d(i) = -Cl_d(i)*alphai;

    CD = CD + 2*Cd_d(i)*cij_d(i)*bij_d(i);

    

end

%% Acople ala derecha e izquierda
N_pc_y = 2*N_paneles_y - 1;

y_pc = zeros(N_pc_y, N_paneles_x);
Cl  = zeros(N_pc_y, N_paneles_x);
Cm_oy = zeros(N_pc_y, N_paneles_x);
Cd = zeros(N_pc_y, N_paneles_x);

for i = 1:N_paneles_y

    y_pc(i) = y_pc_i(N_paneles_y - i + 1);
    y_pc(N_paneles_y + i -1) = y_pc_d(i);

    Cl(i) = Cl_d(N_paneles_y - i + 1);
    Cl(N_paneles_y + i -1) = Cl_d(i);

    Cm_oy(i) = Cm_oy_d(N_paneles_y - i + 1);
    Cm_oy(N_paneles_y + i -1) = Cm_oy_d(i);

    Cd(i) = Cd_d(N_paneles_y - i + 1);
    Cd(N_paneles_y + i -1) = Cd_d(i);

end


%% Figuras
% Representacion Semiala Derecha

figure(1)
hold
axis equal
% Borde Ataque y Borde Salida
plot(x_BA_d, y_d, 'Color', '#0072BD')
plot(x_BA_i, y_i, 'Color', '#0072BD')
plot(x_BS_d, y_d, 'Color', '#0072BD')
plot(x_BS_i, y_i, 'Color', '#0072BD')

% Linea 1/4
plot(x_1_4_d, y_d, 'blue')
plot(x_1_4_i, y_i, 'blue')

% Linea 3/4
plot(x_3_4_d, y_d, 'blue')
plot(x_3_4_i, y_i, 'blue')



% Nodos BA y BS
scatter(x_BA_d, y_d, 'red')
scatter(x_BS_d, y_d, 'red')

scatter(x_BA_i, y_i, 'red')
scatter(x_BS_i, y_i, 'red')


% Puntos de control
scatter(x_pc_d, y_pc_d, 's', 'k')
scatter(x_pc_i, y_pc_i, 's', 'k')


% Esquinas torbellinos
scatter(x_tor1_d, y_tor1_d, 'blue')
scatter(x_tor2_d, y_tor2_d,'s', 'red')

scatter(x_tor1_i, y_tor1_i, 'green')
scatter(x_tor2_i, y_tor2_i,'s', 'yellow')

%legend('Ba','BS','1/4','3/4','Nodos','Puntos de Control')



% Coeficiente Sustentacion
figure(2)

plot(y_pc, Cl);

% Coeficiente Momentos
figure(3)

plot(y_pc, Cm_oy);

% Coeficiente Resistencia
figure(4)

plot(y_pc, Cd);
