%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Tecnologico de Costa Rica
% Escuela de Ingenieria Electronica
% EL-5409 Laboratorio de Control Automatico
%
% Estudiante: Jose David Luna Herrera
% Carne: 2020114728
%
% Profesor: Ing. Luis C. Rosales
% II Semestre 2026
% Proyecto individual 3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ==============================================================
% PROYECTO CORTO #3 - CONTROL AUTOMATICO
%
% El programa:
% - ingresa polos, ceros y ganancia de la planta
% - calcula el lugar de las raices
% - escoge los polos deseados con el mouse
% - obtene un compensador dinamico
% - comprueba los polos y la respuesta al escalon
%
% Para los polos complejos se agrega automaticamente el conjugado
% para trabajar con polinomios de coeficientes reales.
% ==============================================================


clc;
clear;
close all;

fprintf('\n====================================================\n');
fprintf('       PROYECTO CORTO #3 - CONTROL AUTOMATICO\n');
fprintf('       VERSION FINAL - SIN TOOLBOX\n');
fprintf('====================================================\n\n');

%% ==============================================================
% 1. Datos que voy a ingresar

fprintf('Ingrese los CEROS de G(s).\n');
fprintf('Ejemplo: [0 -2] | Si no hay ceros: []\n');
z = input('Ceros = ');

fprintf('\nIngrese los POLOS de G(s).\n');
fprintf('Ejemplo: [-1 -3 -5]\n');
p = input('Polos = ');

fprintf('\nIngrese la ganancia de la planta K [Enter = 1]: ');
Kplanta = input('');

if isempty(Kplanta)
    Kplanta = 1;
end

z = z(:).';
p = p(:).';

if isempty(p)
    error('Debe ingresar por lo menos un polo.');
end

if ~isscalar(Kplanta) || ~isreal(Kplanta)
    error('La ganancia K de la planta debe ser un numero real.');
end

% La funcion de transferencia debe tener coeficientes reales.
if ~sonConjugados(z) || ~sonConjugados(p)
    error(['Los polos y ceros complejos deben aparecer en pares ', ...
           'conjugados para obtener una funcion con coeficientes reales.']);
end

%% ==============================================================
% 2. Construyo la funcion de transferencia
% G(s) = K * N(s) / D(s)

denG = real(poly(p));

if isempty(z)
    numG = Kplanta;
else
    numG = Kplanta * real(poly(z));
end

denG = quitarZeros(denG);
numG = quitarZeros(numG);

% Dejo el denominador con coeficiente principal igual a 1.
numG = numG / denG(1);
denG = denG / denG(1);

% La planta debe ser propia para el procedimiento clasico.
if length(numG) >= length(denG)
    error(['La planta ingresada no es estrictamente propia. ', ...
           'Se requiere grado(N) < grado(D).']);
end

n = length(denG)-1;

fprintf('\n----------------------------------------------------\n');
fprintf('FUNCION DE TRANSFERENCIA DE LA PLANTA\n');
fprintf('----------------------------------------------------\n');
fprintf('G(s) = ');
imprimirFraccion(numG,denG);
fprintf('\nNumerador N(s): ');
imprimirPolinomio(numG);
fprintf('\nDenominador D(s): ');
imprimirPolinomio(denG);
fprintf('\nOrden de la planta: %d\n',n);

%% ==============================================================
% 3. ECUACION CARACTERISTICA INICIAL
%
% 1 + G(s) = 0
%
% Como Kplanta ya esta incluido en N(s):
%
% D(s) + N(s) = 0
% ==============================================================

[dPad,nPad] = igualarLongitud(denG,numG);
charInicial = quitarZeros(dPad+nPad);

fprintf('\n----------------------------------------------------\n');
fprintf('ECUACION CARACTERISTICA INICIAL\n');
fprintf('----------------------------------------------------\n');
fprintf('D(s) + N(s) = ');
imprimirPolinomio(charInicial);
fprintf(' = 0\n');

%% ==============================================================
% 4. Calculo el Root Locus sin usar rlocus()
% Voy cambiando K y calculando las raices de:
% D(s) + K*N(s) = 0

fprintf('\nCalculando Root Locus...\n');

% Pruebo un rango de ganancias para dibujar el lugar de las raices.
kVec = [0 logspace(-6,4,3000)];

R = nan(n,length(kVec));

for ii = 1:length(kVec)

    kk = kVec(ii);

    [dTemp,nTemp] = igualarLongitud(denG,numG);
    polTemp = quitarZeros(dTemp + kk*nTemp);

    rr = roots(polTemp);

    if ii == 1
        rr = ordenarRaices(rr);
    else
        rr = asociarRaices(rr,R(:,ii-1));
    end

    R(:,ii) = rr(:);
end

%% ==============================================================
% 5. Grafica del lugar de las raices

figure('Name','Proyecto Corto #3 - Root Locus', ...
       'NumberTitle','off');

hold on;
grid on;
box on;

for ii = 1:n
    plot(real(R(ii,:)),imag(R(ii,:)), ...
        'b','LineWidth',1.4);
end

plot(real(p),imag(p), ...
    'rx','MarkerSize',11,'LineWidth',2);

if ~isempty(z)
    plot(real(z),imag(z), ...
        'go','MarkerSize',9,'LineWidth',2);
end

xline(0,'k:','LineWidth',1);
yline(0,'k:','LineWidth',1);

xlabel('Parte real \sigma');
ylabel('Parte imaginaria j\omega');
title('Lugar de las raices - Seleccione los polos deseados');

if isempty(z)
    legend('Root Locus','Polos','Location','best');
else
    legend('Root Locus','Polos','Ceros','Location','best');
end

% Ajusto los limites de la grafica para que se pueda seleccionar mejor.
pts = R(isfinite(R));

if ~isempty(pts)

    xr = real(pts);
    yr = imag(pts);

    xmin = min([xr(:);real(p(:));real(z(:))]);
    xmax = max([xr(:);real(p(:));real(z(:))]);
    ymin = min([yr(:);imag(p(:));imag(z(:))]);
    ymax = max([yr(:);imag(p(:));imag(z(:))]);

    dx = max(0.5,0.12*(xmax-xmin+1));
    dy = max(0.5,0.12*(ymax-ymin+1));

    xlim([xmin-dx,xmax+dx]);
    ylim([ymin-dy,ymax+dy]);
end

%% ==============================================================
% 6. Seleccion de los polos deseados
% Los puedo escoger sobre o fuera del Root Locus.
% Si es complejo, agrego tambien el conjugado.

fprintf('\n----------------------------------------------------\n');
fprintf('SELECCION DE POLOS DESEADOS\n');
fprintf('----------------------------------------------------\n');

fprintf(['Ahora puede seleccionar la NUEVA ubicacion deseada de los ', ...
         'polos.\n']);
fprintf(['Los clicks pueden hacerse FUERA del Root Locus: esa es ', ...
         'precisamente la ubicacion que el compensador debe lograr.\n\n']);

fprintf(['La planta es de orden %d. Seleccione %d polos principales.\n'], ...
        n,n);
fprintf(['El programa agregara automaticamente %d polos auxiliares ', ...
         'para sintetizar el compensador dinamico.\n'],max(0,n-1));

fprintf(['Para un polo complejo, haga UN CLICK y el conjugado se ', ...
         'agrega automaticamente.\n']);
fprintf(['Si la planta es de orden impar, el ultimo polo se toma real ', ...
         'automaticamente para conservar coeficientes reales.\n']);
fprintf(['Los polos deseados pueden estar FUERA del Root Locus original.\n\n']);

fprintf('Presione ENTER despues de completar la seleccion.\n\n');

polosDeseados = [];

while length(polosDeseados) < n

    faltan = n-length(polosDeseados);

    fprintf('Faltan %d polo(s). Haga CLICK sobre la ubicacion deseada.\n',faltan);

    [xc,yc] = ginput(1);

    candidato = xc + 1i*yc;

    % Si el click no entrega un numero valido, vuelvo a pedirlo.
    if ~isfinite(candidato)
        fprintf('Click no valido. Intente nuevamente.\n');
        continue;
    end

    % ----------------------------------------------------------
    % Regla para conservar coeficientes reales:
    %
    % - Si quedan 2 o mas posiciones, un click complejo agrega
    %   automaticamente su conjugado.
    % - Si queda EXACTAMENTE 1 posicion, ese ultimo polo debe ser
    %   real. Por eso se acepta el click como polo real aunque el
    %   mouse tenga una pequena desviacion vertical del eje real.
    % ----------------------------------------------------------

    if faltan == 1

        % Ultimo polo: para una ecuacion real de grado impar debe
        % ser real. Se ignora una pequena componente imaginaria
        % accidental del click.
        candidato = real(candidato);

        polosDeseados(end+1) = candidato;

        plot(candidato,0, ...
            'ms','MarkerSize',11,'LineWidth',2);

        fprintf('Agregado ultimo polo real: %.6g\n',candidato);

    elseif abs(imag(candidato)) > 1e-7

        % Par complejo conjugado.
        polosDeseados(end+1) = candidato;
        polosDeseados(end+1) = conj(candidato);

        plot(real(candidato),imag(candidato), ...
            'ms','MarkerSize',11,'LineWidth',2);

        plot(real(conj(candidato)),imag(conj(candidato)), ...
            'ms','MarkerSize',11,'LineWidth',2);

        fprintf('Agregado el par: %.6g %+.6gj y %.6g %+.6gj\n', ...
            real(candidato),imag(candidato), ...
            real(conj(candidato)),imag(conj(candidato)));

    else

        % Polo real.
        polosDeseados(end+1) = real(candidato);

        plot(real(candidato),0, ...
            'ms','MarkerSize',11,'LineWidth',2);

        fprintf('Agregado polo real: %.6g\n',real(candidato));
    end
end

% Ordenar polos para una salida limpia.
polosDeseados = ordenarRaices(polosDeseados).';

fprintf('\nPolos deseados seleccionados:\n');
disp(polosDeseados);

%% ==============================================================
% 7. Armo la ecuacion caracteristica deseada
% Los polos auxiliares se colocan mas a la izquierda para que
% no sean los que dominen la respuesta.

polosAux = [];

if n > 1

    parteRealDeseada = real(polosDeseados);
    margen = max(1,max(abs(parteRealDeseada)));

    % Los polos auxiliares quedan 5 veces mas a la izquierda que
    % la magnitud real maxima de los polos seleccionados.
    baseAux = 5*margen;

    for kk=1:n-1
        polosAux(end+1) = -(baseAux + 2*(kk-1)*margen);
    end
end

polosTotalesDeseados = [polosDeseados(:).' polosAux(:).'];

Pd = real(poly(polosTotalesDeseados));
Pd = quitarZeros(Pd);
Pd = Pd/Pd(1);

fprintf('\n----------------------------------------------------\n');
fprintf('NUEVA ECUACION CARACTERISTICA DESEADA\n');
fprintf('----------------------------------------------------\n');

fprintf('Polos principales seleccionados por el usuario:\n');
disp(polosDeseados);

if ~isempty(polosAux)
    fprintf('Polos auxiliares agregados para el compensador dinamico:\n');
    disp(polosAux);
end

fprintf('Polos totales usados para la sintesis:\n');
disp(polosTotalesDeseados.');

fprintf('P_d(s) = ');
imprimirPolinomio(Pd);
fprintf(' = 0\n');

%% ==============================================================
% 8. Calculo del compensador
% Uso la ecuacion diofantina:
% D(s)A(s) + N(s)B(s) = Pd(s)
%
% Con esto puedo colocar los polos deseados aunque no esten
% sobre el Root Locus original.

fprintf('\nCalculando compensador dinamico...\n');

[A,B,errDioph] = calcularCompensador(DnormSeguro(denG), ...
                                      NnormSeguro(numG),Pd);

if errDioph > 1e-6
    warning(['La solucion de la ecuacion diofantina tiene un error ', ...
             'numerico de %.3e.'],errDioph);
end

A = quitarZeros(real(A));
B = quitarZeros(real(B));

% Normalizo A para que quede monico.
if abs(A(1)) < 1e-12
    error('No fue posible obtener un denominador valido para el compensador.');
end

escala = A(1);
A = A/escala;
B = B/escala;

%% ==============================================================
% 9. Muestro el compensador obtenido

fprintf('\n----------------------------------------------------\n');
fprintf('COMPENSADOR\n');
fprintf('----------------------------------------------------\n');

fprintf('Numerador B(s): ');
imprimirPolinomio(B);

fprintf('\nDenominador A(s): ');
imprimirPolinomio(A);

fprintf('\n\nC(s) = ');
imprimirFraccion(B,A);

%% ==============================================================
% 10. Ecuacion caracteristica del sistema compensado

[Dca,Nca] = igualarLongitud(denG,numG);
[DA,~] = igualarLongitud(conv(denG,A),conv(numG,B));

% DA contiene D*A; construimos por separado para imprimir/verificar.
DA = conv(denG,A);
NB = conv(numG,B);

[DAp,NBp] = igualarLongitud(DA,NB);
charComp = quitarZeros(DAp+NBp);
charComp = charComp/charComp(1);

fprintf('\n----------------------------------------------------\n');
fprintf('ECUACION CARACTERISTICA DEL SISTEMA COMPENSADO\n');
fprintf('----------------------------------------------------\n');
fprintf('D(s)A(s) + N(s)B(s) = ');
imprimirPolinomio(charComp);
fprintf(' = 0\n');

%% ==============================================================
% 11. Calculo y verifico los polos finales

polosComp = roots(charComp);
polosComp = ordenarRaices(polosComp).';

fprintf('\nPolos del sistema compensado:\n');
disp(polosComp);

% Comparo contra todos los polos usados en Pd,
% incluyendo los auxiliares.
polosCompOrdenados = sortComplejo(polosComp);
polosTotalesOrdenados = sortComplejo(polosTotalesDeseados);

errorPolosTotal = max(abs(polosCompOrdenados-polosTotalesOrdenados));

% Tambien reviso por separado los polos principales.
polosPrincipalesOrdenados = sortComplejo(polosDeseados);
erroresPrincipales = zeros(length(polosPrincipalesOrdenados),1);

for ii=1:length(polosPrincipalesOrdenados)
    erroresPrincipales(ii) = min(abs(polosComp-polosPrincipalesOrdenados(ii)));
end

errorPolos = max(erroresPrincipales);

fprintf('Error maximo de TODOS los polos (incluyendo auxiliares): %.3e\n', ...
        errorPolosTotal);
fprintf('Error maximo de los polos principales seleccionados: %.3e\n', ...
        errorPolos);

%% ==============================================================
% 12. Funcion de transferencia de lazo cerrado
% T(s) = N(s)B(s) / [D(s)A(s) + N(s)B(s)]

numT = NB;
denT = charComp;

numT = numT/denT(1);
denT = denT/denT(1);

fprintf('\n----------------------------------------------------\n');
fprintf('FUNCION DE TRANSFERENCIA COMPENSADA\n');
fprintf('----------------------------------------------------\n');

fprintf('T(s) = ');
imprimirFraccion(numT,denT);

%% ==============================================================
% 13. Polos y ceros del sistema final

cerosComp = roots(numT);
cerosComp = quitarRaicesInvalidas(cerosComp);

figure('Name','Polos y ceros - Sistema compensado', ...
       'NumberTitle','off');

hold on;
grid on;
box on;

if ~isempty(polosComp)
    plot(real(polosComp),imag(polosComp), ...
        'rx','MarkerSize',11,'LineWidth',2);
end

if ~isempty(cerosComp)
    plot(real(cerosComp),imag(cerosComp), ...
        'go','MarkerSize',9,'LineWidth',2);
end

xline(0,'k:');
yline(0,'k:');

xlabel('Parte real');
ylabel('Parte imaginaria');
title('Polos y ceros del sistema compensado');

if isempty(cerosComp)
    legend('Polos','Location','best');
else
    legend('Polos','Ceros','Location','best');
end

%% ==============================================================
% 14. Respuestas al escalon
% Comparo la planta con el sistema compensado y,
% aparte, muestro la respuesta del compensador.

fprintf('\nCalculando respuestas al escalon...\n');
fprintf(['Se generaran dos graficas: planta vs sistema compensado, ', ...
         'y la respuesta del compensador como complemento.\n']);

Tfinal = 10;
Npts = 5000;

[tG,yG] = respuestaEscalon(numG,denG,Tfinal,Npts);
[tC,yC] = respuestaEscalon(B,A,Tfinal,Npts);
[tT,yT] = respuestaEscalon(numT,denT,Tfinal,Npts);

% Grafica de la planta y del sistema compensado.
figure('Name','Respuesta Planta vs Sistema Compensado', ...
       'NumberTitle','off');

hold on;
grid on;
box on;

plot(tG,yG,'LineWidth',1.5);
plot(tT,yT,'LineWidth',1.5);

legend('Planta G(s)','Sistema compensado T(s)', ...
       'Location','best');

xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Comparacion de respuesta al escalon');

% Grafica aparte para el compensador.
if all(isfinite(yC))

    figure('Name','Respuesta del Compensador', ...
           'NumberTitle','off');

    plot(tC,yC,'LineWidth',1.5);
    grid on;
    box on;

    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    title('Respuesta al escalon del compensador C(s)');
end

%% ==============================================================
% 15. Resultados que dejo en la ventana de comandos

fprintf('\n====================================================\n');
fprintf('                 RESULTADO FINAL\n');
fprintf('====================================================\n');

fprintf('\n1) Planta:\n');
fprintf('G(s) = ');
imprimirFraccion(numG,denG);
fprintf('\n');

fprintf('\n2) Polos originales:\n');
disp(p);

fprintf('\n3) Ceros originales:\n');
if isempty(z)
    fprintf('    No hay ceros.\n');
else
    disp(z);
end

fprintf('\n4) Polos principales deseados seleccionados por el usuario:\n');
disp(polosDeseados);

fprintf('\n5) Polos adicionales utilizados por el compensador dinamico:\n');
if isempty(polosAux)
    fprintf('    No se requieren polos adicionales.\n');
else
    disp(polosAux);
end

fprintf('\n6) Nueva ecuacion caracteristica deseada:\n');
imprimirPolinomio(Pd);
fprintf(' = 0\n');

fprintf('\n7) Compensador dinamico:\n');
fprintf('C(s) = B(s)/A(s) = ');
imprimirFraccion(B,A);
fprintf('\n');

fprintf('\n8) Ecuacion caracteristica compensada:\n');
imprimirPolinomio(charComp);
fprintf(' = 0\n');

fprintf('\n9) Polos obtenidos del sistema compensado:\n');
disp(polosComp);

fprintf('\n10) Funcion de transferencia compensada:\n');
fprintf('T(s) = ');
imprimirFraccion(numT,denT);
fprintf('\n');

fprintf('\n11) Verificacion de colocacion de polos:\n');
fprintf('    Error maximo de los polos principales = %.3e\n',errorPolos);
fprintf('    Error maximo de todos los polos (incluyendo auxiliares) = %.3e\n', ...
        errorPolosTotal);

if errorPolos < 1e-7
    fprintf('    RESULTADO: los polos principales fueron colocados correctamente.\n');
else
    fprintf('    ADVERTENCIA: revise la seleccion de polos.\n');
end

fprintf('\n====================================================\n');
fprintf('             PROYECTO TERMINADO\n');
fprintf('====================================================\n');

%% ==============================================================
% Funciones auxiliares

function a = quitarZeros(a)
    a = a(:).';

    if isempty(a)
        a = 0;
        return;
    end

    escala = max(1,max(abs(a)));

    while length(a)>1 && abs(a(1)) < 1e-10*escala
        a(1)=[];
    end

    a(abs(a)<1e-12*escala)=0;
end

function [a,b] = igualarLongitud(a,b)
    a=a(:).';
    b=b(:).';

    if length(a)<length(b)
        a=[zeros(1,length(b)-length(a)),a];
    elseif length(b)<length(a)
        b=[zeros(1,length(a)-length(b)),b];
    end
end

function ok = sonConjugados(v)
    v=v(:).';

    if isempty(v)
        ok=true;
        return;
    end

    ok=true;
    usados=false(size(v));

    for i=1:length(v)

        if usados(i)
            continue;
        end

        if abs(imag(v(i)))<1e-10
            usados(i)=true;
        else
            candidatos=find(~usados & abs(v-conj(v(i)))<1e-8);
            if isempty(candidatos)
                ok=false;
                return;
            end
            usados(i)=true;
            usados(candidatos(1))=true;
        end
    end
end

function rr=ordenarRaices(rr)
    rr=rr(:);

    % Primero por parte real y luego por parte imaginaria.
    M=[real(rr),imag(rr)];
    [~,idx]=sortrows(M,[1 2]);
    rr=rr(idx);
end

function rr=asociarRaices(rr,prev)
    rr=rr(:);
    prev=prev(:);

    if length(rr)~=length(prev)
        rr=ordenarRaices(rr);
        return;
    end

    usados=false(length(rr),1);
    salida=zeros(size(rr));

    for i=1:length(prev)

        dist=abs(rr-prev(i));
        dist(usados)=inf;

        [~,idx]=min(dist);

        salida(i)=rr(idx);
        usados(idx)=true;
    end

    rr=salida;
end

function rr=sortComplejo(rr)
    rr=rr(:);
    [~,idx]=sortrows([real(rr),imag(rr)],[1 2]);
    rr=rr(idx);
end

function rr=quitarRaicesInvalidas(rr)
    if isempty(rr)
        return;
    end

    rr=rr(isfinite(rr));
end

function [A,B,err] = calcularCompensador(D,N,P)
    % Resuelve:
    %
    % D*A + N*B = P
    %
    % con:
    % A(s)=s^(n-1)+a_(n-2)s^(n-2)+...+a0
    % B(s)=b_(n-1)s^(n-1)+...+b0
    %
    % D es de grado n.

    D=quitarZeros(D);
    N=quitarZeros(N);
    P=quitarZeros(P);

    n=length(D)-1;

    if length(P)~=2*n
        error('Grado de P_d incorrecto para la sintesis del compensador.');
    end

    % Variables:
    % x = [a0 ... a_(n-2) b0 ... b_(n-1)]'
    nA=n-1;
    nB=n;
    nVar=nA+nB;

    % Todos los polinomios se manejan en orden descendente.
    % Solo necesitamos igualar los coeficientes de grado 0...2n-2.
    M=zeros(2*n-1,nVar);

    % Parte fija de D*A cuando el primer coeficiente de A es 1.
    Afijo=[1 zeros(1,n-1)];
    fijo=conv(D,Afijo);

    fijo=padLeft(fijo,2*n);

    rhs=P-fijo;

    % Terminos que corresponden a los coeficientes de A.
    col=1;

    for k=0:n-2

        Ak=zeros(1,n);
        Ak(end-k)=1;

        base=conv(D,Ak);
        base=padLeft(base,2*n);

        M(:,col)=base(2:end).';
        col=col+1;
    end

    % Terminos que corresponden a los coeficientes de B.
    for k=0:n-1

        Bk=zeros(1,n);
        Bk(end-k)=1;

        base=conv(N,Bk);
        base=padLeft(base,2*n);

        M(:,col)=base(2:end).';
        col=col+1;
    end

    rhs2=rhs(2:end).';

    % Resuelvo el sistema de ecuaciones.
    if rcond(M)<1e-12
        warning(['La matriz de sintesis esta mal condicionada. ', ...
                 'Se utilizara una solucion por minimos cuadrados.']);
    end

    x=M\rhs2;

    % Con los coeficientes obtenidos construyo A y B.
    A=zeros(1,n);
    A(1)=1;

    for k=0:n-2
        A(end-k)=x(k+1);
    end

    B=zeros(1,n);

    for k=0:n-1
        B(end-k)=x(n-1+k+1);
    end

    % Reviso que la ecuacion diofantina se cumpla numericamente.
    DA=conv(D,A);
    NB=conv(N,B);

    [DAp,NBp]=igualarLongitud(DA,NB);
    suma=DAp+NBp;
    suma=padLeft(suma,length(P));

    err=max(abs(suma-P));
end

function a=padLeft(a,L)
    a=a(:).';

    if length(a)<L
        a=[zeros(1,L-length(a)),a];
    elseif length(a)>L
        a=a(end-L+1:end);
    end
end

function D=DnormSeguro(D)
    D=quitarZeros(D);
    D=D/D(1);
end

function N=NnormSeguro(N)
    N=quitarZeros(N);
end

function imprimirPolinomio(a)

    a=quitarZeros(a);
    grado=length(a)-1;

    texto='';

    for i=1:length(a)

        c=a(i);
        g=grado-(i-1);

        if abs(c)<1e-10
            continue;
        end

        if isempty(texto)

            if c<0
                signo='-';
            else
                signo='';
            end

        else

            if c>=0
                signo=' + ';
            else
                signo=' - ';
            end
        end

        c=abs(c);

        if g==0
            parte=sprintf('%.6g',c);

        elseif g==1

            if abs(c-1)<1e-10
                parte='s';
            else
                parte=sprintf('%.6g*s',c);
            end

        else

            if abs(c-1)<1e-10
                parte=sprintf('s^%d',g);
            else
                parte=sprintf('%.6g*s^%d',c,g);
            end
        end

        texto=[texto signo parte]; %#ok<AGROW>
    end

    if isempty(texto)
        texto='0';
    end

    fprintf('%s',texto);
end

function imprimirFraccion(num,den)
    imprimirPolinomio(num);
    fprintf(' / (');
    imprimirPolinomio(den);
    fprintf(')');
end

function [t,y]=respuestaEscalon(num,den,T,Npts)
    % Respuesta al escalon unitario sin Control System Toolbox.
    % Integracion RK4 en forma canonica controlable.

    num=num(:).';
    den=den(:).';

    den=den/den(1);
    num=num/den(1);

    num=quitarZeros(num);
    den=quitarZeros(den);

    n=length(den)-1;

    t=linspace(0,T,Npts);

    if n==0
        y=(num(end)/den(end))*ones(size(t));
        return;
    end

    if length(num)>length(den)
        y=nan(size(t));
        return;
    end

    if length(num)<n+1
        num=[zeros(1,n+1-length(num)),num];
    end

    a=den(2:end);
    b=num;

    A=zeros(n);

    if n>1
        A(1:n-1,2:n)=eye(n-1);
    end

    A(n,:)=-fliplr(a);

    B=zeros(n,1);
    B(n)=1;

    D=b(1);

    C=b(2:end)-D*a;
    C=fliplr(C);

    dt=t(2)-t(1);
    x=zeros(n,1);

    y=zeros(size(t));

    u=1;

    for k=1:length(t)

        y(k)=C*x+D*u;

        if k<length(t)

            f=@(xx) A*xx+B*u;

            k1=f(x);
            k2=f(x+0.5*dt*k1);
            k3=f(x+0.5*dt*k2);
            k4=f(x+dt*k3);

            x=x+(dt/6)*(k1+2*k2+2*k3+k4);
        end
    end
end
