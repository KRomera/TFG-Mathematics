function [E,Enorm,U]=Implicito(alpha,sigma,r,K,T,xmin,xmax,dt,dx,tval,ej)

%Definimos el mallado:

% En tiempo, se representa $t_k$=t(k+1)
N = floor(T/dt);
t = 0:dt:T;

% En espacio, se representa $x_i$=x(i+1)
n = floor((xmax - xmin) / dx);
x = xmin : dx : xmax;

nval = floor(tval/dt)+1;
%Definimos las constantes en funciones de los parámetros de entrada:
a = sigma^2/2;
c = r;
b = c - a;
d = gamma(2-alpha) * dt^alpha;

%Términos constantes de la matriz;

a1 = a*d/ (dx^2);
a2 = b*d/(2*dx);

% Con estas constantes también construiremos G^k:
g0 = a1 - a2;
gk = a1 + a2;

%Extra para comprobación del orden de convergencia:
F=zeros(n+1,N+1);
if nargin > 10
    if ej==1
        fun = @(x,t) (2.*t.^(2-alpha)./gamma(3-alpha)+2.*t.^(1-alpha)./gamma(2-alpha)).*x.^2.*(1-x) - (t+1).^2.*(a.*(-6.*x+2)+b.*(-3.*x.^2+2.*x)-c.*x.^2.*(1-x));
    elseif ej==2
        fun = @(x,t) (2.*t.^(2-alpha)./gamma(3-alpha)+2.*t.^(1-alpha)./gamma(2-alpha)).*(x.^3+x.^2+1) - (t+1).^2.*(a.*(6.*x+2)+b.*(3.*x.^2+2.*x)-c.*(x.^3+x.^2+1));
    end
    for i=1:length(x)
        F(i,:)=fun(x(i),t);
    end
end

% Construimos la matriz:
e = ones(n-1,1);
A = spdiags([-g0 * e, (1 + 2 * a1 + c * d)*e, -gk*e],-1:1,n-1,n-1);

% Vamos a necesitar realmente su inversa para poder calcular:
%B = inv(A); % Matlab indica que no es eficiente computacionalmente



% b_j=(j+1)^(1-alpha)-j^(1-alpha), pero b_j=b_j(j+1):
b_j=zeros(1,N+1);

% Se representa $f_i=f(x_i)$=f(x(i+1))=f(i)
f=zeros(n-1,1);

% $q^k=q(t_k)$=q(t(k+1))=q(k+1)$
q=zeros(1,N+1);

for k=1:N+1
    q(k) = exp(xmax) - K* exp(-r* (T- t(k)));
    b_j(k) = k^(1-alpha)-(k-1)^(1-alpha);
end


% Como f=$(0,f_2,...,f_{n-2},q^0)^T$=(0,f(x(3)),...,f(x(n-1)), q(1))
for i=2:n-2
    f(i) = max (exp(x(i+1))-K,0);
end

f(n-1)=q(1);
% Definimos U, solución aproximada como una matriz que inicializamos el 0,
% establecemos entonces la condiciones iniciales y a continuación
% computamos el resto de valores:

U = zeros(n+1,N+1);
G = zeros(n-1,N+1);

%Condiciones iniciales:



if nargin < 11
    %U(1,:)=0;
    U(n+1,:)=q;
    U(2:n,1)=f;
else
    if ej==1
        %U(1,:)=0;
        %U(n+1,:)=0;
        U(:,1)=x.^2.*(1-x);
    elseif ej==2
        U(1,:)=(t+1).^2;
        U(n+1,:)=3.*(t+1).^2;
        U(:,1)=x.^3+x.^2+1;
    end
end

G(1,:) = g0 * U(1,:);
G(n-1,:) = gk * U(n+1,:);

% Computación del resto de valores:

% Inicio:
U(2:n,2) = A \ ( U(2:n,1)+G(:,2)+d*F(2:n,2) );

% Resto 1 <= k <= N-1
for k=2:N
    aux=0;
    for j=1:k-1
        aux =  aux + (b_j(j)-b_j(j+1))*U(2:n,1+k-j);
    end
    aux = aux + b_j(k)* U(2:n,1) + G(:,k+1)+d*F(2:n,k+1);
    U(2:n,k+1)= A \ aux;
end

if nargin < 11
    E=NaN;
    Enorm=NaN;
    %f1=figure;
    plot(exp(x),U(:,nval),'DisplayName',num2str(alpha))
    legend('-DynamicLegend','Location','northwest')
    xlabel('Z - Precio de la acción'), ylabel('C - Precio de la opción call')
    hold on
else
    if ej==1
        funsol= @(x,t) (t+1).^2.*x.^2.*(1-x);
    elseif ej==2
        funsol= @(x,t) (t+1).^2.*(x.^3+x.^2+1);
    end
    Usol=zeros(n+1,N+1);
    for i=1:length(x)
        Usol(i,:)=funsol(x(i),t);
    end
    E=max(max(abs(U-Usol)));
    Enorm=0;
    for k=1:N+1
        Enorm=max(Enorm,sqrt(sum(dx*(U(:,k)-Usol(:,k)).^2)));
    end
end

% [tt,X]=meshgrid(t,exp(x));
% f2=figure;
% s=surf(X,tt,U,'FaceAlpha',0.9);
% hold on
% s.EdgeColor = 'none';
% xlabel('Z - Precio de la acción'), ylabel('t - tiempo (años)'), zlabel('C - Precio de la opción call')
end