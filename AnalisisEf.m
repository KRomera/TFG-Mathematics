function [x,xnorm,E,Enorm]=AnalisisEf(alpha,N,n,ej)
%Esta funcion devuelve el valor de la eficiencia
%alpha =orden de derivacion, Por otro lado se debe indicar
% el maximo n o N tal que dx=2^(-n) o dt=2^(-N). El analisis
% va desde dx=2^(-2) hasta dx=2^(-n) con dt=1/N o 
% de dt=2^(-2) hasta dt=2^(-N) con dx=1/n

%Para ello se introduce m:
m=min(n,N);
if N>n
    dt=1/N;
    paso='dx';
else
    dx=1/n;
    paso='dt';
end


%Variables para registrar los errores
E=zeros(1,m-1); Enorm=zeros(1,m-1); Time=zeros(1,m-1);
H=zeros(1,m-1); X=zeros(1,m-2); Xnorm=zeros(1,m-2);

% Hay dos posibles ejemplos con los que verificar
if ej==1
    sigma=0.25; r=0.05; K=0; T=1; xmin=0; xmax=1; xval=0.01;
elseif ej==2
    sigma=sqrt(2); r=0.5; K=0; T=1; xmin=0; xmax=1; xval=0.01;
end

for i=2:m
    H(i-1)=2^(-i); h=H(i-1);
    if N>n
        tic
        [E(i-1),Enorm(i-1)]=Implicito(alpha,sigma,r,K,T,xmin,xmax,dt,h,xval,ej);
        Time(i-1)=toc;
    else
        tic
        [E(i-1),Enorm(i-1)]=Implicito(alpha,sigma,r,K,T,xmin,xmax,h,dx,xval,ej);
        Time(i-1)=toc;
    end
    if i>2
        X(i-2)= log2(E(i-2)/E(i-1));
        Xnorm(i-2)=log2(Enorm(i-2)/Enorm(i-1));
    end
end
Np=1./H;
%x=mean(X);
x=X';
%xnorm=mean(Xnorm);
xnorm=Xnorm';
E=E';
Enorm=Enorm';
%Representaciones graficas:
figure
subplot(2,2,1);
    loglog(H,E,'-*')
    title(['Error en funcion del paso ', paso])
subplot(2,2,2);
    loglog(H,Enorm,'-*')
    title(['Error ||\cdot||_{dx} en funcion del paso ', paso])
subplot(2,2,3); 
    loglog(Np,Time,'-*')
    title('Tiempo computacional vs pasos n \cdot N')
    
%Error vs tiempo
subplot(2,2,4);
loglog(Time,E,'-*');
title('Error vs Tiempo computacional')
    
end