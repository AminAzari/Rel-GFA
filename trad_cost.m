
%IHN
clc
clear all 
close all
% load('datan1.mat')
c1=1;c2=1/2000; 
for j=2:20
 
sBW=10000*j;
        gamma_th=1;
% % % % % % % % % %     Ra=20000;
% % % % % % % % %     lambda_a=1/(pi*1500^2);
    lambda_1=1000/(pi*1000^2)/200;
    BW=10000;
% % % % % % % % %     
    Ndc=200;
    Tt=1;
    Rp=300;
% % % % % % % % %     ups_21=(BW/sBW)*Ndc*Tt/Rp;
% % % % % % % % % %     
    Omega=1;
    N=10^(-20.4)*BW;
    P_1=10^((21-30)/10);
% % % % % % % % %     alpha=10^(-13.3);
% % % % % % % % % %     
% % % % % % % % %     Red=sqrt(1/lambda_a);
% % % % % % % % %     AA=lambda_1*ups_21*((gamma_th)/(Omega))^(0.5)*0.5*(pi^2 )*csc(pi/2);
% % % % % % % % %     
% % % % % % % % %     
% % % % % % % % %     Psm=exp(-ups_21)*exp(-((Red/1000)^2)*(N*gamma_th)/(Omega*P_1*alpha))*exp(-AA*Red^2);
% % % % % % % % %     Psm; %%maximum distance
% % % % % % % % % %     %%--------
% % % % % % % % %     b=pi*lambda_a;
% % % % % % % % %     syms r
% % % % % % % % %     Rea=2*b*double(limit(- (r*exp(-r^2*b))/(2*b) - (pi^(1/2)*erfi(r*(-b)^(1/2)))/(4*(-b)^(3/2)), r, Inf));
% % % % % % % % %     AA=lambda_1*ups_21*((gamma_th)/(Omega))^(0.5)*0.5*(pi^2)*csc(pi/2);
% % % % % % % % %     Psa=exp(-ups_21)*exp(-((Rea/1000)^2)*(N*gamma_th)/(Omega*P_1*alpha))*exp(-AA*Rea^2);
% % % % % % % % %     Psa; %%avg distance
% % % % % % % % % %     
% % % % % % % % %     syms r
% % % % % % % % %     % clear r
% % % % % % % % %     fun=@(r)(2*lambda_a*pi.*(r)).*exp(-lambda_a*pi*(r).^2)...
% % % % % % % % %         .*exp(-ups_21)...
% % % % % % % % %         .*exp(-((r/1000).^2).*(N*gamma_th)/(Omega*P_1*alpha))...
% % % % % % % % %         .*exp(-AA*(r).^2);
% % % % % % % % %     Psni=integral(fun,0,Inf);
% % % % % % % % % %     %%---
% % % % % % % % %     alpha=10^(-(133-60)/10);%kiloometer
% % % % % % % % %     ups_21=(BW/sBW)*Ndc*Tt/Rp;
% % % % % % % % %     AA=lambda_1*ups_21*((gamma_th)/(Omega))^(0.5)*0.5*(pi)*csc(pi/2);
% % % % % % % % %     Nu=0.5*sqrt(pi)*(lambda_a)*exp(-ups_21 );
% % % % % % % % %     Nd=AA+ lambda_a+(  N*gamma_th)/(Omega*P_1*alpha*pi);
% % % % % % % % %     Psi=Nu/Nd; %average ps from int
% % % % % % % % % %     %%%-----
% % % % % % % % %   
% % % % % % % % % %     %% ----
% % % % % % % % % %     
% % % % % % % % % %     a0=0.5*sqrt(pi)/Psi;
    w=sBW;
 alpha=10^(-(133-60)/10);%kiloometer
    a1=(BW)*Ndc*Tt/Rp; 
    a2=lambda_1*a1*((gamma_th)/(Omega))^(0.5)*0.5*(pi )*csc(pi/2);
    a3=(N*gamma_th)/(Omega*P_1*alpha*pi); 
% % % % % % % % % %     a3=0; 
% % % % % % % % % %     % floor([lambda_a,l_a2]*Ra^2)
% % % % % % % % % %      
% % % % % % % % %     Psu=0.5*sqrt(pi)*lambda_a*exp(-a1/sBW)/...
% % % % % % % % %         (a2/sBW+lambda_a+a3);
% % % % % % % % % %       [Psm,Psa,Psi,Psni,Psu]
    %%
    % clear all
    % syms Psu a1 a2 a3  lambda_a w
    % solve(Psu-0.5*sqrt(pi)*lambda_a*exp(-a1/w)/(a2/w+lambda_a+a3),lambda_a)
% %     save('datan.mat','a1','a3','a2','Psu')
Psu=0.4;  
    Nu=Psu*exp(a1/w)*(a3 + a2/w);
    Nd=-Psu*exp(a1/w) + 0.5*sqrt(pi);
    La(j)=Nu/Nd;  
    La1(j)=(a3 + a2/w)/(0.5*sqrt(pi)/Psu*exp(-a1/w)-1);
    co(j)=c1*La1(j)*10000*10000+c2*w;
end 
yyaxis left

plot(2:20,La1(2:end))
yyaxis right
plot(2:20,co(2:end))

%%
% syms Psu a1 a2 a3  lambda_a w
% Nu=(a3 + a2/w);
%     Nd=-1 + 0.5*sqrt(pi)/(Psu*exp(a1/w));
% solve(diff(a4*Nu/Nd+w)==0,w)