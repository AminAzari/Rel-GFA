%IHN
clc
close all
clear all
A=20*20;
i=0;
W=10000;
% lan=9000/(pi*3000^2);
lan=3000/(pi*1000^2)/200 ;
la1=1*lan;



Ndc=200; %device per cluster
Rc=200;
Rp      =     300;                           % Reporting period 
Ps      =     10*8;                         % Bits
PtR=21-30;
% c1=2000; %euro per bs per year-->0.02 meuro /bs annual
% c2=0.876;%euro-->euro per watt per year--> eur876/kWh annual
% c3=210000000/190000000; %euro per Hz --> 210 meuro per 190 Mhz
% c3=166000000/10000000;
% c31=16/2000;
% c31=1/20000;
c311=1/2000;
c312=1/4000;
c313=1/10000;
c21=0.876/2000;

% 2/2000 1/2000 1/5000  1/20000; 

% TauF=0.01;
Pt=1.5;
Pr=0.5;
A1=logspace(log10(10.5),log10(0.5),64);%[10.5:-.15:5,5:-0.05:3,2.75:-0.025:1.5,1.4:-0.01:0.79]
for nW=1:.3:20
    i=i+1;
    j=0; 
    for din=A1
        j=j+1; 
        lap(j)=1/(pi*din^2);
        Ndpb=lan*(pi*(1000*din)^2);
        E=Pt*0.1+Pr; 
%         p(i,j)= Fcov(din*1000,nW,W,lan,    Ndc,Rc,Rp,Ps,PtR);
        p(i,j)=FcovT(din*1000,nW,W,la1,la1,Ndc,Rc,Rp,Ps,PtR,1,1,0,1,1 );
         
        cost1(i,j)=A*lap(j)+nW*W*c311+A*lap(j)*E*c21;
        cost2(i,j)=A*lap(j)+nW*W*c312+A*lap(j)*E*c21; 
        cost3(i,j)=A*lap(j)+nW*W*c313+A*lap(j)*E*c21;
    end    
end



pTh=0.7;
ocov=nan(length(lap));
In=p>=pTh;
ocov(In)=p(In);

ocov1=nan(length(lap));
In1=p<pTh;
ocov1(In1)=p(In1); 

icost2New=nan(length(lap));
icost2New(In)=cost2(In);
icost2New1=nan(length(lap));
icost2New1(In1)=cost2(In1);
 

[X,Y]=meshgrid(lap,[1:.3:20]);
figure(1)
surf(X,Y, ocov1)
hold on
surf(X,Y, ocov)  %pTh*ones(length(lap))
grid on
xlabel('Density of BSs (Km^{-2})')
ylabel('Bandwidth (Khz)')
zlabel('Coverage probability')
view(3)
F=nan(1,length(lap));
for i=1:length(lap)
    for j=1:length(lap)
        if(p(i,j)>=pTh)
            F(i)=j;
            break
        end
    end
end





icost1=nan(length(lap));
for i=1:length(lap)
    if(F(i)>0)
        icost1(i,F(i):end)=cost1(i,F(i):end);
    end
end

icost2=nan(length(lap));
for i=1:length(lap)
    if(F(i)>0)
        icost2(i,F(i):end)=cost2(i,F(i):end);
    end
end

icost3=nan(length(lap));
for i=1:length(lap)
    if(F(i)>0)
        icost3(i,F(i):end)=cost3(i,F(i):end);
    end
end

figure(2)
% surf(X,Y, icost1)
% grid on
% hold on
surf(X,Y, icost2New)
hold on
surf(X,Y, icost2New1)
% surf(X,Y, icost3)
% hold off
xlabel('Density of BSs (Km^{-2})')
ylabel('Bandwidth (Khz)')
zlabel('Cost (\times c_1)')
view(3)
% [m1,sah]=nanmin(icost1);
% [~ ,so1]=nanmin(m1);
% sa1=sah(so1);
[m2,sah]=nanmin(icost2);
[~ ,so2]=nanmin(m2);
sa2=sah(so2);
% [m3,sah]=nanmin(icost3);
% [~ ,so3]=nanmin(m3);
% sa3=sah(so3);
% [sa1,so1;sa2,so2;sa3,so3]
[sa2,so2]