%IHN
clc
close all
clear all
A=20*20;
i=0;
W=10000;
lan=4000/(pi*1000^2)/200;
la1=lan;
Ndc=200; %device per cluster
Rc=200;
Rp      =     300;                           % Reporting period
Ps      =     10*8;                        % Bits
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
PtR=21-30;
% 2/2000 1/2000 1/5000  1/20000;

TauF=0.01;

Pr=0.7;
E0=1000; Est=.1;Pc=0.01; et=2;
Pt=1.5;
Pr=0.5;
for nW=10
    j=0;
    for din=[8:-1.5:5.5,5:-0.5:3,2.75:-0.25:1.5,1.4:-0.1:0.5]
        j=j+1;
        lap(j)=1/(pi*din^2);
        Ndpb=lan*(pi*(1000*din)^2);
        E=Pt*0.1+Pr; 
        %         E=Pt*min(1,TauF+Ndpb/Rp*0.001)+Pr+Pt*min(1,Ndpb/Rp*0.005); %A*lap(j)*c2*E+
        %         p(j)=  Fcov (din*1000,nW,W,lan*180,    Ndc,Rc,Rp,Ps);
        p1(j)=FcovT(din*1000,nW,W,la1,la1,Ndc,Rc,Rp,Ps,PtR,1,1,0,1,1 );
        CR=3.65/5;Rd=2*7.*W.*CR./(2.^7); Tt = floor(Ps./(Rd)*1000);                      % Transmission time
        cost1(j)=A*lap(j)+nW*W*c311+A*lap(j)*E*c21;
        if(1/p1(j)<8)
            Nre(j)=1/p1(j);
            LiF(j)=E0/(Est+(Nre(j)-1)*2*Est+Nre(j)*(Pc+et*10^((PtR)/10)*1)*Tt/1000 );
        else
            Nre(j)=nan;
            LiF(j)=E0/(Est+(8-1)*2*Est+8*(Pc+et*10^((PtR)/10)*1)*Tt/1000 );
        end
    end
end
yyaxis left
plot(lap,cost1*10,'-b*')
ylabel('Cost, Lifetime')
hold on
plot(lap,LiF ,'-g.')
yyaxis right
ylabel('Nr. of transmissions')
plot(lap, 1-p1 ,'--r')
% plot(lap,Nre ,'--k')
hold off
xlabel('Density of APs (\times Km^2)')
legend('Cost (\times 0.1 c_1)','Battery lifetime (  reporting period)','Outage probability')
 

