%IHN
clc
close all
clear all
A=20*20;
i=0;
W=10000;
lan=2000/(pi*1000^2)/200;
Ndc=200; %device per cluster
Rc=200;
PtR=21-30;
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
% 2/2000 1/2000 1/5000  1/20000;

TauF=0.01;


E0=1000; Est=.1;Pc=0.01; et=2;

CovR=[0.3:0.1:0.9,0.91:0.01:0.99];

ik=0;
rw=nan(1,length(CovR));
for CovRi=CovR
    ik=ik+1;
    j=0;
    for nW=1:20
        din=2.4;
        j=j+1;
        lap(j)=1/(pi*din^2);
        %         p(j)=  Fcov(din*1000,nW,W,lan,Ndc,Rc,Rp,Ps);
        p(j)=FcovT(din*1000,nW,W,lan ,lan ,Ndc,Rc,Rp,Ps,PtR,1,1,0,1,1 );
        if(p(j)>CovRi)
            rw(ik)=nW;
            break
        end
    end
end



ik=0;
rd=nan(1,length(CovR));
for CovRi=CovR
    ik=ik+1;
    nW=10;
    j=0;
    for din=[8:-1.5:5.5,5:-0.5:3,2.75:-0.25:1.5,1.4:-0.1:0.5]
        j=j+1;
        lap(j)=1/(pi*din^2);
        p(j)=FcovT(din*1000,nW,W,lan ,lan ,Ndc,Rc,Rp,Ps,PtR,1,1,0,1,1 );
        if(p(j)>CovRi)
            rd(ik)=lap(j);
            break
        end
    end
    
    
end



ik=0;
rr=nan(1,length(CovR));
for CovRi=CovR
    ik=ik+1;
    nW=10;
    j=0;
    din=2.4;
    for ret=1:30
        j=j+1;
        lap(j)=1/(pi*din^2);
        p1(j)=FcovT(din*1000,nW,W,lan ,lan ,Ndc,Rc,Rp,Ps,PtR,1,1,0,1,1 );
        p(j)=funP(p1(j),ret);
        
        if(p(j)>CovRi)
            rr(ik)=ret;
            break
        end
    end
    
    
    
end





ik=0;
rp=nan(1,length(CovR));
for CovRi=CovR
    ik=ik+1;
    nW=10;
    j=0;
    din=2.4;
    for pow1=1:30
        j=j+1;
        lap(j)=1/(pi*din^2);
        p1(j)=FcovT(din*1000,nW,W,lan ,lan ,Ndc,Rc,Rp,Ps,PtR,1,pow1,0,1,1 );
        ret=1;
        p(j)=funP(p1(j),ret);
        
        if(p(j)>CovRi)
            rp(ik)=pow1; 
            ik
            break
        end
    end
    
    
    
end




yyaxis left
plot(CovR,rd,'-b','LineWidth',1.5)
ylabel('Required AP density (\times Km^{-2})')
hold on
yyaxis right
ylabel('Required BW,  Nr transmissions, transmit power')
plot(CovR,rw,'-k','LineWidth',1.5)
plot(CovR,rr,'--r','LineWidth',1.5)
plot(CovR,rp,'--bd','LineWidth',1.5)
hold off
xlabel('Required coverage probability')
legend('Required AP density (\times Km^{-2})','Required BW (\times 10 KHz)','Required Nr transmissions','Required transmit power')

