%IHN
clc
close all
clear all

TR=1:40;
Po1=5-30;
Po2=0;
% kl=50;

A=20*20;
i=0;
W=10000;
lan=100/(pi*1000^2)/200;
Ndc=200; %device per cluster
Rc=200;
Rp      =     300;                           % Reporting period
Ps      =     10*8;                        % Bits
la1=1*lan;
la2=1*lan;

nW= 10;
din=2.4;


E0=1000; Est=.1;Pc=0.01; et=2;
CovRi=0.5;

ik1=0;
rr=nan(1,length(TR));

for tr=TR
    
    ik1=ik1+1;
    
    pow1=tr;
    
    ik=0;
    for tr1=TR
        
        ik=ik+1;
        pow2=0;nt1=tr1;nt2=1;
        p31(ik)=FcovT(din*1000,nW,W,la1,la2,Ndc,Rc,Rp,Ps,Po1,Po2,pow1,pow2,nt1,nt2 );
        p3(ik)=funP(p31(ik),nt1);
        CR=3.65/5;Rd=2*7.*W.*CR./(2.^7); Tt = floor(Ps./(Rd)*1000);                      % Transmission time
        Nre(ik)=(1/p3(ik))*nt1;
        
        if(p3(ik)>CovRi)
            rr(ik1)=tr1;
            L(ik1)=E0/(Est+((1/p3(ik))-1)*2*Est+Nre(ik)*(Pc+et*pow1*10^((Po1)/10)*1)*Tt/1000);
            break
            

        end
    end
    
    
end
yyaxis left
plot(rr)
yyaxis right
plot(L)