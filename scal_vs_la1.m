%IHN
clc
% close all
clear all
 
TR=1:10;
Po1=21-30;
Po2=14-30;
% kl=50;

A=20*20; 
i=0;
W=10000;
lan=200/(pi*1000^2)/200;
Ndc=200; %device per cluster
Rc=200;
Rp      =     300;                           % Reporting period 
Ps      =     10*8;                        % Bits
 
la2=1*lan;

% PtR=21-30; 
 
E0=1000; Est=.1;Pc=0.01; et=2;
CovR=[0.5];

ik=0;
rw=nan(1,length(TR));
for tr=TR
     
    CovRi=CovR;
    ik=ik+1; 
    rw; 
    for nW=1:30
        j=0;
        din=2.4;
        j=j+1;
        lap(j)=1/(pi*din^2);
        Ndpb=lan*(pi*(1000*din)^2);
        %         E=Pt*min(1,TauF+Ndpb/Rp*0.001)+Pr+Pt*min(1,Ndpb/Rp*0.005); %A*lap(j)*c2*E+
        pow1=1;pow2=0;nt1=1;nt2=1;
        la1=tr*lan;
        p(j)=FcovT(din*1000,nW,W, la1,la2,Ndc,Rc,Rp,Ps,Po1,Po2,pow1,pow2,nt1,nt2 );
 
        if(p(j)>CovRi)
            rw(ik)=nW;
            break
        end
        
        
    end
end



ik=0;
rd=nan(1,length(TR));
for tr=TR
%     lan=kl*100/(pi*1000^2);
    CovRi=CovR;
    ik=ik+1;
    nW=10;
    j=0;
    for din=[8:-1.5:5.5,5:-0.5:3,2.75:-0.25:1.5,1.4:-0.1:0.5]
        j=j+1;
        lap(j)=1/(pi*din^2);
        Ndpb=lan*(pi*(1000*din)^2);
        %         E=Pt*min(1,TauF+Ndpb/Rp*0.001)+Pr+Pt*min(1,Ndpb/Rp*0.005); %A*lap(j)*c2*E+
        pow1=1;pow2=0;nt1=1;nt2=1;
        la1=tr*lan;
        p(j)=FcovT(din*1000,nW,W,la1,la2,Ndc,Rc,Rp,Ps,Po1,Po2,pow1,pow2,nt1,nt2);
        if(p(j)>CovRi)
            rd(ik)=lap(j);
            break
        end
    end
    
    
end 


 
ik=0;
rr=nan(1,length(TR));
for tr=TR
%     lan=kl*100/(pi*1000^2);
    CovRi=CovR;
    ik=ik+1;
    nW=10;
    j=0;
    din=2.4;
    for ret=1:30
        j=j+1;
        lap(j)=1/(pi*din^2);
        Ndpb=lan*(pi*(1000*din)^2);
        %         E=Pt*min(1,TauF+Ndpb/Rp*0.001)+Pr+Pt*min(1,Ndpb/Rp*0.005); %A*lap(j)*c2*E+
        pow1=1;pow2=0;nt1=ret;nt2=1;
        la1=tr*lan;
        p(j)=FcovT(din*1000,nW,W,la1,la2,Ndc,Rc,Rp,Ps,Po1,Po2,pow1,pow2,nt1,nt2);
         pp(j)=funP(p(j),ret);     
        
        if(pp(j)>CovRi)
            rr(ik)=ret;
            break
        end
    end
    
    
    
end






ik=0;
rp=nan(1,length(TR));
for tr=TR
%     lan=kl*100/(pi*1000^2);
    CovRi=CovR;
    ik=ik+1;
    nW=10;
    j=0;
    din=2.4;
    ret=1;
    for pow1=1:30
         pow2=0;nt1=1;nt2=1;
        j=j+1;
        lap(j)=1/(pi*din^2);
        Ndpb=lan*(pi*(1000*din)^2);
        la1=tr*lan;
        %         E=Pt*min(1,TauF+Ndpb/Rp*0.001)+Pr+Pt*min(1,Ndpb/Rp*0.005); %A*lap(j)*c2*E+
        p(j)=FcovT(din*1000,nW,W,la1,la2,Ndc,Rc,Rp,Ps,Po1,Po2,pow1,pow2,nt1,nt2);
        
        if(p(j)>CovRi)
            rp(ik)=pow1; 
            break 
        end
    end
end



figure(2)
yyaxis left
plot(TR,rd,'-b','LineWidth',1.5)
ylabel('Required AP density (\times Km^{-2})')
hold on
yyaxis right
ylabel('Required BW,  Nr transmissions, transmit power')
plot(TR,rw,'-k','LineWidth',1.5)
plot(TR,rr,'--r','LineWidth',1.5)
plot(TR,rp,'--bd','LineWidth',1.5)
hold off
xlabel('Density of nodes (\times \lambda_1\upsilon_1)')
legend('Required AP density (\times Km^{-2})','Required BW (\times 10 KHz)','Required Nr transmissions','Required transmit power')
grid on