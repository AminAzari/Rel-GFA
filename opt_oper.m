%IHN
clc
close all
clear all
 
TR=1:10;
Po1=21-30;
Po2=21-30;
% kl=50;

A=20*20; 
i=0;
W=10000;
lan=1500/(pi*1000^2)/200;
Ndc=200; %device per cluster
Rc=200;
Rp      =     300;                           % Reporting period 
Ps      =     10*8;                        % Bits
la1=1*lan;
la2=1*lan;

nW= 10;
din=2.4;

 
E0=1000; Est=.1;Pc=0.01; et=2;
 
ik=0;
for tr=TR
     
    ik=ik+1; 
    pow1=tr;pow2=1;nt1=1;nt2=1;  
    p1(ik)=FcovT(din*1000,nW,W,la1,la2,Ndc,Rc,Rp,Ps,Po1,Po2,pow1,pow2,nt1,nt2 );
    p2(ik)=FcovT(din*1000,nW,W,la2,la1,Ndc,Rc,Rp,Ps,Po2,Po1,pow2,pow1,nt2,nt1 );
    CR=3.65/5;Rd=2*7.*W.*CR./(2.^7); Tt = floor(Ps./(Rd)*1000);                      % Transmission time
    Nre(ik)=(1/p1(ik));  
    LiF1(ik)=E0/(Est+(Nre(ik)-1)*2*Est+Nre(ik)*(Pc+et*10^((Po1)/10)*pow1)*Tt/1000 );
end
 
 
ik=0; 
for tr=TR
     
    ik=ik+1;
    pow1=1;pow2=1;nt1=tr;nt2=1;
    p31(ik)=FcovT(din*1000,nW,W,la1,la2,Ndc,Rc,Rp,Ps,Po1,Po2,pow1,pow2,nt1,nt2 );
    p3(ik)=funP(p31(ik),nt1);
    p4(ik)=FcovT(din*1000,nW,W,la2,la1,Ndc,Rc,Rp,Ps,Po2,Po1,pow2,pow1,nt2,nt1 );
    CR=3.65/5;Rd=2*7.*W.*CR./(2.^7); Tt = floor(Ps./(Rd)*1000);                      % Transmission time
    Nre(ik)=(1/p3(ik))*nt1;
    LiF2(ik)=E0/(Est+((1/p3(ik))-1)*2*Est+Nre(ik)*(Pc+et*10^((Po1)/10)*1)*Tt/1000);
end

 



la1=0.5*lan;
la2=lan;
 
ik=0;
for tr=TR
     
    ik=ik+1;
    pow1=tr;pow2=1;nt1=1;nt2=1;
    p5(ik)=FcovT(din*1000,nW,W,la1,la2,Ndc,Rc,Rp,Ps,Po1,Po2,pow1,pow2,nt1,nt2 );
    p6(ik)=FcovT(din*1000,nW,W,la2,la1,Ndc,Rc,Rp,Ps,Po2,Po1,pow2,pow1,nt2,nt1 );
    CR=3.65/5;Rd=2*7.*W.*CR./(2.^7); Tt = floor(Ps./(Rd)*1000);                      % Transmission time
    Nre(ik)=(1/p5(ik));
    LiF3(ik)=E0/(Est+(Nre(ik)-1)*2*Est+Nre(ik)*(Pc+et*10^((Po1)/10)*pow1)*Tt/1000);
end



ik=0;
for tr=TR
     
    ik=ik+1;
    pow1=1;pow2=1;nt1=tr;nt2=1;
    p71(ik)=FcovT(din*1000,nW,W,la1,la2,Ndc,Rc,Rp,Ps,Po1,Po2,pow1,pow2,nt1,nt2 );
    p7(ik)=funP(p71(ik),nt1);
    p8(ik)=FcovT(din*1000,nW,W,la2,la1,Ndc,Rc,Rp,Ps,Po2,Po1,pow2,pow1,nt2,nt1 );
    CR=3.65/5;Rd=2*7.*W.*CR./(2.^7); Tt = floor(Ps./(Rd)*1000);                      % Transmission time
    Nre(ik)=(1/p7(ik))*nt1;
    LiF4(ik)=E0/(Est+((1/p7(ik))-1)*2*Est+Nre(ik)*(Pc+et*10^((Po1)/10)*1)*Tt/1000);
end

 





figure(1)
yyaxis left
plot(TR,p1)
hold on
plot(TR,p3)
plot(TR,p5)
plot(TR,p7)
ylabel('Probability of success (Type 1)')
yyaxis right
plot(TR,p2)
plot(TR,p4)
plot(TR,p6)
plot(TR,p8)
ylabel('Probability of success (Type 2)')
hold off
legend('p_s(1,d_{eg}) vs. {\it P}_1; Sc1',...
    'p_s(1,d_{eg}) vs. {\it n}_1; Sc1',...
    'p_s(1,d_{eg}) vs. {\it P}_1; Sc2',...
    'p_s(1,d_{eg}) vs. {\it n}_1; Sc2',...
    'p_s(2,d_{eg}) vs. {\it P}_1; Sc1',...
    'p_s(2,d_{eg}) vs. {\it n}_1; Sc1',...
    'p_s(2,d_{eg}) vs. {\it P}_1; Sc2',...
    'p_s(2,d_{eg}) vs. {\it n}_1; Sc2')
xlabel('Nr of transmissions, Transmit power (\times 126 mW)')

figure(2)
plot(TR,LiF1)
hold on
plot(TR,LiF2)
plot(TR,LiF3)
plot(TR,LiF4)
hold off
xlabel('Nr of transmissions, Transmit power (\times 126 mW)')
ylabel('Battery lifetime (Type 1)')
legend('Lifetime vs. {\it P}_1 (Sc1)','Lifetime vs. {\it n}_1 (Sc1)','Lifetime vs. {\it P}_1 (Sc2)','Lifetime vs. {\it n}_1 (Sc2)')

