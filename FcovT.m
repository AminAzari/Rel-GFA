%IHN
function cove=FcovT(dis,nW,BW,la1,la2,Ndc,Rc,Rp,Ps,Po1,Po2,pow1,pow2,nt1,nt2 )
 si=100;
Pu=10^(Po1/10)*pow1;
Po=10^(Po2/10)*pow2;
sBW=nW*BW;

CR      =     3.65/5;                            % Code rate
Rd      =     2*7.*BW.*CR./(2.^7);         % Data rate
% Ps      =     10*8;                        % Bits
Tt      =     floor(Ps./(Rd)*1000);                      % Transmission time
% Rp      =     300;                           % Reporting period Ps      =     10*8;                        % Bits
Pot=10*log10(10^(Po1/10)*pow1);
Plu=133+38.3*log10(dis/1000);
PiGz=(10.^((Pot-Plu)/10));
No=10^(-20.4)*BW;
Th=1;
Pn=exp(-No*Th/PiGz );
del=3.83; 



lam1=(BW/sBW)*Ndc*la1*0.001*Tt*nt1/Rp; %5.3e-6
Ca=(2*(pi^2))/del*csc(2*pi/del); %pi/sinc(2/del)  %5.16
Pp1=exp(-lam1*(dis^2)*(Th^(2/del))*Ca);  %0.9534

lam2=(BW/sBW)*Ndc*0.001*Tt*nt1/Rp;
fun = @(r) r.*exp((-(r).^2)./(4*si^2))./(1+(r.^del)/((dis.^del)*(Th))) ;
H2=integral(fun,0,Inf)*2*pi/(4*pi*si^2);
Pp2=exp(-lam2*H2);

lam12=(BW/sBW)*Ndc*la2*0.001*Tt*nt2/Rp; %5.3e-6
Pp1=Pp1*(exp(-lam12*(dis^2)*(((Po/Pu)*Th)^(2/del))*Ca));  %0.9534

cove=Pn*Pp1*Pp2;
%%-----
 

 
