function [dmsSK,LL]=FunCs(in,dis,disD, Nd1,Nd2,Ra,BW,sBW,Rp,CR,Ps,PtR1,PtR2,kt,Rc,Ndc1,Ndc2)
Nsm=1;
for cont=1:Nsm
    % % Parameters
    %     P0      =     0;
    %     Nd1      =     4500;                         % # of devices
    %     Nd2=4500;
    %     Ra      =     3000;                          % Radius
    %     BW      =     10000;
    %     sBW     =     10*10000;                       % System BW
    Nc1=floor(Nd1/Ndc1); 
    Nc2=floor(Nd2/Ndc2);
    Nd1=Nc1*Ndc1;
    Nd2=Nc2*Ndc2;
    Nc1=floor(Nd1/Ndc1);
    Nc2=floor(Nd2/Ndc2);
    
    BWi     =     BW/5;
    BWp     =     1:sBW/1000;
    BWpn     =     length(BWp);
    Action     =     BW/2:BWi:sBW-BW/2;
    Acn=length(Action); 
    %     Rp      =     300;                           % Reporting period
    %     CR      =     3.65/5;                           % Code rate
    Rd      =     2*7.*BW.*CR./(2.^7);         % Data rate
    %     Ps      =     10*8;                        % Bits
    Tt1      =     floor(Ps./(Rd)*1000);                      % Transmission time
    Tt2      =     kt*floor(Ps./(Rd)*1000);                      % Transmission time
    %     PtR1     =     [21]-30;                          % Trans Pow. range
    %     PtR2     =     [21]-30;                          % Trans Pow. range
    IPAI    =     2;                            % Inv. Pow. Amp. Eff.
    %     Pc      =     0.01;                         % Circuit power
    %     PT      =     (IPAI*(10.^(PtR/10))/1000+Pc);  % Total power Cons.
    SimT    =     5000;     % Symulation time
    Np      =     floor(SimT/Rp);
    
    DataI   =     zeros(Nd1+Nd2,Np);
    Th=  1;%3 dB
    
    %     Nc=25;
%     Nc1=Nd1/Ndc1;
%     Nc2=Nd2/Ndc2;
    %     Rc=50;
    %  figure(1)
    ij=0;
    clear x y Rfu1
    for n=1:Nc1
        Rfc=(Ra-Rc/2)*sqrt(rand(1,1));
        % Rfc=0;
        Phc=360*rand(1,1);
        for m=1:Ndc1
            ij=ij+1;
            Rfcd=Rc*sqrt(rand(1,1)); 
            Phcd=360*rand(1,1);
            x(ij)=(Rfc*cos(Phc)+Rfcd*cos(Phcd));
            y(ij)=(Rfc*sin(Phc)+Rfcd*sin(Phcd));
            Rfu1(ij)=sqrt((Rfc*cos(Phc)+Rfcd*cos(Phcd))^2+(Rfc*sin(Phc)+Rfcd*sin(Phcd))^2);
        end
    end
%         scatter(x,y)
    %     hold on
    ij=0;  
    clear x y Rfu2
    for n=1:Nc2
        Rfc=(Ra-Rc/2)*sqrt(rand(1,1));
        % Rfc=0;
        Phc=360*rand(1,1);
        for m=1:Ndc2
            ij=ij+1;
            Rfcd=Rc*sqrt(rand(1,1));
            Phcd=360*rand(1,1);
            x(ij)=(Rfc*cos(Phc)+Rfcd*cos(Phcd));
            y(ij)=(Rfc*sin(Phc)+Rfcd*sin(Phcd));
            Rfu2(ij)=sqrt((Rfc*cos(Phc)+Rfcd*cos(Phcd))^2+(Rfc*sin(Phc)+Rfcd*sin(Phcd))^2);
        end
    end
    
    %     scatter(x,y,'r')
    
    %           Rfu=Ra*sqrt(rand(1,Nd));
    Rfu=[Rfu1,Rfu2];
    In=[ones(1,Nd1),2*ones(1,Nd2)];
    [Rf,sI]=sort(Rfu);
    sIn=In(sI);
    Dt=[];
    clear DataI 
    for j=1:Nd1+Nd2
        DataI(j,:)=sort(randperm(SimT*1000,Np));
        A=DataI(j,:);
        for ij=2:Np
            if(A(ij)-A(ij-1)<2*max(Tt1,Tt2))
                A(ij)=A(ij-1)+2*max(Tt1,Tt2);
            end
        end
        DataI(j,:)=A;
        Dt=[Dt,DataI(j,:)];
        
    end
    Rf=sort(Rf);
    sDt=sort(Dt);
    HoQ=DataI(:,1);
    CoI=ones(1,Nd1+Nd2);
    
    % %     load dataee.mat
    %     EC=PT*Tt/1000;
    
    PI=zeros(BWpn,SimT*1000);
    
    HoQ=DataI(:,1);
    CoI=ones(1,Nd1+Nd2);
    
    PI=zeros(BWpn,(SimT+1000)*1000);
    AKi=zeros(Nd1+Nd2,Np);
%     Tri=zeros(Nd1+Nd2,Nd1+Nd2);
    ACT=zeros(Nd1+Nd2,Np);
    co=0;
    AK=0;
    clear AK AKi ACT 
    cop=0;
    for j=1:floor(0.5*Np*(Nd1+Nd2))
        if(mod (j,500)==0)
%             [cont,  j/(0.6*Np*(Nd1+Nd2))]
        end
        fH=find(HoQ==sDt(j));
        if(isempty(fH))
            fH1=randperm(Nd1+Nd2,1);
            cop=cop+1;
            cop
        else
            fH1=fH(1);
        end
        if(CoI(fH1)==1)
            aci=randperm(Acn,1);
        elseif(CoI(fH1)>1)
            N=10^(-20.4)*BW;
            co=co+1;
            
            AK(Jj(fH1))=abs(Ts(fH1)*0.001* PIs(fH1)/...
                (  Ts(fH1)*0.001*N+...
                sum(sum( 0.001*PI(Action(acs(fH1))/1000-5+1:...
                Action(acs(fH1))/1000+5,...
                Tss(fH1):Tss(fH1)+Ts(fH1)-1)  ))...
                -Ts(fH1)*0.001*(PIs(fH1)) ))...
                >Th;
            
            
            
            
            AKi(fH1,CoI(fH1)-1)=AK(Jj(fH1));
            ACT(fH1,CoI(fH1)-1)=acs(fH1);
            %                 ECi=EC;
            aci=randperm(Acn,1);
            
        end
        %             Pot=(10.^(PtR/10))/1000;
         
        if(sIn(fH1)==1)
            Pot=PtR1;
            tau=Tt1;
        elseif(sIn(fH1)==2)
            Pot=PtR2;
            tau=Tt2;
        end
        
        
        
        %             f=868;hB=10;hM=1;
        %             CH= 0.8 + (1.1*log10(f - 0.7))*hM - 1.56*log10(f);
        %             Plu=P0+69.55 + 26.16*log10(f) - 13.82*log10(hB) - CH + (44.9 - 6.55*log10(hB))*log10(Rf(fH1)/1000) +sqrt(6)*randn(1,1)+raylrnd(1);
        Plu=133+38.3*log10(Rf(fH1)/1000)+raylrnd(1);
        PI(Action(aci)/1000-5+1:Action(aci)/1000+5,DataI(fH1,CoI(fH1)):...
            DataI(fH1,CoI(fH1))+tau-1)=...
            PI(Action(aci)/1000-5+1:Action(aci)/1000+5,DataI(fH1,CoI(fH1)):...
            DataI(fH1,CoI(fH1))+tau-1)+...
            (10.^((Pot-Plu)/10))/10*ones(10,tau);
         
        PIs(fH1)=10^((Pot-Plu)/10);
        Ts(fH1)=tau;
        Tss(fH1)=DataI(fH1,CoI(fH1));
        acs(fH1)=aci;
        Jj(fH1)=j;
        if(CoI(fH1)<Np)
            CoI(fH1)=CoI(fH1)+1;
            HoQ(fH1)=DataI(fH1,CoI(fH1));
        end
        
    end
    
    
    
    
    
    
    Rf1=Rf;
    isIn= find(sIn~=in);
    Rf1(isIn)=100000;
    A=find(abs(Rf1-dis)<disD);
    L =length(A);
    AA=find(CoI(A)>1);
    LL (cont)=length(AA);
%     MC=min(CoI(A))-1;
    ij=0;
    clear sSK
    for k=1:length(Rf);
        co=min(10,CoI(k)-2);
        if(co>1)
%             [k,co,size(ACT(k,:)) ]
            A3=find(ACT(k,1:co)>5 & ACT(k,1:co)<40);
            if(length(A3)>1)
                sSK(k)=sum(AKi(k,A3))/length(A3);
            else
                sSK(k)=nan;
            end
        else
            sSK(k)=nan;
        end
    end
    dmsSK(cont)=nanmean(sSK(A));
    
    
    
end

% pr=nanmean(dmsSK);
% Lr=nanmean(LL);
