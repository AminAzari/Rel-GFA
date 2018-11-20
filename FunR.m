%IHN

function pr=FunR(in,dis,disD, Nd1,Nd2,Ra,BW,sBW,Rp,CR,Ps,PtR1,PtR2,kt,Rc,Ndc1,Ndc2)
% LLr=1;
k=0;
LLs=[];
Dd=zeros(1,2);
while(sum(Dd(:,2))<500)
    [dmsSK,LL]=FunCs(in,dis,disD, Nd1,Nd2,Ra,BW,sBW,Rp,CR,Ps,PtR1,PtR2,kt,Rc,Ndc1,Ndc2);
    if(LL>1)
        k=k+1;
        LLs=[LLs;LL];
        LLr=nanmean(LLs);
        Dd(k,:)=[dmsSK,LL];
    end
end
pr=sum(Dd(:,1).*Dd(:,2))/sum(Dd(:,2));


