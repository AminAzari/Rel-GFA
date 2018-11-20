function P=funP(p11,nr)
p=0;
for i=0:nr-1
    p=p+(1-p11)^i*p11;
end
P=p;
