clear all;
load x0.0_U2.xy;
M = x0_0_U2;
index = 234;

Uup = M(1:index,2);
yup = M(1:index,1);

Ulow = M(index:end,2);
ylow = M(index:end,1);

UupAvg = trapz(yup',Uup')/(yup(end)-yup(1))*0.028

UlowAvg = trapz(ylow',Ulow')/(ylow(end)-ylow(1))*0.028