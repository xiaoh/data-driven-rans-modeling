load RANScoord;
load x0.0_U2.xy;

Ub = 0.028;
xq = RANScoord;
M = x0_0_U2;

x = flipud(M(:,1));
Udns = flipud(M(:,2:4))*Ub;
Urans = interp1(x,Udns,xq);

fid = fopen('Urans','w');
fprintf(fid, '(%f %f %f)\n',Urans');
fclose(fid);

hold on;
plot(x,Udns(:,1),'b','LineWidth',2);
plot(xq,Urans(:,1),'r');
hold off;