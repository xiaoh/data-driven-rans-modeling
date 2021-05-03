for iter = 0:8
    RANScoord = load(['RANScoord',num2str(iter)]);
    data = load(['x',num2str(iter),'.0_normalized_Results.xy']);

    Ub = 0.028;
    xq = RANScoord(:,1);
    M = data;

    x = M(:,1);
    uudns = M(:,4)*Ub^2;
    uurans = interp1(x,uudns,xq,'spline');

    vvdns = M(:,5)*Ub^2;
    vvrans = interp1(x,vvdns,xq,'spline');

    wwdns = (M(:,7)*2-M(:,4)-M(:,5))*Ub^2;
    wwrans = interp1(x,wwdns,xq,'spline');

    uvdns = M(:,6)*Ub^2;
    uvrans = interp1(x,uvdns,xq,'spline');

    zeroR = zeros(size(uurans));
    zeroD = zeros(size(uudns));

    ransTau = [xq'; uurans'; uvrans'; zeroR'; vvrans'; zeroR'; wwrans']';
    dnsTau = [x'; uudns'; uvdns'; zeroD'; vvdns'; zeroD'; wwdns']';

    fid = fopen(['ransTauFine',num2str(iter)],'w');
    fprintf(fid, '%f %e %e %e %e %e %e\n',ransTau');
    fclose(fid);

    hold on;
    plot(x,vvdns(:,1),'b','LineWidth',2);
    plot(xq,vvrans(:,1),'r');
    hold off;
end