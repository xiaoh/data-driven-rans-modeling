clear all;
Ub = 0.028;
for hill = 0:1
    fidObsL = fopen('obsCoord','w');
    fidObsV = fopen('obsVelocity','w');
    fidObsVPlot = fopen('obsVPlot', 'w');
    for iter = [1,2,3,4,5,6]
        fid = fopen(['Coord',num2str(iter)],'w');
        fid2 = fopen(['U',num2str(iter)],'w');
        dnsData = load(['x',num2str(iter),'.0_Ureal.xy']);

        M = flipud(dnsData);

        y = [M(end/2,1)];
        x = ones(size(y)) * iter;
        z = ones(size(y)) * 0.05;

        Coords = [x'; y'; z']';
        U = [M(end/2,2:4)];
       
        fprintf(fid, '%e %e %e\n',Coords');
        fprintf(fidObsL, '%e %e %e\n',Coords');
        fprintf(fidObsVPlot, '%e %e %e %e %e %e\n',[Coords';U']);
        fprintf(fid2,'%e\n',U');
        fprintf(fidObsV,'%e\n',U');
        fclose(fid);
        fclose(fid2);
    end
    for iter = [1,2,7,8]
        fid = fopen(['Coord',num2str(iter)],'w');
        fid2 = fopen(['U',num2str(iter)],'w');
        dnsData = load(['x',num2str(iter),'.0_Ureal.xy']);

        M = flipud(dnsData);

        y = [M(end-270,1);M(end-20,1)];
        x = ones(size(y)) * iter;
        z = ones(size(y)) * 0.05;

        Coords = [x'; y'; z']';
        U = [M(end-270,2:4);M(end-20,2:4)];
       
        fprintf(fid, '%e %e %e\n',Coords');
        fprintf(fidObsL, '%e %e %e\n',Coords');
        fprintf(fidObsVPlot, '%e %e %e %e %e %e\n',[Coords';U']);
        fprintf(fid2,'%e\n',U');
        fprintf(fidObsV,'%e\n',U');
        fclose(fid);
        fclose(fid2);
    end
%     for iter = [1]
%         fid = fopen(['Coord',num2str(iter)],'w');
%         fid2 = fopen(['U',num2str(iter)],'w');
%         dnsData = load(['x',num2str(iter),'.0_Ureal.xy']);
% 
%         M = flipud(dnsData);
% 
%         y = [M(end-270,1)];
%         x = ones(size(y)) * iter;
%         z = ones(size(y)) * 0.05;
% 
%         Coords = [x'; y'; z']';
%         U = [M(end-270,2:4)];
%        
%         fprintf(fid, '%e %e %e\n',Coords');
%         fprintf(fidObsL, '%e %e %e\n',Coords');
%         fprintf(fidObsVPlot, '%e %e %e %e %e %e\n',[Coords';U']);
%         fprintf(fid2,'%e\n',U');
%         fprintf(fidObsV,'%e\n',U');
%         fclose(fid);
%         fclose(fid2);
%     end
    fclose(fidObsL);
    fclose(fidObsV);
    fclose(fidObsVPlot);
end
