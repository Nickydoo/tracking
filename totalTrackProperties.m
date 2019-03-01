function [speed,vx,vy,averageSpeed,stdSpeed,lengthTrack,d,dtot,dnet,dmax,MSD,MI,OR,vAC,alpha,phi,DA] = totalTrackProperties(tracksIn, nbTracks, dt)

iout=1;
for iPos = 1:length(nbTracks)
for iNbTr = 1:nbTracks(iPos)
    iTrack = tracksIn{iPos}(:,end) == iNbTr;
    track = tracksIn{iPos}(iTrack,1:2);
    track(:,1:2) = track(~isnan(track(:,1)),1:2);
    
    x = track(:,1);
    y = track(:,2);
        d{iout} = sqrt(diff(x).^2+diff(y).^2); %distance
        dtot(iout) = sum(d{iout})/length(x);
        dnet(iout) = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2)/length(x);
        dmax(iout) = max(d{iout});
        MSD(iout) = mean(d{iout}.^2); %mean square distance
        MI(iout) = dnet(iout)/dtot(iout); %meandering index
        OR(iout) = dmax/dtot; %outreach ratio
        
        vx{iout} = diff(x)/dt; %instantaneous velocity
        vy{iout} = diff(y)/dt;
        v = [vx{iout} vy{iout}];
        speed{iout} = sqrt(vx{iout}.^2+vy{iout}.^2);
        averageSpeed(iout) = mean(speed{iout});
        stdSpeed(iout) = std(speed{iout});
        
        vAC{iout} = v*v';
        alpha{iout} = atan(diff(y)./diff(x)); %global turning angle
        phi{iout} = acos(diag(vAC{iout},1)./(speed{iout}(1:end-1).*speed{iout}(2:end)));
        DA{iout} = cos(diff(alpha{iout})); %direction autocorrelation
        lengthTrack(iout) = length(track)*dt;
    iout = iout+1;
     clear iTrack track x y v
end
end
    