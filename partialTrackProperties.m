function [averageSpeed,stdSpeed,dtot,dnet,dmax,MSD,MI,OR] = partialTrackProperties(tracksIn, nbTracks, dt, period)

nbPts = period/dt;

iout=1;
for iPos = 1:length(nbTracks)
for iNbTr = 1:nbTracks(iPos)
        iTrack = tracksIn{iPos}(:,end) == iNbTr;
        track = tracksIn{iPos}(iTrack,1:3);
    
        x = track(:,1);
        y = track(:,2);
        t = track(:,3);
        d = sqrt(diff(x).^2+diff(y).^2); %distance
        
        vx = diff(x)/dt; %instantaneous velocity
        vy = diff(y)/dt;
        v = [vx vy];
        speed = sqrt(vx.^2+vy.^2);
        
        for idxAv = 1:length(track)-nbPts-1
            speed_tmp = speed(idxAv:idxAv+nbPts);
            averageSpeed{iout}(idxAv) = mean(speed_tmp(~isnan(speed_tmp));
            stdSpeed{iout}(idxAv) = std(speed(idxAv:idxAv+nbPts));
            MovMaxVelocity{iout}(idxAv) = max(speed(idxAv:idxAv+nbPts));
            
            dtot{iout}(idxAv) = sum(d(idxAv:idxAv+nbPts));
            dnet{iout}(idxAv) = sqrt((x(idxAv+nbPts)-x(idxAv))^2 + (y(idxAv+nbPts)-y(idxAv))^2);
            dmax{iout}(idxAv) = max(d(idxAv:idxAv+nbPts));
            MSD{iout}(idxAv) = mean(d(idxAv:idxAv+nbPts).^2); %mean square distance         
        end
            MI(iout) = dnet{iout}/dtot{iout}; %meandering index
            OR(iout) = dmax{iout}/dtot{iout}; %outreach ratio
            
    iout = iout+1;
     clear iTrack track x y v d speed vx vy
end
end
    