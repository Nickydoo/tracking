function [averageSpeed,stdSpeed,dtot,dnet,dmax,MSD,MI,OR] = partialTrackProperties(tracksIn, nbTracks, dt, period)

nbPts = period/dt;

iout=1;
for it = 1:length(nbTracks)
for idx = 1:nbTracks(it)
        iTrack = tracksIn{it}(:,end) == idx;
        track = tracksIn{it}(iTrack,1:2);
    
        x = track(:,1);
        y = track(:,2);
        d = sqrt(diff(x).^2+diff(y).^2); %distance
        
        vx = diff(x)/dt; %instantaneous velocity
        vy = diff(y)/dt;
        v = [vx vy];
        speed = sqrt(vx.^2+vy.^2);
        
        for idxAv = 1:length(track)-nbPts  
            averageSpeed{iout}(idxAv) = mean(speed(idxAv:idxAv+nbPts-1));
            stdSpeed{iout}(idxAv) = std(speed(idxAv:idxAv+nbPts-1));
            MovMaxVelocity{iout}(idxAv) = max(speed(idxAv:idxAv+nbPts-1));
            
            dtot{iout}(idxAv) = sum(d(idxAv:idxAv+nbPts-1));
            dnet{iout}(idxAv) = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);
            dmax{iout}(idxAv) = max(d(idxAv:idxAv+nbPts-1));
            MSD{iout}(idxAv) = mean(d(idxAv:idxAv+nbPts-1).^2); %mean square distance         
        end
            MI{iout} = dnet{iout}/dtot{iout}; %meandering index
            OR{iout} = dmax{iout}/dtot{iout}; %outreach ratio
            
    iout = iout+1;
     clear iTrack track x y v d speed vx vy
end
end
    