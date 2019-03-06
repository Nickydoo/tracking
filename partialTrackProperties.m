function [averageSpeed,stdSpeed,dtot,dnet,dmax,MSD,MI,OR] = partialTrackProperties(tracksIn, nbTracks, dt, period)

%colonne = # de période
%ligne = track
nbPts = period/dt +1;


iNbTr=1;

for iNbTr = 1:nbTracks
        iTrack = tracksIn(:,end) == iNbTr;
        track = tracksIn(iTrack,1:3);
    
        for idxAv = 1:length(track)-nbPts-1
            iNoNAN = [];
            iPeriod = idxAv:nbPts+(idxAv-1);
            iNoNAN = iPeriod(~isnan(track(iPeriod,1)));
            x = track(iNoNAN,1);
            y = track(iNoNAN,2);
            t = track(iNoNAN,3);
            d = sqrt(diff(x).^2+diff(y).^2); %distance

            vx = diff(x)./(diff(t)*dt); %instantaneous velocity
            vy = diff(y)./(diff(t)*dt);
            v = [vx vy];
            speed = sqrt(vx.^2+vy.^2);
         
            averageSpeed(iNbTr,idxAv) = mean(speed);
            stdSpeed(iNbTr,idxAv) = std(speed);
            if ~isempty(speed)
                MovMaxSpeed(iNbTr,idxAv) = max(speed);
                dmax(iNbTr,idxAv) = max(d);
                dnet(iNbTr,idxAv) = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2)/length(x);
                dtot(iNbTr,idxAv) = sum(d)/length(x);
            else
                MovMaxSpeed(iNbTr,idxAv) = NaN;
                dmax(iNbTr,idxAv) = NaN;
                dnet(iNbTr,idxAv) = NaN;
                dtot(iNbTr,idxAv) = NaN;
            end
            
            MSD(iNbTr,idxAv) = mean(d.^2); %mean square distance
            MI(iNbTr,idxAv) = dnet(iNbTr,idxAv)./dtot(iNbTr,idxAv); %meandering index
            OR(iNbTr,idxAv) = dmax(iNbTr,idxAv)./dtot(iNbTr,idxAv); %outreach ratio
        end
            
     clear iTrack track x y v d speed vx vy
end
end
    