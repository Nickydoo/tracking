function [MovAverageVelocity,MovStdVelocity,MovMaxVelocity,velocity,vx,vy,averageVelocity,stdVelocity] = movingAverageVelocity(tracksIn, nbTracks, period, dt)

nbPts = period/dt;

for idx = 1:nbTracks
    iTrack = tracksIn(:,end) == idx;
    track = tracksIn(iTrack,1:2);
        vx{idx} = diff(track(:,1))/dt;
        vy{idx} = diff(track(:,2))/dt;
        velocity{idx} = sqrt(vx{idx}.^2+vy{idx}.^2);
        averageVelocity(idx) = mean(velocity{idx});
        stdVelocity(idx) = std(velocity{idx});
    for idxAv = 1:length(track)-nbPts  
        MovAverageVelocity{idx}(idxAv) = mean(velocity{idx}(idxAv:idxAv+nbPts-1));
        MovStdVelocity{idx}(idxAv) = std(velocity{idx}(idxAv:idxAv+nbPts-1));
        MovMaxVelocity{idx}(idxAv) = max(velocity{idx}(idxAv:idxAv+nbPts-1));
    end
     iTrack = [];
     track = [];
end
    