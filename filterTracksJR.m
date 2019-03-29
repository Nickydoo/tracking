% Get Subset of tracks at least 20 frames long

function outTracks = filterTracksJR(inTracks, minTrackLength,start,jump,smoothFactor)

%minTrackLength = 20;

outTracks = [];
if isempty(inTracks), return, end

% take few frames
pos = logical(sum(inTracks(:,3) == start:jump:length(inTracks(:,3)),2));
inTracks = inTracks(pos,:);

% Get track IDs and their number of frames
[id,idCnt,~] = countEntries(inTracks(:,4),0, 0);

% Select IDs
goodIDs = id(idCnt >= minTrackLength);

% Get tracks with selected IDs
outTracks  = inTracks(ismember(inTracks(:,4),goodIDs),:);

% smooth tracks
for idx = goodIDs'
    posi = outTracks(:,4)==idx;
    outTracks(posi,1)  = smooth(outTracks(posi,1),smoothFactor);
    outTracks(posi,2)  = smooth(outTracks(posi,2),smoothFactor);
end

end