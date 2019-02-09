%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title : trackAnalysis
% description : takes tracks from SegmentOnline and analyses them
% author : Nicolas Desjardins-Lecavalier
% last edition : 2019-02-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imagesFolder = fullfile('G:' , 'Nicolas' , '20190127scan35mm');
numPositions = 100;
tmin = 6; %durée d'une track minimale (h)
dt = 2/60; %durée entre deux acquisition (h)
nt = tmin/dt;
%theseFileNames = fullfile(rawDir,{theseFileNames(:).name});

for idx = 1:numPositions
    resultsFolder{idx} = fullfile(imagesFolder,['Results' num2str(idx-1)],'trackResult.mat');
    if ~exist(resultsFolder{idx})
        tracks_tmp = [];
        tracks{idx} = tracks_tmp;
    else
        tracks_tmp = load(resultsFolder{idx});
        tracks{idx} = tracks_tmp.tracksuT;
        N1(idx) = tracks{idx}(end,end); %nombre de tracks initiale
        list1 = 1:N1(idx);
        L{idx} = sum(tracks{idx}(:,end) == list1,1)'; %longueur des tracks
        selectPos1 = L{idx}>nt; %tracks plus longue
        N2(idx) = sum(selectPos1); %nombre de tracks plus longues que 2h
        selectedTrackNb{idx} = list1(selectPos1); %numéro des tracks conservées
        selectPos2 = sum(tracks{idx}(:,end) == selectedTrackNb{idx},2); %positions des tracks conservées
        tracks{idx}(~selectPos2,:) = []; %supression des tracks trop courtes
        
        %renumérotation des tracks
        for idx2 = 1:length(selectedTrackNb{idx})
            pos2Chn = tracks{idx}(:,end)==selectedTrackNb{idx}(idx2);
            tracks{idx}(pos2Chn,end)=idx2;
        end
    end
end
%% analyse des tracks

%A2, vitesse moyenne, vitesse max, vitesse vectorielle, paramètres de
%marche aléatoire

%pour des tranches de 2h :) et totale évidemment

%% évaluation du tracking

%% figures
figure
subplot(1,2,1)
bar([N1',N2'])
xlabel('position')
ylabel('nombre de tracks')
legend('initial','plus longues que 2h')
subplot(1,2,2)
bar(N2./N1)
xlabel('position')
ylabel('proportion de tracks conservées')
legend('initial','plus longues que 2h')