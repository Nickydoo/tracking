%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title : trackAnalysis
% description : takes tracks from SegmentOnline and analyses them
% author : Nicolas Desjardins-Lecavalier
% last edition : 2019-02-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

addpath('C:\work\nicolas\kakearney-boundedline-pkg-50f7e4b\Inpaint_nans')
addpath('C:\work\nicolas\kakearney-boundedline-pkg-50f7e4b\boundedline')
addpath('C:\work\nicolas\kakearney-boundedline-pkg-50f7e4b\catuneven')
addpath('C:\work\nicolas\kakearney-boundedline-pkg-50f7e4b\singlepatch')
addpath('C:\work\nicolas\kakearney-boundedline-pkg-50f7e4b\boundedline')

imagesFolder = fullfile('G:' , 'Nicolas' , '20190127scan35mm');
numPositions = 12;
tmin = [0.5 1 2 3 6 8 12]; %durée d'une track minimale (h)
dt = 2/60; %durée entre deux acquisition (h)
ntvect = tmin/dt;
period = 12;
nt = period/dt;
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
        
        for it = 1:length(ntvect)
            selectPos_tmp = L{idx}>ntvect(it); %tracks plus longue
            NperH(idx,it) = sum(selectPos_tmp); %nombre de tracks plus longues que nt
        end
        
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

[MovAverageVelocity,MovStdVelocity,MovMaxVelocity,velocity,vx,vy,averageVelocity,stdVelocity,lengthTrack] = movingAverageVelocity(tracks, N2, 1, dt);

%A2, paramètres de marche aléatoire
%pour des tranches de 2h :) et totale évidemment

%% évaluation du tracking

%% figures

%analyse de la longueur des tracks
lnplot = 2;
for it2 = 1:length(ntvect)
    nameLeg{it2} = ['plus longues que' num2str(tmin(it2)) 'h'];
end

figure
subplot(1,2,1)
bar([N1',NperH])
xlabel('position')
ylabel('nombre de tracks')
legend(['initial',nameLeg])
subplot(1,2,2) 
bar(NperH./N1')
xlabel('position')
ylabel('proportion de tracks conservées')
legend(nameLeg)

figure;
plot(averageVelocity,stdVelocity,'o')
xlabel('average Velocity')
 ylabel('std')
 
 figure;
plot(lengthTrack,averageVelocity,'o')
ylabel('average Velocity')
 xlabel('length (h)')

[avePlot,iav] = sort(averageVelocity);
 figure;
 errorbar(1:sum(N2),avePlot,stdVelocity(iav),'o')
 ylabel('average Velocity + std')
 xlabel('movie frame')
 
 
 
 figure
  for iplot = 1:sum(N2)
    boundedline(1:length(MovAverageVelocity{iplot}),MovAverageVelocity{iplot}, MovStdVelocity{iplot},'alpha','cmap',rand(1,3))
  end
 
 ylabel('Ave vel 2h')

figure
 for iplot = 1:sum(N2)
     hold on
     plot(MovAverageVelocity{iplot})
 end
ylabel('max vel 2h')
%figure des moyennes
