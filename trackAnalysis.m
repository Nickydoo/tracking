%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title : trackAnalysis
% description : takes tracks from SegmentOnline and analyses them
% author : Nicolas Desjardins-Lecavalier
% last edition : 2019-02-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

addpath('C:\work\nicolas\segementation-de-cellules-master\kakearney-boundedline-pkg-50f7e4b\Inpaint_nans')
addpath('C:\work\nicolas\segementation-de-cellules-master\kakearney-boundedline-pkg-50f7e4b\boundedline')
addpath('C:\work\nicolas\segementation-de-cellules-master\kakearney-boundedline-pkg-50f7e4b\catuneven')
addpath('C:\work\nicolas\segementation-de-cellules-master\kakearney-boundedline-pkg-50f7e4b\singlepatch')
addpath('C:\work\nicolas\segementation-de-cellules-master\kakearney-boundedline-pkg-50f7e4b\boundedline')

imagesFolder = fullfile('G:' , 'Nicolas' , '20190127scan35mm');
numPositions = 12;
tmin = [0.5 1 2 3 6 8 12]; %durée d'une track minimale (h)
dt = 2/60; %durée entre deux acquisition (h)
ntvect = tmin/dt;
period = 8;
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

% pour la track totale
[speed,vx,vy,averageSpeed,stdSpeed,lengthTrack,d,dtot,dnet,dmax,MSD,MI,OR,vAC,alpha,phi,DA] = totalTrackProperties(tracks, N2, dt);
selectAverageSpeed = find( averageSpeed >(mean(averageSpeed)+std(averageSpeed)));
selectMSD = find( MSD >(mean(MSD)+std(MSD)));
selectDtot = find( dtot >(mean(dtot)+std(dtot)));
selectDnet = find( dnet >(mean(dnet)+std(dnet)));
selectDmax = find( dmax >(mean(dmax)+std(dmax)));
selectMI = find( MI >(mean(MI)+std(MI)));
selectOR = find( OR >(mean(OR)+std(OR)));
%%
periodmov = 2:2:8;
for iPeriod = 1:length(periodmov)


[MovAverageVelocity,MovStdVelocity,MovMaxVelocity,velocity,vx,vy,averageVelocity,stdVelocity,lengthTrack] = movingAverageVelocity(tracks, N2, periodmov(iPeriod), dt);
[MOVaverageSpeed,MOVstdSpeed,MOVdtot,MOVdnet,MOVdmax,MOVMSD,MOVMI,MOVOR] = partialTrackProperties(tracks, N2, dt, periodmov(iPeriod));
[movR12,movR22,movRG2,movang,mova2,movA2,R12,R22,RG2,ang,a2,A2] = movingGyrationTensor(tracks, N2, periodmov(iPeriod), dt);

for idx = 1:sum(N2);
    conservationA2(idx) = mean(diff(movA2{idx})).^2;
    conservationV(idx) = mean(diff(MovAverageVelocity{idx})/max(max([MovAverageVelocity{:}]))).^2;
end


%% évaluation du tracking

%% figures

%analyse de la longueur des tracks
% lnplot = 2;
% for it2 = 1:length(ntvect)
%     nameLeg{it2} = ['plus longues que' num2str(tmin(it2)) 'h'];
% end

% figure
% subplot(1,2,1)
% bar([N1',NperH])
% xlabel('position')
% ylabel('nombre de tracks')
% legend(['initial',nameLeg])
% subplot(1,2,2) 
% bar(NperH./N1')
% xlabel('position')
% ylabel('proportion de tracks conservées')
% legend(nameLeg)

% figure;
% subplot(2,2,1)
% plot(averageVelocity,stdVelocity,'o')
% xlabel('average Velocity')
%  ylabel('std')
%  
% subplot(2,2,2)
% plot(lengthTrack,averageVelocity,'o')
% ylabel('<v>')
%  xlabel('length (h)')
% 
%  subplot(2,2,3)
% [avePlot,iav] = sort(averageVelocity);
%  errorbar(1:sum(N2),avePlot,stdVelocity(iav),'o')
%  ylabel('average Velocity + std')
%  xlabel('movie frame')
%  
%  subplot(2,2,4)
%  plot(averageVelocity,A2,'o')
%  xlabel('<v>')
%  ylabel('A2')
%  
%  title(['periodmov = ' num2str(periodmov(iPeriod))])
 
%  
%  figure
%   for iplot = 1:sum(N2)
%     boundedline(1:length(MovAverageVelocity{iplot}),MovAverageVelocity{iplot}, MovStdVelocity{iplot},'alpha','cmap',rand(1,3))
%   end
%  
%  ylabel('Ave vel 2h')


figure
 for iplot = 1:sum(N2)
     
     subplot(1,4,1)
     hold on
     plot(MovAverageVelocity{iplot})
     ylabel('<v> 2h')
     title(['periodmov = ' num2str(periodmov(iPeriod))])
     legend(cellfun(@num2str,mat2cell(1:sum(N2),1,ones(1,sum(N2))),'UniformOutput',0))
     
     subplot(1,4,2)
     hold on
     plot(movA2{iplot})
     ylabel('A2 2h')
     title(['periodmov = ' num2str(periodmov(iPeriod))])
     
    subplot(1,4,3)
    hold on
    plot(MovAverageVelocity{iplot},movA2{iplot},'o')
    xlabel('<v> 2h')
    ylabel('A2 2h')
    title(['periodmov = ' num2str(periodmov(iPeriod))])
    
    subplot(1,4,4)
    hold on
    plot(conservationA2,'o')
    plot(conservationV,'o')
    legend('conserv A2','conserv V')
    title(['periodmov = ' num2str(periodmov(iPeriod))])
 end
 


clear  MovAverageVelocity MovStdVelocity MovMaxVelocity velocity vx vy averageVelocity stdVelocity lengthTrack conservationA2  conservationV
end
%figure des moyennes
