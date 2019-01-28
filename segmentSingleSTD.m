function [imOut, cellStats]=segmentSingleSTD(imInMt, imInAd, minObjectSize, regionMask)

    %imageSTD_Mt  = stdfilt(imInMt, ones(21)); 
    %%maskSTD_Mt   = im2bw(imageSTD_Mt, graythresh(imageSTD_Mt)); 
    %maskClean    = bwareaopen(maskSTD_Mt, (2*21)^2);
%     maskClean    = imclose(maskClean, strel('disk', 5));
    %maskClean    = maskClean & logical(regionMask);
    %cellStats_Mt = regionprops('table', maskClean, 'Area', 'PixelIdxList');
    
    %imInAd = adapthisteq(imInAd,'cliplimit',0.0001,'NumTiles',[4 4]);%,'Distribution','exponential');
    imageSTD_Ad  = stdfilt(imInAd, strel('disk',2,0).Neighborhood); 
%     imageSTD_Ad  = stdfilt(imInAd, strel('disk',10,0).Neighborhood); 
%     maskSTD_Ad   = imbinarize(mat2gray(imageSTD_Ad),graythresh(mat2gray(imageSTD_Ad))*2); 
    maskSTD_Ad   = imbinarize(mat2gray(imageSTD_Ad)); 
    maskClean = imfill(maskSTD_Ad,'holes');
    %for itCell=1:size(cellStats_Mt, 1)
        %if cellStats_Mt.Area(itCell) > 30*minObjectSize % >5*mean(cellStats_Mt.Area) 
            %maskClean(cellStats_Mt.PixelIdxList{itCell})=maskSTD_Ad(cellStats_Mt.PixelIdxList{itCell});
        %end
    %end
    imOut=bwareaopen(maskClean, minObjectSize);
    
    cellStats_Cl= regionprops('table', imOut, 'Area', 'Centroid', 'PixelIdxList', 'Solidity');
    
%     for itCell=1:size(cellStats_Cl, 1)
%         if cellStats_Cl.Area(itCell)>5*mean(cellStats_Cl.Area)
%             maskClean(cellStats_Cl.PixelIdxList{itCell})=0;
%         end
%     end
    
    for itCell=1:size(cellStats_Cl, 1)
        if cellStats_Cl.Solidity(itCell) < 0.85 && cellStats_Cl.Area(itCell)<5*mean(cellStats_Cl.Area) && cellStats_Cl.Area(itCell) > 20*minObjectSize
            gluedCells=zeros(size(imInAd));
            gluedCells(cellStats_Cl.PixelIdxList{itCell})=1;
            myDist=-(bwdist(~gluedCells));
            k = round(cellStats_Cl.Area(itCell)/10/minObjectSize);
            if k > 1 && k <4
                cellKMeansCentroids=getCellCentroidsKmeans(gluedCells, imInAd,k);
                centroidsImage=zeros(size(imInAd));
                centroidsImage(sub2ind(size(imInAd), round(cellKMeansCentroids(:,1)), round(cellKMeansCentroids(:,2))))=1;
                myDistExtMinima=imimposemin(myDist, imdilate(centroidsImage, strel('disk', 5)));
                wshTheseCells = watershed(myDistExtMinima);
                %wshTheseCells = watershed(myDist);
                gluedCells(wshTheseCells==0)=0;
            end
            maskClean(cellStats_Cl.PixelIdxList{itCell})=gluedCells(cellStats_Cl.PixelIdxList{itCell});
        end
    end
    
    imOut=bwareaopen(maskClean, minObjectSize);
    cellStats=struct2table(regionprops(imOut, 'Area', 'Centroid'));
    
end
    