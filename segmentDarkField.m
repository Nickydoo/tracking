function [imOut, cellStats]=segmentDarkField(imInMt, minObjectSize)
    I = adapthisteq(mat2gray(imInMt));
    mask = imbinarize(I);
    mask = imfill(mask,'holes');
    imOut=bwareaopen(mask, minObjectSize);
    cellStats=struct2table(regionprops(imOut, 'Area', 'Centroid'));    
end
    