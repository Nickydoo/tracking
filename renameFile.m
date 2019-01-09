function outName = renameFile(inName)
    
    position = regexp(inName,'(?<=_s)(.[0-9]?)(?=\_t)','match');
    myTime   = regexp(inName,'(?<=_t)([0-9]+?)(?=\.)','match');
    
    if isempty(position) || isempty(myTime)
        error('Unable to extract position and time from Metamorph file name')
    end
    
    %outName = ['N' num2str(str2double(position)-1, '%03d') 'T' num2str(str2double(myTime), '%06d') '.jpg'];
    outName = ['N' num2str(str2double(position)-1, '%03d') 'T' num2str(str2double(myTime), '%06d') '.TIF']; % change made to keep high quality

end
