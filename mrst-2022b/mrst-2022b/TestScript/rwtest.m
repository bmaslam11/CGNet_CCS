%function rwtest(fname)
fname='caseTest.txt';

file = fullfile(getcurrentdir,fname);
fileID = fopen(file);

deck=textscan(fileID, '%s');
deck=deck{:};
fclose(fileID);

nParams = 2;

ind_grid=find(contains(deck,'grid'));
ind_nw=find(contains(deck,'nw'));


gsz=vertcat(deck{ind_grid(1)+1:ind_grid(2)-1});
gsz=str2num(gsz);

nw=vertcat(deck{ind_nw(1)+1:ind_nw(2)-1});
nw=str2num(nw);

data=[gsz*10;nw*10];


casename=erase(fname,".txt");
outID=fopen(['out-',casename,'.txt'],'w');
fprintf(outID,'%d\n',data);
fclose(outID);


function currentDir=getcurrentdir()
    if isdeployed %stand-alone mode
        [status, result] = system('path');
        currentDir = char(regexpi(result, 'Path=(.*?);','tokens','once'));
    else
        currentDir = pwd;
    end
end




%end
