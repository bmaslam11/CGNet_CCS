fname='caseTest.txt';

load("pv.mat")
load("permIJ.mat")
load("permK.mat")

file = fullfile(getcurrentdir,fname);
fileID = fopen(file);

deck=textscan(fileID, '%s');
deck=deck{:};
fclose(fileID);

pv=num2cell(pv);
permIJ=num2cell(permIJ);
permK=permIJ;


deck=vertcat(deck,'pv', pv,'\pv','permIJ',permIJ,'\permIJ','permK',permK,'\permK');

writecell(deck,'caseInput.txt','Delimiter','space')

function currentDir=getcurrentdir()
    if isdeployed %stand-alone mode
        [status, result] = system('path');
        currentDir = char(regexpi(result, 'Path=(.*?);','tokens','once'));
    else
        currentDir = pwd;
    end
end