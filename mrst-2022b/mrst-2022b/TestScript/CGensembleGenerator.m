%load template files and get variable

fname='caseTemplate.txt';
%!!important note!! make sure array dimension for pv, perm consisten with
%active cell num in coarse grid


file = fullfile(getcurrentdir,fname);
fileID = fopen(file);

deck=textscan(fileID, '%s');
deck=deck{:};
fclose(fileID);

%% #############------------ ADD PARAMETER ------------##################
%nParams = 9;

ind_tTrain=find(contains(deck,'tTrain'));
tTrain = deck((ind_tTrain(1)+1));
tTrain = cellfun(@(a)str2double(a), tTrain);

ind_tRun=find(contains(deck,'tRun'));
tRun = deck((ind_tRun(1)+1));
tRun = cellfun(@(a)str2double(a), tRun);

ind_grid=find(contains(deck,'grid'));
ind_nw=find(contains(deck,'nw'));
ind_ng=find(contains(deck,'ng'));
ind_kw=find(contains(deck,'kw'));
ind_kg=find(contains(deck,'kg'));
ind_srw=find(contains(deck,'srw'));
ind_wi=find(contains(deck,'wi'));
ind_pv=find(contains(deck,'pv'));
ind_permIJ=find(contains(deck,'permIJ'));
ind_permK=find(contains(deck,'permK'));

gsize=deck((ind_grid(1)+1):(ind_grid(2)-1));
gsize=cellfun(@(a)str2double(a), gsize);

nw=deck((ind_nw(1)+1):(ind_nw(2)-1));
nw=cellfun(@(a)str2double(a), nw);

ng=deck((ind_ng(1)+1):(ind_ng(2)-1));
ng=cellfun(@(a)str2double(a), ng);

kw=deck((ind_kw(1)+1):(ind_kw(2)-1));
kw=cellfun(@(a)str2double(a), kw);

kg=deck((ind_kg(1)+1):(ind_kg(2)-1));
kg=cellfun(@(a)str2double(a), kg);

srw=deck((ind_srw(1)+1):(ind_srw(2)-1));
srw=cellfun(@(a)str2double(a), srw);

wi=deck((ind_wi(1)+1):(ind_wi(2)-1));
wi=cellfun(@(a)str2double(a), wi);

pv=deck((ind_pv(1)+1):(ind_pv(2)-1));
pv=cellfun(@(a)str2double(a), pv);

permIJ=deck((ind_permIJ(1)+1):(ind_permIJ(2)-1));
permIJ=cellfun(@(a)str2double(a), permIJ);

permK=deck((ind_permK(1)+1):(ind_permK(2)-1));
permK=cellfun(@(a)str2double(a), permK);

%% Modify Parameters 

%% Write Data
ensSize=50;

nmin=1; nmax=6; 
kmin=0.1; kmax=1;

for i=1:ensSize

    %change params here
    nw=rand(4,1)*(nmax-nmin)+nmin;
    ng=rand(4,1)*(nmax-nmin)+nmin;
    kw=rand(4,1)*(kmax-kmin)+kmin;
    kg=rand(4,1)*(kmax-kmin)+kmin;

    ensname=['ens-',num2str(i),'.txt'];
    outID=fopen(ensname,'w');

    fprintf(outID,'%s\n','tTrain');
    fprintf(outID,'%d\n',tTrain);
    fprintf(outID,'%s\n','\tTrain');
    fprintf(outID,'%s\n','tRun');
    fprintf(outID,'%d\n',tRun);
    fprintf(outID,'%s\n','\tRun');
    fprintf(outID,'%s\n','grid');
    fprintf(outID,'%d\n',gsize);
    fprintf(outID,'%s\n','\grid');
    fprintf(outID,'%s\n','nw');
    fprintf(outID,'%d\n',nw);
    fprintf(outID,'%s\n','\nw');
    fprintf(outID,'%s\n','ng');
    fprintf(outID,'%d\n',ng);
    fprintf(outID,'%s\n','\ng');
    fprintf(outID,'%s\n','kw');
    fprintf(outID,'%d\n',kw);
    fprintf(outID,'%s\n','\kw');
    fprintf(outID,'%s\n','kg');
    fprintf(outID,'%d\n',kg);
    fprintf(outID,'%s\n','\kg');
    fprintf(outID,'%s\n','srw');
    fprintf(outID,'%d\n',srw);
    fprintf(outID,'%s\n','\srw');
    fprintf(outID,'%s\n','wi');
    fprintf(outID,'%d\n',wi);
    fprintf(outID,'%s\n','\wi');
    fprintf(outID,'%s\n','pv');
    fprintf(outID,'%d\n',pv);
    fprintf(outID,'%s\n','\pv');
    fprintf(outID,'%s\n','permIJ');
    fprintf(outID,'%d\n',permIJ);
    fprintf(outID,'%s\n','\permIJ');
    fprintf(outID,'%s\n','permK');
    fprintf(outID,'%d\n',permK);
    fprintf(outID,'%s\n','\permK ');


    fclose(outID);


end

    function currentDir=getcurrentdir()
        if isdeployed %stand-alone mode
            %[status, result] = system('path');
            [~, result] = system('path');
            currentDir = char(regexpi(result, 'Path=(.*?);','tokens','once'));
        else
            currentDir = pwd;
        end
    end
