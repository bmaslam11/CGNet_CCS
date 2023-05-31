clear; clc; close all;
%% CGnetESMDA
casename ='test_randGaussPolar_oneInj_preK_Mgridb';
obs_dir  ='CGnet_obsdata_oneInj_preloadK';


rnd_seed=9005;
rng(rnd_seed);

%period of training and forecast (years)
t_train = 30;
t_train_plume = 8;
t_run   = 50;

%get basemodel
load("CGnet_obsdata_oneInj_preloadK\base_fine\model.mat")
load("CGnet_obsdata_oneInj_preloadK\base_fine\q.mat")
load("CGnet_obsdata_oneInj_preloadK\base_fine\all_G.mat")
load("CGnet_obsdata_oneInj_preloadK\base_fine\all_rock.mat")
load("CGnet_obsdata_oneInj_preloadK\base_fine\W_fine.mat")

load("CGnet_obsdata_oneInj_preloadK\base_cg\cModel.mat")
load("CGnet_obsdata_oneInj_preloadK\base_cg\cPredSched.mat")
load("CGnet_obsdata_oneInj_preloadK\base_cg\cState0.mat")

load("PERM\PERMALL_KL_n_100_100_11_101.mat"); %PERMALL
W_fine=W;
%% Check Reference Model
conv_md=0.9869233*10^-16;
%fine model
figure; plotToolbar(all_G, all_rock.perm(:,1)/conv_md); c=colorbar; colormap jet% perm
clim([min(all_rock.perm(:,1)/conv_md) max(all_rock.perm(:,1)/conv_md)])
%CG model
figure; plotToolbar(cModel.G, cModel.rock.perm(:,1)/conv_md); c=colorbar; colormap jet% perm
clim([min(all_rock.perm(:,1)/conv_md) max(all_rock.perm(:,1)/conv_md)])

%% History Match Setup

cModelOri=cModel;

W=cPredSched.control.W;
n_well = length(W);
n_well_prd = get_nWellPrd(W);
n_cell = cModelOri.G.cells.num;

opt=[];
hm =[];

opt.NI = 15;
opt.NJ = 15;
opt.NK = 1;

%paramToCalibrate
opt.permIJ.flag  = 1;
opt.permK.flag   = 0;
opt.pv.flag      = 1;
opt.wi.flag      = 1;
opt.nw.flag      = 1;
opt.ng.flag      = 1;
opt.kw.flag      = 1;
opt.kg.flag      = 1;
opt.srw.flag     = 0;

nParams = get_nParam(cModelOri,W,opt);

%params range
opt.permIJ.min = 1*min(cModelOri.rock.perm(:,1));
opt.permIJ.max = 1*max(cModelOri.rock.perm(:,1));
opt.permK.min = 0.1*min(cModelOri.rock.perm(:,3));
opt.permK.max = 10*max(cModelOri.rock.perm(:,3));
opt.pv.min = 1*min(cModelOri.operators.pv);
opt.pv.max = 1*max(cModelOri.operators.pv);

wi=vertcat(W.WI);
opt.wi.size=size(wi);
opt.wi.min = min(wi);
opt.wi.max = max(wi);


opt.nw.min = 1.5; 
opt.nw.max = 2;

opt.ng.min = 1.5; 
opt.ng.max = 2;

opt.kw.min = 0.7; 
opt.kw.max = 1;

opt.kg.min = 0.7; 
opt.kg.max = 1;

opt.srw.min = 0.1; 
opt.srw.max = 0.4;

%obs data to match
hm.qbhp.flag    = 1;
hm.plum.flag    = 0;
hm.plum.freq    = 20; %plume boundary 
hm.sgas.flag    = 0;

hm.nWell        = 1; %I1 I2

%Get Obs Data
[y_obs,col,T_obs]=get_obsData(t_train,t_train_plume,hm); %put text data in folder spec by get_obsData

hm.qbhp.idx=col;
hm.qbhp.size=[t_train nnz(col)];

angSample=15;
xq=-180:angSample:180;
hm.plum.size=length(xq)*hm.nWell*t_train_plume;



%% ESMDA Parameters

nEnsemble   = 100;
Na          = 4; %data assimilation iteration

alpha=[28/3; 7; 4; 2];

%alpha   = (1/Na)*ones(Na,1);
% alpha=[100; 75; 50; 20; 1.103];
%alpha=[6; 6; 6; 6];

obs_dim = length(y_obs);
obs_sigma = getSigma(t_train,t_train_plume,hm);
obs_var = obs_sigma.^2;
Cd = diag(obs_var);
Cds = sparse(Cd);
%spy(Cds)

y_obs=repmat(y_obs,1,nEnsemble);

%%
%get prior model parameter
%[m_prior,opt]=initializeEnsemble(cModelOri,nEnsemble,opt,casename);
[m_prior,opt]=initializeEnsemblePreloadK(model,q,cModelOri,nEnsemble,opt,casename,all_G,all_rock,PERMALL);
m_prior=cell2mat(m_prior);
m_prior0=m_prior;



%% test scaling

% m_priorT=transpose(m_prior0);
% m_prior_scaled = scaleMdl(transpose(m_prior0),opt);
% m_prior_unscaled = unscaleMdl(m_prior_scaled,opt);
% histogram(m_prior_scaled((end-7):end,:))
% delta=m_priorT-m_prior_unscaled;
%% run prior ensemble-test
prior_ens = ['prior_',casename];

priordir=['priorTest_',casename,'_result'];

tic
for i = 1:nEnsemble %can be updated to parfor
    CGsimEnsembleFwd(i,m_prior,opt,priordir,t_train,cModelOri,W_fine,cState0)
end
toc

%% Inspect Ensemble

ens_id = 11;
CGsimEnsembleFwdChk(ens_id,m_prior,opt,priordir,t_train,cModelOri,W_fine,cState0)
%%
S_thresh = 0.1;
nPlume = 10; %all plume that does not merge

for i = 1:nEnsemble
    writePlume(i,priordir,model,cModel,W,S_thresh,hm.nWell,t_train_plume)
end

%% get y_prior
ens_list=1:nEnsemble;
y_prior=get_yEns(priordir,t_train,t_train_plume,nEnsemble,ens_list,hm);
y_0=y_prior;
%% Check Prior Ensemble Result 
k=8;
plotPlumeProb(hm,obs_dir,priordir,t_train_plume,nEnsemble,k)
%ylim([0 14000])
%%
plotEnsemble(col,T_obs,priordir,t_train,t_train,nEnsemble)
sgtitle('Prior Dist','FontWeight','bold')
%% ESMDA Iteration
clc;
tic
%error matrix data
error_en=zeros(Na,nEnsemble);
f = waitbar(0,'Starting ESMDA iteration');


for i = 1:Na
    
    waitbar(i/Na,f,['ESMDA iteration # ',num2str(i),' of ',num2str(Na)])
    iterdir=['esmda_',casename,'_iter_',num2str(i)];
    %evaluate prior params
    
    status = ones(nEnsemble,Na);
    status = logical(status);
    msg = cell(nEnsemble,Na);

    parfor iter = 1:nEnsemble %can be updated to parfor
        try
            j=ens_list(iter);
            CGsimEnsembleFwd(j,m_prior,opt,iterdir,t_train,cModelOri,W,cState0)
            status(iter,i)=true;
            msg{iter,i}=[];
        catch ME
            status(iter,i)=false;
            msg{iter,i}=ME;
        end
    end

    nEnsemble = nnz(status(:,i)==1);   %surviving ensemble count
    ens_list=ens_list(status(:,i)==1); %surviving ensemble ID

    if hm.plum.flag
        for id = 1:nEnsemble
            writePlume(ens_list(id),iterdir,model,cModel,W,S_thresh,hm.nWell,t_train)
        end
    end

    y_prior=get_yEns(iterdir,t_train,t_train_plume,nEnsemble,ens_list,hm);
    
    %filter y_obs
    y_obs=y_obs(:,ens_list);
    
    %add noise to d_obs
    y_uc=y_obs+sqrt(alpha(i)).*(sqrt(Cd)*randn(size(y_obs)));
    
    %filter m_prior    
    m_prior=m_prior(ens_list,:);
    m_prior=transpose(m_prior);

    %scale m_prior
    m_prior=scaleMdl(m_prior,opt);
    %get MSE
    for k=1:nEnsemble
        
        error_en(i,k)=norm(y_prior(:,k)-y_obs(:,k))/norm(y_obs(:,1));
    end

    mu_m = mean(m_prior,2); %average params of all ensemble
    mu_d = mean(y_prior,2);

    Cmd = mtimes(m_prior-mu_m,transpose(y_prior-mu_d))/(nEnsemble-1);
    Cdd = mtimes(y_prior-mu_d,transpose(y_prior-mu_d))/(nEnsemble-1);

    %Kalman Gain calculation by pseudo-inverse
    K = mtimes(Cmd,pinv(Cdd + alpha(i)*Cd));
    m_pred=m_prior+mtimes(K,y_uc-y_prior);
    
    %unscale m
    m_pred=unscaleMdl(m_pred,opt);
    %update m
    m_prior=transpose(m_pred);

end

toc
%%
figure
boxplot(transpose(error_en))
xlabel('ESMDA Iteration','FontWeight','bold')
ylabel('RMSE',FontWeight='bold')
%ylim([0 0.025])
%% run posterior ensemble-test prediction
prior_ens = ['post_',casename];

priordir=['postTest_',casename,'_result'];

status = ones(nEnsemble,1);
status = logical(status);
msg = cell(nEnsemble,1);


tic
parfor i = 1:nEnsemble %can be updated to parfor
    try
        CGsimEnsembleFwd(i,m_prior,opt,priordir,t_run,cModelOri,W,cState0)
        status(i)=true;
        msg{i}=[];
    catch ME
        status(i)=false;
        msg{i}=ME;
    end
end
toc
%%
for id = 1:nEnsemble
    writePlume(id,priordir,model,cModel,W,S_thresh,hm.nWell,t_train)
end
%% run prior ensemble-test prediction
prior_ens = ['prior_',casename];

priordir=['priorTest_',casename,'_result'];



tic
parfor i = 1:nEnsemble %can be updated to parfor
    CGsimEnsembleFwd(i,m_prior0,opt,priordir,t_run,cModelOri,W,cState0)
end
toc
%%
figure
i=4;
%iterdir=['esmda_',casename,'_iter_',num2str(i)];
 iterdir=['postTest_',casename,'_result'];
plotEnsemble(col,T_obs,iterdir,t_run, t_train,nEnsemble)
sgtitle(['ESMDA iteration #',num2str(i)],'FontWeight','bold')
%% Plot Plume
k=8;
plotPlumeProb(hm,obs_dir,iterdir,t_train_plume,nEnsemble,k)
ylim([0 14000])
%% PermIJ Hist Compare
%load("CGnet_basemodel\cModelOri_ens1.mat")
permIJ_prior_mean=mean(m_prior0(:,1:n_cell));
permIJ_post_mean=mean(m_prior(:,1:n_cell));
%%

convDarcy=9.869233*10^-13;
figure
hold on

h1=histogram(cModel.rock.perm(:,1)*1000/convDarcy);
h2=histogram(permIJ_prior_mean*1000/convDarcy);
h3=histogram(permIJ_post_mean*1000/convDarcy);
h1.Normalization = 'probability';
h1.BinWidth = 10;
h2.Normalization = 'probability';
h2.BinWidth = 10;
h3.Normalization = 'probability';
h3.BinWidth = 10;

xline(mean(permIJ_prior_mean*1000/convDarcy),'LineStyle','--','LineWidth',1.5,'Color','blue','Label',{'prior mean'})
xline(mean(permIJ_post_mean*1000/convDarcy),'LineStyle','--','LineWidth',1.5,'Color','red','Label',{'post mean'})
xline(mean(cModelOri.rock.perm(:,1)*1000/convDarcy),'LineStyle','-','LineWidth',1.5,'Label',{'truth mean'})

xlabel('Permeability (milliDarcy)','FontWeight','bold')
ylabel('Probability','FontWeight','bold')
box on
legend('truth','mean prior','mean post')
hold off


% figure
% hold on
% nbin=15;
% histogram(permIJ_prior_mean,nbin);
% histogram(permIJ_post_mean,nbin);
% histogram(cModelOri.rock.perm(:,1),nbin)

% legend('prior','posterior','truth')
% hold off
%%
permXref = cModelOri.rock.perm(:,1);
permXsim = permIJ_post_mean';
permXprior = permIJ_prior_mean';
%log transform


max_lim=0.9*max(permXref*1000/convDarcy);
min_lim=min(permXref*1000/convDarcy);

figure
tiledlayout(1,3)

nexttile
plotToolbar(cModelOri.G, permXref*1000/convDarcy);
colorbar;
colormap jet
caxis([min_lim max_lim])
title('Upscaled Reference Model (Kij)')
%view(-55, 60);

nexttile
plotToolbar(cModelOri.G, permXprior*1000/convDarcy);
%colorbar;
colormap jet
%caxis([min_lim max_lim])
title('Mean Prior Model (Kij)')

nexttile
plotToolbar(cModelOri.G, permXsim*1000/convDarcy);
%colorbar;
colormap jet
%caxis([min_lim max_lim])
title('Mean Posterior Model (Kij)')
%view(-55, 60);

%% CellPV Hist Compare
%load("CGnet_basemodel\cModelOri_ens1.mat")

pv_prior_mean=m_prior0(2,2*n_cell+1:3*n_cell);
pv_post_mean=mean(m_prior(:,2*n_cell+1:3*n_cell),1);

figure
hold on
nbin=15;
histogram(pv_prior_mean,nbin);
histogram(pv_post_mean,nbin);
histogram(cModelOri.operators.pv,nbin)
xline(mean(pv_prior_mean),'LineStyle','--','LineWidth',1.5,'Color','blue','Label',{'prior mean'})
xline(mean(pv_post_mean),'LineStyle','--','LineWidth',1.5,'Color','red','Label',{'post mean'})
xline(mean(cModelOri.operators.pv),'LineStyle','-','LineWidth',1.5,'Label',{'truth mean'})
legend('prior','posterior','truth')
hold off
%%
pvref = cModelOri.operators.pv;
pvsim = pv_post_mean';
pvprior = pv_prior_mean';
%log transform


max_lim=max([pvref;pvsim]);
min_lim=min([pvref;pvsim]);

figure
tiledlayout(1,3)

nexttile
plotToolbar(cModelOri.G, pvref);
colorbar;
colormap jet
% clim([min_lim max_lim])
title('Upscaled Reference Model (pv)')
%view(-55, 60);

nexttile
plotToolbar(cModelOri.G, pvprior);
colorbar;
colormap jet
% clim([min_lim max_lim])
title('Mean Prior Model (pv)')

nexttile
plotToolbar(cModelOri.G, pvsim);
colorbar;
colormap jet
% clim([min_lim max_lim])
title('Mean Posterior Model (pv)')
%view(-55, 60);

%% Helper Functions

function [y_obs,varargout]=get_obsData(t_train,t_train_plume,hm)
% [y_obs,col,T]=get_obsData(t_train,hm)
% col id to reshape array
% T is well prod table 

%read well data (in file)
T          =readtable("CGnet_obsdata_oneInj_preloadK\wellResult_Obs50_M.txt",'VariableNamingRule','preserve');
wellResult =table2array(T);
wellResult =wellResult(1:t_train,:);
varargout{2}=T;


wellObs=[];
if hm.qbhp.flag
    %select well obs data used in matching
    wellPropName = T.Properties.VariableNames;

    colPrdQw = contains(wellPropName,"qwS_P");
    colPrdQg = contains(wellPropName,"qgS_P");
    colPrdBHP = contains(wellPropName,"BHP_");

    col = colPrdQw+colPrdQg+colPrdBHP;
    col(col>1)=1;
    col=logical(col);

    wellObs = wellResult(:,col);
    wellObs = wellObs(:);
    varargout{1}=col;
end

sgasObs=[];
if hm.sgas.flag
    %read sat data
    S = readmatrix("CGnet_obsdata_oneInj_preloadK\SgasObs_.txt");
    S = transpose(S); %to get Sgas per cell ordered with time in row
    sgasObs = S(1:t_train,:);
    sgasObs = sgasObs(:);
end

splumeObs=[];
if hm.plum.flag
    %read sat data
    S = readmatrix("CGnet_obsdata_oneInj_preloadK\plumeData.txt");
    S = reshape(S,[],2,hm.nWell,hm.plum.freq);
    S = S(:,2,:,1:t_train_plume);

    splumeObs = S(:);
    
end

%append all obs data
y_obs=[wellObs;sgasObs;splumeObs];

assert(~isempty(y_obs),'Obs Data is not defined')

end

function y_ens=get_yEns(ens_dir,t_train,t_train_plume,nEnsemble,ens_list,hm)

col=hm.qbhp.idx;

wellResult=[];
if hm.qbhp.flag

    %preallocate array for ensemble obs data
    sz=[hm.qbhp.size nEnsemble];
    wellResult = zeros(sz);

    %read well data (in ensemble)
    for i=1:nEnsemble
        id=ens_list(i);
        data = fullfile(pwd,ens_dir,['wellResult_ens_',num2str(id),'.txt']);
        T_ens =readtable(data,'VariableNamingRule','preserve');
        T_ens =table2array(T_ens);
        T_ens =T_ens(1:t_train,col);
        wellResult(:,:,i)=T_ens;
    end
    wellResult=reshape(wellResult,[],nEnsemble);
end

sgasResult=[];
if hm.sgas.flag

    %preallocate array for ensemble obs data
    sz=[hm.sgas.size nEnsemble];
    sgasResult = zeros(sz);

    %read sat data (in ensemble)
    for i=1:nEnsemble
        id=ens_list(i);
        data = fullfile(pwd,ens_dir,['Sgas_ens-',num2str(id),'.txt']);
        S_ens =readmatrix(data);
        sgasResult(:,:,i)=S_ens;
    end
    sgasResult=reshape(sgasResult,[],nEnsemble); %need to be vectorized
    
    %modify for plume polar transform
    %...

end

plumeResult=[];
if hm.plum.flag

sz=[hm.plum.size nEnsemble];
plumeResult=zeros(sz);

for i=1:nEnsemble
    id=ens_list(i);
    data = fullfile(pwd,ens_dir,['plumeData_',num2str(id),'.txt']);
    P_ens= readmatrix(data);
    P_ens= reshape(P_ens,[],2,hm.nWell,t_train_plume);
    P_ens= P_ens(:,2,:,:);
    P_ens= P_ens(:);
    plumeResult(:,i)=P_ens;
end


end

%collect obs data
y_ens=[wellResult; sgasResult; plumeResult];

assert(~isempty(y_ens),'ensemble result is not defined');

end

function nParams=get_nParam(cModelOri,W,opt)

    nParams = 0;
    
    if opt.permIJ.flag
        nPermIJ = numel(cModelOri.rock.perm(:,1));
        nParams = nParams + nPermIJ;
    end
    
    if opt.permK.flag
        nPermK = numel(cModelOri.rock.perm(:,3));
        nParams = nParams + nPermK;
    end

    if opt.pv.flag
        nPV = numel(cModelOri.operators.pv);
        nParams = nParams + nPV;
    end

    if opt.wi.flag
        nWI = length(W);
        nParams = nParams + nWI;
    end

    if opt.nw.flag
        nnW = numel(cModelOri.fluid.nW);
        nParams = nParams + nnW;
    end

    if opt.ng.flag
        nnG = numel(cModelOri.fluid.nG);
        nParams = nParams + nnG;
    end

    if opt.kw.flag
        nkW = numel(cModelOri.fluid.kW);
        nParams = nParams + nkW;
    end

    if opt.kg.flag
        nkG = numel(cModelOri.fluid.kG);
        nParams = nParams + nkG;
    end

    if opt.srw.flag
        nsr = numel(cModelOri.fluid.srW);
        nParams = nParams + nsr;
    end

end

function n_well_prd=get_nWellPrd(W)

    wname = vertcat(W.name);
    n_well_prd = 0;
    for i =1:length(W)
        n_well_prd = n_well_prd + contains(wname(i),"P");
    end

end

function [params,opt]=initializeEnsemble(cModelOri,nEnsemble,opt,casename)
prior_ens = ['prior_',casename];

params = cell(nEnsemble,1);

if ~exist(prior_ens, 'dir')
    mkdir(prior_ens)
end

pvTot_ori = sum(cModelOri.operators.pv);

for i = 1: nEnsemble

    %write ensemble
    ensname=['ens-',num2str(i),'.txt'];
    outID=fopen([pwd,'\',prior_ens,'\',ensname],'w');

    %set initial ensemble with homogen properties 
    if opt.permIJ.flag
        permIJ=cModelOri.rock.perm(:,1);
        pert=0.1*mean(permIJ).*randn(size(permIJ));
        
        permIJ=permIJ+pert;
        fprintf(outID,'%s\n','permIJ');
        fprintf(outID,'%d\n',permIJ);
        fprintf(outID,'%s\n','\permIJ');
                    
        params{i}=horzcat(params{i},permIJ');
    end

    if opt.permK.flag

        permK=cModelOri.rock.perm(:,3);
        pert=0.1*mean(permK).*randn(size(permK));
        permK=permK+pert;
        
        fprintf(outID,'%s\n','permK');
        fprintf(outID,'%d\n',permK);
        fprintf(outID,'%s\n','\permK ');
        
        params{i}=horzcat(params{i},permK');
    end

    if opt.pv.flag 
  
        pv=cModelOri.operators.pv;
        pert=0.1*mean(pv).*randn(size(pv));
        pv=pv+pert;
        %pvTot_ens = sum(pv);
        
        %guarantee same total pv
        %pv=pv.*(pvTot_ori/pvTot_ens);

        fprintf(outID,'%s\n','pv');
        fprintf(outID,'%d\n',pv);
        fprintf(outID,'%s\n','\pv');
        
        params{i}=horzcat(params{i},pv');
    end

    if opt.wi.flag
        minval =opt.wi.min;
        maxval =opt.wi.max;
        sz_wi =opt.wi.size;
        %pert=0.1*mean(wi).*rand(sz_wi);
        %wi = wi+pert;
        wi = minval+(maxval-minval).*rand(sz_wi);
        fprintf(outID,'%s\n','wi');
        fprintf(outID,'%d\n',wi);
        fprintf(outID,'%s\n','\wi');

        params{i}=horzcat(params{i},wi');
    end

    if opt.nw.flag
        minval =opt.nw.min;
        maxval =opt.nw.max;
        nw=cModelOri.fluid.nW;
        nw=minval+(maxval-minval).*rand(size(nw));

        fprintf(outID,'%s\n','nw');
        fprintf(outID,'%d\n',nw);
        fprintf(outID,'%s\n','\nw');
                
        params{i}=horzcat(params{i},nw');
    end

    if opt.ng.flag
        minval =opt.ng.min;
        maxval =opt.ng.max;
        ng=cModelOri.fluid.nG;
        ng=minval+(maxval-minval).*rand(size(ng));

        fprintf(outID,'%s\n','ng');
        fprintf(outID,'%d\n',ng);
        fprintf(outID,'%s\n','\ng');

        params{i}=horzcat(params{i},ng');
    end

    if opt.kw.flag
        minval =opt.kw.min;
        maxval =opt.kw.max;
        kw=cModelOri.fluid.kW;
        kw=minval+(maxval-minval).*rand(size(kw));

        fprintf(outID,'%s\n','kw');
        fprintf(outID,'%d\n',kw);
        fprintf(outID,'%s\n','\kw');

        params{i}=horzcat(params{i},kw');
    end

    if opt.kg.flag
        minval =opt.kg.min;
        maxval =opt.kg.max;
        kg=cModelOri.fluid.kG;
        kg=minval+(maxval-minval).*rand(size(kg));

        fprintf(outID,'%s\n','kg');
        fprintf(outID,'%d\n',kg);
        fprintf(outID,'%s\n','\kg');

        params{i}=horzcat(params{i},kg');
    end

    if opt.srw.flag
        minval =opt.srw.min;
        maxval =opt.srw.max;
        srw=cModelOri.fluid.srW;
        srw=minval+(maxval-minval).*rand(size(srw));

        fprintf(outID,'%s\n','srw');
        fprintf(outID,'%d\n',srw);
        fprintf(outID,'%s\n','\srw');

        params{i}=horzcat(params{i},srw');
    end

end

if opt.permIJ.flag
    opt.permIJ.size=size(permIJ);
end

if opt.permK.flag
    opt.permK.size=size(permK);
end

if opt.pv.flag
    opt.pv.size=size(pv);
end

if opt.nw.flag
    opt.nw.size=size(nw);
end

if opt.ng.flag
    opt.ng.size=size(ng);
end

if opt.kw.flag
    opt.kw.size=size(kw);
end

if opt.kg.flag
    opt.kg.size=size(kg);
end

if opt.srw.flag
    opt.srw.size=size(srw);
end

end

function plotEnsemble(col, T_obs, ens_dir, t_train,t_hm,nEnsemble)

%tic
%ens_dir=['prior_',casename,'_result'];

conv_qw=-24*60*60*6.289;
conv_qg=-24*60*60*35.3/(1e6);
conv_p =1/6895; 


wellPropName = T_obs.Properties.VariableNames;
wellPropName = wellPropName(col); 
wellObs =table2array(T_obs);
wellObs =wellObs(:,col);

%preallocate array for ensemble obs data
wellResult = zeros(t_train,nnz(col),nEnsemble);

%read well data (in ensemble)
for i=1:nEnsemble
    data = fullfile(pwd,ens_dir,['wellResult_ens_',num2str(i),'.txt']);
    T_ens =readtable(data,'VariableNamingRule','preserve');
    T_ens =table2array(T_ens);
    T_ens =T_ens(1:t_train,col);
    wellResult(:,:,i)=T_ens;
end

t=1:length(wellObs);

%plot data
tiledlayout('flow')

%water rate well-1
nexttile
hold on

plot(t,wellObs(:,1)*conv_qw,'LineWidth',1.5,'Color','k')
plot(t(1:t_train),squeeze(wellResult(:,1,:))*conv_qw,'Color',[.7 .7 .7])
scatter(t(1:t_train),wellObs(1:t_train,1)*conv_qw,"blue","filled")
box on
ylim([0 3e5])

xline(t_hm,'LineStyle','--','LineWidth',1,'Color','r')

title(wellPropName{1},'Interpreter','none')
xlabel('time (years)')
ylabel('water rate (bwpd)')
hold off

%gas rate well-1
nexttile
hold on
plot(t,wellObs(:,4)*conv_qg,'LineWidth',1.5,'Color','k')
plot(t(1:t_train),squeeze(wellResult(:,4,:))*conv_qg,'Color',[.7 .7 .7])
scatter(t(1:t_train),wellObs(1:t_train,4)*conv_qg,"red","filled")
box on
ylim([0 1.5])

xline(t_hm,'LineStyle','--','LineWidth',1,'Color','r')

title(wellPropName{4},'Interpreter','none')
xlabel('time (years)')
ylabel('gas rate (MMscfd)')
hold off

%water rate well-2
nexttile
hold on
plot(t,wellObs(:,2)*conv_qw,'LineWidth',1.5,'Color','k')
plot(t(1:t_train),squeeze(wellResult(:,2,:))*conv_qw,'Color',[.7 .7 .7])
scatter(t(1:t_train),wellObs(1:t_train,2)*conv_qw,"blue","filled")
box on
ylim([0 3e5])

xline(t_hm,'LineStyle','--','LineWidth',1,'Color','r')

title(wellPropName{2},'Interpreter','none')
xlabel('time (years)')
ylabel('water rate (bwpd)')
hold off

%gas rate well-2
nexttile
hold on
plot(t,wellObs(:,5)*conv_qg,'LineWidth',1.5,'Color','k')
plot(t(1:t_train),squeeze(wellResult(:,5,:))*conv_qg,'Color',[.7 .7 .7])
scatter(t(1:t_train),wellObs(1:t_train,5)*conv_qg,"red","filled")
box on
ylim([0 1.5])

xline(t_hm,'LineStyle','--','LineWidth',1,'Color','r')

title(wellPropName{5},'Interpreter','none')
xlabel('time (years)')
ylabel('gas rate (MMscfd)')
hold off
%toc

%water rate well-3
nexttile
hold on
plot(t,wellObs(:,3)*conv_qw,'LineWidth',1.5,'Color','k')
plot(t(1:t_train),squeeze(wellResult(:,3,:))*conv_qw,'Color',[.7 .7 .7])
scatter(t(1:t_train),wellObs(1:t_train,3)*conv_qw,"blue","filled")
box on
ylim([0 3e5])

xline(t_hm,'LineStyle','--','LineWidth',1,'Color','r')

title(wellPropName{3},'Interpreter','none')
xlabel('time (years)')
ylabel('water rate (bwpd)')
hold off

%gas rate well-3
nexttile
hold on
plot(t,wellObs(:,6)*conv_qg,'LineWidth',1.5,'Color','k')
plot(t(1:t_train),squeeze(wellResult(:,6,:))*conv_qg,'Color',[.7 .7 .7])
scatter(t(1:t_train),wellObs(1:t_train,6)*conv_qg,"red","filled")
box on
ylim([0 1.5])

xline(t_hm,'LineStyle','--','LineWidth',1,'Color','r')

title(wellPropName{6},'Interpreter','none')
xlabel('time (years)')
ylabel('gas rate (MMscfd)')
hold off



%BHP well-I1
nexttile
hold on
plot(t,wellObs(:,7)*conv_p,'LineWidth',1.5,'Color','k')
plot(t(1:t_train),squeeze(wellResult(:,7,:))*conv_p,'Color',[.7 .7 .7])
scatter(t(1:t_train),wellObs(1:t_train,7)*conv_p,"black","filled")
box on
%ylim([4000 7000])

xline(t_hm,'LineStyle','--','LineWidth',1,'Color','r')

title(wellPropName{7},'Interpreter','none')
xlabel('time (years)')
ylabel('BHP (psi)')
hold off

%BHP well-M1
nexttile
hold on
plot(t,wellObs(:,8)*conv_p,'LineWidth',1.5,'Color','k')
plot(t(1:t_train),squeeze(wellResult(:,8,:))*conv_p,'Color',[.7 .7 .7])
scatter(t(1:t_train),wellObs(1:t_train,8)*conv_p,"black","filled")
box on
%ylim([4000 6000])

xline(t_hm,'LineStyle','--','LineWidth',1,'Color','r')

title(wellPropName{8},'Interpreter','none')
xlabel('time (years)')
ylabel('BHP (psi)')
hold off

%BHP well-M2
nexttile
hold on
plot(t,wellObs(:,9)*conv_p,'LineWidth',1.5,'Color','k')
plot(t(1:t_train),squeeze(wellResult(:,9,:))*conv_p,'Color',[.7 .7 .7])
scatter(t(1:t_train),wellObs(1:t_train,9)*conv_p,"black","filled")
box on

xline(t_train,'LineStyle','--','LineWidth',1,'Color','r')

title(wellPropName{9},'Interpreter','none')
xlabel('time (years)')
ylabel('BHP (psi)')
hold off


end

function [sigma_obs] = getSigma(t_train,t_train_plume,hm)

    obs_qbhp=[];
    if hm.qbhp.flag
        %read well data (in file)
        T          =readtable("CGnet_obsdata_oneInj_preloadK\wellResult_Obs50_M.txt",'VariableNamingRule','preserve');
        wellResult =table2array(T);
        wellResult =wellResult(1:t_train,:);
        
        wellPropName = T.Properties.VariableNames;
        
        %select well obs data used in matching
        colPrdQw = contains(wellPropName,"qwS_P");
        colPrdQg = contains(wellPropName,"qgS_P");
        colPrdBHP = contains(wellPropName,"BHP_");
        
        %add noise matrix
        qW=wellResult(:,colPrdQw);
        sigma_qW =max(abs(qW(:)))*0.05; % 5-10 percent of data
        sigma_qW =sigma_qW*ones(length(qW(:)),1);
        
        qG=wellResult(:,colPrdQg);
        sigma_qG =max(abs(qG(:)))*0.05; % 5-10 percent of data
        sigma_qG =sigma_qG*ones(length(qG(:)),1);
        
        bhp=wellResult(:,colPrdBHP);
        sigma_bhp =max(abs(bhp(:)))*0.05; % 5-10 percent of data
        sigma_bhp =sigma_bhp*ones(length(bhp(:)),1);

        obs_qbhp = [sigma_qW; sigma_qG; sigma_bhp];
    end
    
    obs_sgas=[];
    if hm.sgas.flag
        %read sat data
        S = readmatrix("CGnet_obsdata_oneInj_preloadK\SgasObs_.txt");
        S = transpose(S); %to get Sgas per cell ordered with time in row
        sgasObs = S(1:t_train,:);
        sgasObs = sgasObs(:);

        sigma_sG = max(abs(sgasObs))*0.05; % 5-10 percent of data
        sigma_sG = sigma_sG*ones(length(sgasObs),1);
        obs_sgas = sigma_sG;
    end
    
    obs_plum=[];
    if hm.plum.flag
        %read sat data
        S = readmatrix("CGnet_obsdata_oneInj_preloadK\plumeData.txt");
        S = reshape(S,[],2,hm.nWell,hm.plum.freq);
        S = S(:,2,:,1:t_train_plume);

        splumeObs = S(:);
        sigma_plume = max(abs(splumeObs))*0.07;
        sigma_plume = sigma_plume*ones(length(splumeObs),1);
        obs_plum = sigma_plume;

    end

    sigma_obs=[obs_qbhp; obs_sgas; obs_plum];
    


end

function m_prior_scaled = scaleMdl(m_prior,opt)

m_prior_scaled=[];

if opt.permIJ.flag
    minval =opt.permIJ.min;
    maxval =opt.permIJ.max;
    ind = 1:opt.permIJ.size(1);
    tmp=(m_prior(ind,:)-minval)./(maxval-minval);
    m_prior_scaled=[m_prior_scaled; tmp];

    m_prior(ind,:)=[];
end

if opt.permK.flag
    minval =opt.permK.min;
    maxval =opt.permK.max;
    ind = 1:opt.permK.size(1);
    tmp=(m_prior(ind,:)-minval)./(maxval-minval);
    m_prior_scaled=[m_prior_scaled; tmp];

    m_prior(ind,:)=[];
end

if opt.pv.flag
    minval =opt.pv.min;
    maxval =opt.pv.max;
    ind = 1:opt.pv.size(1);
    tmp=(m_prior(ind,:)-minval)./(maxval-minval);
    m_prior_scaled=[m_prior_scaled; tmp];

    m_prior(ind,:)=[];
end

if opt.wi.flag
    minval =opt.wi.min;
    maxval =opt.wi.max;
    ind = 1:opt.wi.size(1);
    tmp=(m_prior(ind,:)-minval)./(maxval-minval);
    m_prior_scaled=[m_prior_scaled; tmp];

    m_prior(ind,:)=[];
end

if opt.nw.flag
    minval =opt.nw.min;
    maxval =opt.nw.max;
    ind = 1:opt.nw.size(1);
    tmp=(m_prior(ind,:)-minval)./(maxval-minval);
    m_prior_scaled=[m_prior_scaled; tmp];

    m_prior(ind,:)=[];

end

if opt.ng.flag

    minval =opt.ng.min;
    maxval =opt.ng.max;
    ind = 1:opt.ng.size(1);
    tmp=(m_prior(ind,:)-minval)./(maxval-minval);
    m_prior_scaled=[m_prior_scaled; tmp];

    m_prior(ind,:)=[];
end

m_prior_scaled=[m_prior_scaled; m_prior];

end

function m_prior_unscaled = unscaleMdl(m_prior,opt)

m_prior_unscaled=[];

if opt.permIJ.flag
    minval =opt.permIJ.min;
    maxval =opt.permIJ.max;
    ind = 1:opt.permIJ.size(1);
    tmp=(m_prior(ind,:).*(maxval-minval))+minval;
    m_prior_unscaled=[m_prior_unscaled; tmp];

    m_prior(ind,:)=[];
end

if opt.permK.flag
    minval =opt.permK.min;
    maxval =opt.permK.max;
    ind = 1:opt.permK.size(1);
    tmp=(m_prior(ind,:).*(maxval-minval))+minval;
    m_prior_unscaled=[m_prior_unscaled; tmp];

    m_prior(ind,:)=[];
end

if opt.pv.flag
    minval =opt.pv.min;
    maxval =opt.pv.max;
    ind = 1:opt.pv.size(1);
    tmp=(m_prior(ind,:).*(maxval-minval))+minval;
    m_prior_unscaled=[m_prior_unscaled; tmp];

    m_prior(ind,:)=[];
end

if opt.wi.flag
    minval =opt.wi.min;
    maxval =opt.wi.max;
    ind = 1:opt.wi.size(1);
    tmp=(m_prior(ind,:).*(maxval-minval))+minval;
    m_prior_unscaled=[m_prior_unscaled; tmp];

    m_prior(ind,:)=[];
end

if opt.nw.flag
    minval =opt.nw.min;
    maxval =opt.nw.max;
    ind = 1:opt.nw.size(1);
    tmp=(m_prior(ind,:).*(maxval-minval))+minval;
    m_prior_unscaled=[m_prior_unscaled; tmp];

    m_prior(ind,:)=[];

end

if opt.ng.flag

    minval =opt.ng.min;
    maxval =opt.ng.max;
    ind = 1:opt.ng.size(1);
    tmp=(m_prior(ind,:).*(maxval-minval))+minval;
    m_prior_unscaled=[m_prior_unscaled; tmp];

    m_prior(ind,:)=[];
end

m_prior_unscaled=[m_prior_unscaled; m_prior];

end

function writePlume(ens_id,indir,model,cModel,W,S_thresh,nWell,nPlume)

data =fullfile(pwd,indir,['Sgas_ens_',num2str(ens_id),'.txt']);
S_ens=readmatrix(data);
S_ens_bin=imbinarize(S_ens,S_thresh);

%nPlume = size(S_ens,2);
%nPlume =10;
%remapping Sbin
S_bin_fine=nan([length(cModel.G.partition) nPlume]);

for i =1:nPlume
    S_bin_fine(:,i) = S_ens_bin(cModel.G.partition,i);
end

[ix ,jx, ~] = ind2sub(cModel.G.parent.cartDims, cModel.G.parent.cells.indexMap);

[plumeXY,plumeLbl,plumeXYtot,plumeBoundary]=deal(nan(max(ix),max(jx),nPlume));

for k=1:nPlume
    for i=1:length(S_bin_fine)
        plumeXY(ix(i),jx(i),k)=S_bin_fine(i,k);
    end
    plumeXYtot(:,:,k)=plumeXY(:,:,k);
    plumeXYtot(isnan(plumeXYtot))=0;
    plumeLbl(:,:,k)=bwlabel(plumeXYtot(:,:,k),4);
end

plumeLblTot=plumeLbl;
plumeLbl(isnan(plumeXY))=[];
plumeLbl=reshape(plumeLbl,[length(cModel.G.partition) nPlume]);


for k=1:nPlume
    plumeBoundary(:,:,k) = edge(plumeLblTot(:,:,k),'canny');
end
plumeBoundaryReg=plumeLblTot;
plumeBoundaryReg(plumeBoundary==0)=0;

plumeBoundaryReg(isnan(plumeXY))=[];
plumeBoundaryReg=reshape(plumeBoundaryReg,[],nPlume);

maxBndPts=round(0.05*cModel.G.parent.cells.num); %potential bug

[xW, yW]=deal(zeros(nWell,1));
[xB, yB, xBdiff, yBdiff, theta, dist]=deal(nan(maxBndPts,nWell,nPlume));
%% sampling
angSample=15;
xq=-180:angSample:180;
xq=deg2rad(xq)';
% plumeSample=cell(nPlume,1);
% plumeSampleTmp=[];

plumeData = nan(length(xq),2,nWell,nPlume);

% clf;
% hold on

for k=1:nPlume
% k=5;
    for i=1:nWell
        xW(i)=model.G.cells.centroids(W(i).cells,1);
        yW(i)=model.G.cells.centroids(W(i).cells,2);

        xBtemp=model.G.cells.centroids(plumeBoundaryReg(:,k)==i,1);
        xB(1:length(xBtemp),i,k)=xBtemp;

        yBtemp=model.G.cells.centroids(plumeBoundaryReg(:,k)==i,2);
        yB(1:length(yBtemp),i,k)=yBtemp;

        xBdiff(:,i,k)=xB(:,i,k)-xW(i);
        yBdiff(:,i,k)=yB(:,i,k)-yW(i);

        xTemp=xBdiff(:,i,k);
        yTemp=yBdiff(:,i,k);
        xTemp(isnan(xTemp))=[];
        yTemp(isnan(yTemp))=[];
        bndsz=length(xTemp);
        [theta(1:bndsz,i,k),dist(1:bndsz,i,k)]=cart2pol(xTemp,yTemp);

        thetaPlume=theta(:,i,k);
        distPlume=dist(:,i,k);
        thetaPlume(isnan(thetaPlume))=[];
        distPlume(isnan(distPlume))=[];

        vq=interp1(thetaPlume,distPlume,xq,'linear','extrap');

        %write plume data to matrix

        plumeData(:,1,i,k)=xq;
        plumeData(:,2,i,k)=vq;
 
    end
    
end
writematrix(plumeData,[pwd,'\',indir,'\','plumeData_',num2str(ens_id),'.txt'])

end

%% run prior ensemble

% prior_ens = ['prior_',casename];
% 
% priordir=['prior_',casename,'_result'];
% 
% parfor i = 1:nEnsemble %can be updated to parfor
%     fname=['ens-',num2str(i),'.txt'];
%     CGsimEnsemble(fname,prior_ens,priordir,t_train,cModelOri,cPredSched,cState0)
% end