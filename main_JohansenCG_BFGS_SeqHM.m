%% Geosequestration Simulation in CGNet

%**************************************************************************
% Purpose           : Simulate CO2 Sequestration Case in Aquifer using 3D
%                     fine-grid model and coarse grid using CG net as proxy
%                     (Implemented using MRST)
% Program Features  : #Create proxy of fine grid model using CG net (faster)
%                     #BFGS solver is used to minimize misfit between CG and
%                     fine mode
%                     #Use ESMDA for uncertainty quantification
%                     (Next Step: Optimized Rate Allocation in CG model)
% Spv               : Dr. Bicheng Yan, KAUST DSFT PI
% Code Writer       : Billal Aslam (billal.aslam@itb.ac.id)- KAUST VS
%**************************************************************************
%% Dataset : Johansen Field Case
% modified from co2lab/example3D.m
clear; clc; close all;
% Load Modules

mrstModule add ad-core ad-blackoil deckformat ...
               agglom upscaling coarsegrid book ...
               mrst-gui ad-props incomp optimization...
               network-models test-suite linearsolvers
% Grid and Rock Properties Setup

% Load Johansen model

[all_G, all_rock, bcIx] = makeJohansenVEgrid('do_plot', false);

%perm_gauss
load("PERM\PERMALL_KL_n_100_100_11_101.mat")
nID = 51;
perm_KL=PERMALL(:,:,:,nID);
perm_KL=perm_KL(:);
perm_KL=perm_KL(all_G.cells.indexMap);

perm_min = 100;      %md
perm_max = 2000;    %md

perm_KL=perm_min+((perm_max-perm_min)*(perm_KL-min(perm_KL))...
    /(max(perm_KL)-min(perm_KL)));

poro_gauss=(1/0.283)*log(perm_KL/0.7583); 

% figure
% scatter(poro_gauss,perm_KL)
% set(gca,'YScale','log')

%conv mD to SI
conv_md=0.9869233*10^-16;

perm_KL=perm_KL*conv_md;


all_rock.perm=repmat(perm_KL,1,3);
%all_rock.poro=poro_gauss;

%% 
figure;
histogram(all_rock.poro)

figure;
histogram(all_rock.perm(:,1)/conv_md)

figure; plotToolbar(all_G, all_rock.perm(:,1)/conv_md); c=colorbar; colormap jet% porosity
c.Label.String='Permeability (md)';
figure; plotToolbar(all_G, all_rock.poro); colorbar; colormap jet% porosity


%% select layers
target_res = false(all_G.cartDims);
target_res(:, :, 6) = true;
ind_res = find(target_res(all_G.cells.indexMap));

G=extractSubgrid(all_G,ind_res);

rock.perm=all_rock.perm(ind_res,:);
rock.poro=all_rock.poro(ind_res);

% adjust rock prop
rock.poro=0.033*rock.poro;
rock.perm=1*rock.perm;

% min_K = min(rock.perm(:,1));
% max_K = max(rock.perm(:,1));
% 
% K=gaussianField(G.cartDims,[min_K max_K],[15,15,15],15);
% K=K(:);
% K=K(G.cells.indexMap);
% 
% rock.perm=repmat(K,1,3);


%%  Plot grid and rock properties (uncomment if needed)
close all
figure; plotToolbar(G, rock.perm(:,1)/conv_md); 
colorbar; colormap jet

figure; plotToolbar(G, rock.poro); 
colorbar; colormap jet

figure; scatter(rock.poro,rock.perm, 'filled')
set(gca,'YScale','log')
xlabel('porosity')
ylabel('permeability')


% figure; plotCellData(G, rock.poro); view(-35,15); colorbar; % porosity
% set(gcf, 'position', [531   337   923   356]); axis tight;


%% Prop Initial Condition Setup
gravity on; % tell MRST to turn on gravity
g = gravity; % get the gravity vector
rhow = 1000; % density of brine corresponding to 94 degrees C and 300 bar
initState.pressure = rhow * g(3) * G.cells.centroids(:,3);
initState.s = repmat([1, 0], G.cells.num, 1);
initState.sGmax = initState.s(:,2);
% Fluid (PVT) Model Setup

co2     = CO2props(); % load sampled tables of co2 fluid properties
p_ref   = 30 * mega * Pascal; % choose reference pressure
t_ref   = 94 + 273.15; % choose reference temperature, in Kelvin
rhoc    = co2.rho(p_ref, t_ref); % co2 density at ref. press/temp
cf_co2  = co2.rhoDP(p_ref, t_ref) / rhoc; % co2 compressibility
cf_wat  = 0; % brine compressibility (zero)
cf_rock = 4.35e-5 / barsa; % rock compressibility
muw     = 8e-4 * Pascal * second; % brine viscosity
muco2   = co2.mu(p_ref, t_ref) * Pascal * second; % co2 viscosity

mrstModule add ad-props; % The module where initSimpleADIFluid is found

% Use function 'initSimpleADIFluid' to make a simple fluid object
fluid = initSimpleADIFluid('phases', 'WG'           , ...
                           'mu'  , [muw, muco2]     , ...
                           'rho' , [rhow, rhoc]     , ...
                           'pRef', p_ref            , ...
                           'c'   , [cf_wat, cf_co2] , ...
                           'cR'  , cf_rock          , ...
                           'n'   , [2 2]); %see Corey model
% Rock Function Adjustment

% % relperm curves adjust end point
% srw = 0.27;
% src = 0.10;
% fluid.krW = @(s) fluid.krW(max((s-srw)./(1-srw), 0));
% fluid.krG = @(s) fluid.krG(max((s-src)./(1-src), 0));



%% Plot relperm curves
% figure; hold on;
% sw = linspace(0, 1, 200);
% plot(sw, fluid.krW(sw),   'b', 'linewidth', 1.5);
% plot(sw, fluid.krG(1-sw), 'r', 'linewidth', 1.5);
% line([srw, srw], [0 1], 'color', 'black', 'linestyle', ':', 'linewidth', 1);
% xlabel('brine saturation'); ylabel('relative permeability')
% set(gca, 'fontsize', 14, 'xlim', [0 1]);
% xticks((0:0.1:1))

%% Add capillary pressure curve
pe = 5 * kilo * Pascal;
pcWG = @(sw) pe * sw.^(-1/2);
fluid.pcWG = @(sg) pcWG(max((1-sg-srw)./(1-srw), 1e-5)); 
%% Add Wells

%Injector-1
% Well cell indices in 'global' grid: 48, 48, 6:10
wc_global = false(G.cartDims);
wc_global(48, 48, 6:10) = true;
wcI1 = find(wc_global(G.cells.indexMap));

%Injector-2
% wc_global = false(G.cartDims);
% wc_global(35, 68, 6:10) = true;
% wcI2 = find(wc_global(G.cells.indexMap));

%Monitor-1
wc_global = false(G.cartDims);
wc_global(25, 25, 2:10) = true;
wcM1 = find(wc_global(G.cells.indexMap));

%Monitor-2
wc_global = false(G.cartDims);
wc_global(15, 75, 2:10) = true;
wcM2 = find(wc_global(G.cells.indexMap));

%Monitor-3
wc_global = false(G.cartDims);
wc_global(52, 25, 2:10) = true;
wcM3 = find(wc_global(G.cells.indexMap));

%PressReliefProducer-1
wc_global = false(G.cartDims);
wc_global(5, 5, 2:6) = true;
wcP1 = find(wc_global(G.cells.indexMap));

%PressReliefProducer-2
wc_global = false(G.cartDims);
wc_global(5, 90, 2:6) = true;
wcP2 = find(wc_global(G.cells.indexMap));

%PressReliefProducer-3
wc_global = false(G.cartDims);
wc_global(53, 3, 2:6) = true;
wcP3 = find(wc_global(G.cells.indexMap));

%% Well Constraints

% Calculate the injection rate
inj_rate1 = 3.7 * mega * 1e3 / year / fluid.rhoGS;
%inj_rate2 = 3 * mega * 1e3 / year / fluid.rhoGS;

% Calculate the min bhp for producer
bhp1 = 0.9*min(initState.pressure);
bhp2 = 0.7*min(initState.pressure);
% Start with empty set of wells
W = [];

% Add a well to the set
W = addWell([], G, rock, wcI1, ...
            'name','I1', ...     % injector-1   
            'type', 'rate', ...  % inject at constant rate
            'val', inj_rate1, ...% volumetric injection rate
            'comp_i', [0 1]);    % inject CO2, not water

% W = addWell(W, G, rock, wcI2, ...
%             'name','I2', ...     % injector-2   
%             'type', 'rate', ...  % inject at constant rate
%             'val', inj_rate2, ...% volumetric injection rate (smaller)
%             'comp_i', [0 1]);    % inject CO2, not water

W = addWell(W, G, rock, wcM1, ...
            'name','M1', ...     % injector-3   
            'Sign',1, ...       % production ...
            'type', 'rate', ...  % inject at constant rate
            'val', 0, ...% shut in monitor
            'comp_i', [0 1]);    % inject CO2, not water

W = addWell(W, G, rock, wcM2, ...
            'name','M2', ...     % injector-3   
            'Sign',1, ...       % production ...
            'type', 'rate', ...  % inject at constant rate
            'val', 0, ...% shut in monitor
            'comp_i', [0 1]);    % inject CO2, not water

W = addWell(W, G, rock, wcM3, ...
            'name','M3', ...     % injector-3   
            'Sign',1, ...       % production ...
            'type', 'rate', ...  % inject at constant rate
            'val', 0, ...% shut in monitor
            'comp_i', [0 1]);    % inject CO2, not water

W = addWell(W, G, rock, wcP1, ...
            'name','P1', ...     % producer-1   
            'Sign',-1, ...       % production ...
            'type', 'bhp', ...   % produce at constant bhp
            'comp_i', [1 0], ...
            'val', bhp1);        % min bhp

W = addWell(W, G, rock, wcP2, ...
            'name','P2', ...     % producer-2   
            'Sign',-1, ...       % production  
            'type', 'bhp', ...   % produce at constant bhp
            'comp_i', [1 0], ...
            'val', bhp2);        % min bhp

W = addWell(W, G, rock, wcP3, ...
            'name','P3', ...     % producer-2   
            'Sign',-1, ...       % production  
            'type', 'bhp', ...   % produce at constant bhp
            'comp_i', [1 0], ...
            'val', bhp2);        % min bhp

%% Plot Well and Grid
clf;
plotGrid(G, 'facecolor', 'none', 'edgealpha', 0.1);
plotWell(G,W,'color','blue','color2','black','colorM','red')
view(30,15);

%%
convDarcy=9.869233*10^-13;
% porosity
figure; plotCellData(G, rock.perm(:,1)*1000/convDarcy); view(-35,15); colorbar;
plotWell(G,W,'color','blue','color2','black','colorM','red')
set(gcf, 'position', [531   337   923   356]); axis tight;

%% Get Cell Regions for Well Cell and Reservoir Cell

%reservoir RTYPE
satnum=ones(G.cells.num,1); 

%well RTYPE
% satnum(wcP1)=2;
% satnum(wcP2)=3;
% satnum(wcP3)=4;

rtype=numel(unique(satnum));
% assign different relperm using Corey utility function fn = coreyPhaseRelpermAD(n, sr, kwm, sr_tot)

% fluid.nW=[2 2 2 2]';
% fluid.nG=[2 2 2 2]';

fluid.nW=2;
fluid.nG=2;

% fluid.srW=[0.1 0.12 0.16 0.18]';
% fluid.srG=[0 0 0 0]';

fluid.srW=0.1;
fluid.srG=0;

% fluid.kW = [1 1 1 1]';
% fluid.kG = [0.8 0.7 0.6 0.5]';

fluid.kW = 1;
fluid.kG = 0.8;

[kRw,kRg]=deal(cell([1,rtype]));

for i=1:rtype
 kRw{i}=coreyPhaseRelpermAD(fluid.nW(i), fluid.srW(i), fluid.kW(i), fluid.srW(i));
 kRg{i}=coreyPhaseRelpermAD(fluid.nG(i), fluid.srG(i), fluid.kG(i), fluid.srW(i));
end

fluid.krG=kRg;
fluid.krW=kRw;
%  fluid.krG={@(S) S.^1.2, @(S) S.^1.5, @(S) S.^1.2, @(S) S.^2};
%  fluid.krW={@(S) S.^1.2, @(S) S.^1.5, @(S) S.^1.2, @(S) S.^2};

rock.regions = struct('saturation', satnum);

%% plot relperm
figure; hold on;
sw = linspace(0, 1, 100);

C={'k','b','r','g'};
for i=1:rtype
    fnW=fluid.krW{i};
    fnG=fluid.krG{i};

    plot(sw, fnW(sw),   C{i}, 'linewidth', 1.5,'LineStyle','-','DisplayName',['krw RT-',num2str(i)]);
    plot(sw, fnG(1-sw), C{i}, 'linewidth', 1.5,'LineStyle','--','DisplayName',['krg RT-',num2str(i)]);
    %line([srw, srw], [0 1], 'color', 'black', 'linestyle', ':', 'linewidth', 1);
end
legend(Location="best")

xlabel('brine saturation'); ylabel('relative permeability')
%% plot grid RTYPE
clf;
plotCellData(G,satnum)
plotWell(G,W,'color','blue','color2','black','colorM','red')

%% Simulation Schedule Setup
%%Schedule

schedule_obs.control    = struct('W', W);  %training data
schedule_full.control    = struct('W', W); %full data

% Specifying length of simulation timesteps 
tObs = 30;
tAll = 50;

schedule_obs.step.val = repmat(year, tObs, 1);
schedule_full.step.val = repmat(year, tAll, 1); %simulate for 50 years

%timesteps value
schedule_obs.step.control = ones(tObs, 1); %same constrain for tObs years
schedule_full.step.control = ones(tAll, 1); %same constrain for tAll years

%Generate Reservoir Model
model = GenericBlackOilModel(G, rock,...
    fluid, 'gas', true, 'water', true, 'oil', false);

%% Simulate Model
%% training
[wellSolObs, statesObs] = simulateScheduleAD(initState, model, schedule_obs);


%% reference
[wellSolFull, statesFull] = simulateScheduleAD(initState, model, schedule_full);

%% Check Result


plotToolbar(G,statesFull)


fh = plotWellSols({wellSolFull}, ...
    {schedule_full.step.val}, ...
    'datasetnames', {'Ori'}, 'zoom', true, ...
    'field', 'qWs');
%% Load Data if you already have saved previous simulation result!
%% Load Sim Full Data
%load("Fine_Basemodel_use\fine_basemodel.mat")
%% CGNet Model Conversion (Upscaling)
statesObs=statesFull(1:50);


predStates=statesObs; %training data


% make a coarse grid defined by a uniform NI x NJ x NK partition and
% perform a simple upscaling to obtain a coarse model
NI = 15;
NJ = NI;
NK = 1;


q = processPartition(model.G, partitionUI(model.G,[NI NJ NK]));
q = compressPartition(q);
cModel = upscaleModelTPFA(model, q,'transFromRock',false);

% We want to include rel-perm scaling as tunable parameters, so include
% these for the coarse model. These parameters have no effect for the
% initial coarse model (they are set equal to the ones given by the
% rel-perm curves).
% pts = cModel.fluid.krPts;
% scaling = {'SWL',   pts.w(1,1), 'SWCR', pts.w(1,2), 'SGU', pts.w(1,3), ...
%            'SGCR', pts.g(1,2), 'KRW',  pts.w(1,4), 'KRG', pts.g(1,4)};
% 
% cModel = imposeRelpermScaling(cModel, scaling{:});
% cModel = cModel.setupOperators();

cModel.toleranceCNV = 1e-6;  % tighter tolerance to improve gradient accuracy

%upscaling of initial condition
cState0 = upscaleState(cModel, model, initState);

%upscaling of reference state (for plume data)
cStateObs={};
for i=1:length(statesObs)
    cStateObs{i}=upscaleState(cModel,model,statesObs{i});
end 
cStateObs=transpose(cStateObs);

%modify CG model class to generic BO model
cModel = GenericBlackOilModel(cModel.G, cModel.rock,...
    cModel.fluid, 'gas', true, 'water', true, 'oil', false);

cModelOri = cModel;
%% check perm 
figure
plotToolbar(cModel.G,cModel.rock.perm(:,1)/conv_md); colorbar; colormap jet
caxis([perm_min perm_max])

figure
plotToolbar(G,model.rock.perm(:,1)/conv_md); colorbar; colormap jet
caxis([perm_min perm_max])
%%
convDarcy=9.869233*10^-13;
figure
hold on
h1=histogram(model.rock.perm(:,1)*1000/convDarcy);
h2=histogram(cModel.rock.perm(:,1)*1000/convDarcy);
h1.Normalization = 'probability';
h1.BinWidth = 10;
h2.Normalization = 'probability';
h2.BinWidth = 10;
xlabel('Permeability (milliDarcy)','FontWeight','bold')
ylabel('Probability','FontWeight','bold')
box on
legend('finegrid','upscaled CG')
hold off

%% Specify training schedule 

cPredSched  = upscaleSchedule(cModel, schedule_full);
cPredProbl  = struct('model', cModel, 'schedule', cPredSched, 'state0', cState0);
cPredProbl.state0.sGmax=cPredProbl.state0.s(:,2);

%% redo relperm region calc
%reservoir RTYPE
satnumC=ones(cModel.G.cells.num,1); 

% wellName=['P1';'P2';'P3'];
% 
% for i=1:length(wellName)
% 
%     ind = find(strcmp({cPredSched.control.W.name}, wellName(i,:))==1);
%     ind = cPredSched.control.W(ind).cells;
%     satnumC(ind)=i+1;
% end

cPredProbl.model.rock.regions = struct('saturation', satnumC);  
cModelOri.rock.regions = struct('saturation', satnumC);  
%% Run CGnet Basemodel
[wellSolCG, statesCG] = simulateScheduleAD(cState0, cModel, cPredSched);

%% Visualize CGNet and Fine Grid Model (Optional)

%Plot CGNet Model Conversion
plotCGnet(model,cModel,schedule_obs,cPredSched,q)

%% Compare Base and Upscaled
%Compare Prediction and Reference Result
fh = plotWellSols({wellSolFull,wellSolCG}, ...
    {schedule_full.step.val,cPredSched.step.val}, ...
    'datasetnames', {'reference','CGnet'}, 'zoom', true, ...
    'field', 'qWs');
%% 
figure
plotToolbar(G,statesFull)

figure
plotToolbar(cModel.G,statesCG)
%% Setup Calibration Parameters

% Choose parameters to match by include. See ModelParameter.m for details
config = {
    ...%name      include     scaling    boxlims  lumping   subset  relativeLimits  type
    'porevolume',       1,   'linear',       [],    [],      [],    [0.1   10],     'value'
    'conntrans',        1,   'log',          [],    [],      [],    [.001 100],     'value'
    'transmissibility', 1,   'log'   ,       [],    [],      [],    [.001 100],     'value'
    'permx',            1,   'log'   ,       [],    [],      [],    [.001 100],     'value'
    'permy',            0,   'log'   ,       [],    [],      [],    [.001 100],     'value'
    'permz',            0,   'log'   ,       [],    [],      [],    [.01  100],     'value'
    'kw',              0,   'linear',        [0.5 1],    [],      [],    [],        'value' 
    'kg',              0,   'linear',        [0.5 1],    [],      [],    [],        'value'  
    'srw',             0,   'linear',        [0 0.5],    [],      [],    [],        'value'
    'srg',             0,   'linear',        [],    [],      [],    [],             'value'
    'nw',              0,   'linear',        [1 6],    [],      [],    [],          'value'
    'ng',              0,   'linear',        [1 6],    [],      [],    [],          'value'};

predPrms = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    predPrms = addParameter(predPrms, cPredProbl, ...
        'name',    config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'lumping', config{k,5}, ...
        'subset',  config{k,6}, 'relativeLimits',config{k,7},...
        'type', config{k,8});
end
% History Match #1 - Matching Rates and Pressure Data

%History Match #1

% Deterministic History Matching using L-BFGS:

weighting  = {'WaterRateWeight',  day/150, ...
              'GasRateWeight',    day/150, ...
              'BHPWeight',        1/(400*barsa)};
pinit = getScaledParameterVector(cPredProbl, predPrms);


%% LM 
% mismatchFn = @(model, states, schedule, states_ref, tt, tstep, state) ...
%     matchObservedGWn(model, states, schedule, states_ref, ...
%                    'computePartials', tt, 'tstep', tstep, weighting{:},...
%                    'state', state, 'from_states', false, 'mismatchSum', false);
% objh = @(p) evaluateMatchSummands(p,mismatchFn,cPredProbl,predPrms,predStates);
% [v, popt, history] = unitBoxLM(pinit, objh, 'MaxIt', 11);

%% BFGS iteration loop
%error function
mismatchFn = @(xmodel, xstates, xschedule, states_ref, tt, tstep, state) ...
    matchObservedGW(xmodel, xstates, xschedule, states_ref,...
                   'computePartials', tt, 'tstep', tstep, weighting{:},...
                   'state', state, 'from_states', false);

objh = @(p)evaluateMatch(p,mismatchFn,cPredProbl,predPrms,predStates,'Gradient','AdjointAD');

[v, popt, history] = unitBoxBFGS(pinit, objh, 'objChangeTol', 1e-8, ...  %solving optimization step and return optimal parameter popt
    'gradTol', 1e-5, 'maxIt',25, 'lbfgsStrategy', 'dynamic', ...
    'lbfgsNum', 5, 'outputHessian',true, 'logPlot', true);
%% Run Prediction from HM#1 and compare to reference data

%update schedule
cPredSchedFull  = upscaleSchedule(cModel, schedule_full);
cPredProblFull  = struct('model', cModel, 'schedule', cPredSchedFull, 'state0', cState0);


[HMmodel,misfitP,~,wellSolP, HMstates] = ...
    evaluateMatchFullAll(popt,mismatchFn,cPredProblFull,predPrms,statesFull,'Gradient','none');
% QC: HM and prediction result with reference
%%
figure
plotToolbar(HMmodel.G,HMstates)
%% Get Observed data from HM-CG model

%%
%Compare properties with histogram

%Pore Volume
compareCellPV(cModelOri,HMmodel)
%%
%PermX
comparePermX(cModelOri,HMmodel)
%%
%PermZ
comparePermZ(cModelOri,HMmodel)
%% History Match Phase 2
% %Plume Compare
% sThresh=0.1;
% tObs=30;
% comparePlume(cModel,cStateObs,HMstates,tObs,sThresh)
% History Match #2 - Matching Plume Data

cPredSched  = upscaleSchedule(HMmodel, schedule_obs);
cPredProbl  = struct('model', HMmodel, 'schedule', cPredSched, 'state0', cState0);
cPredProbl.state0.sGmax=cPredProbl.state0.s(:,2);

% Choose parameters to match by include. See ModelParameter.m for details
config = {
    ...%name      include     scaling    boxlims  lumping   subset  relativeLimits  type
    'porevolume',       0,   'linear',       [],    [],      [],    [0.1   10],     'value'
    'conntrans',        0,   'log',          [],    [],      [],    [.001 100],     'value'
    'transmissibility', 0,   'log'   ,       [],    [],      [],    [.001 100],     'value'
    'permx',            0,   'log'   ,       [],    [],      [],    [.001 100],     'value'
    'permy',            0,   'log'   ,       [],    [],      [],    [.001 100],     'value'
    'permz',            0,   'log'   ,       [],    [],      [],    [.01  100],     'value'
    'kw',              1,   'linear',        [0.5 1],    [],      [],    [],        'value' 
    'kg',              1,   'linear',        [0.5 1],    [],      [],    [],        'value'  
    'srw',             1,   'linear',        [0 0.5],    [],      [],    [],        'value'
    'srg',             0,   'linear',        [],    [],      [],    [],             'value'
    'nw',              1,   'linear',        [1 6],    [],      [],    [],          'value'
    'ng',              1,   'linear',        [1 6],    [],      [],    [],          'value'};

predPrms = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    predPrms = addParameter(predPrms, cPredProbl, ...
        'name',    config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'lumping', config{k,5}, ...
        'subset',  config{k,6}, 'relativeLimits',config{k,7},...
        'type', config{k,8});
end

pinit = getScaledParameterVector(cPredProbl, predPrms);

% Deterministic History Matching using L-BFGS:

weighting  = {'WaterRateWeight',  day/150, ...
              'GasRateWeight',    day/150, ...
              'BHPWeight',        1/(400*barsa)};

%error function
mismatchFn = @(xmodel, xstates, xschedule, states_ref, tt, tstep, state) ...
    matchObservedGW(xmodel, xstates, xschedule, states_ref,...
                   'computePartials', tt, 'tstep', tstep, weighting{:},...
                   'state', state, 'from_states', false);

sThresh = 0.01;

objh = @(p)evaluateMatchWithPlumeAll(p,mismatchFn,cPredProbl,predPrms,predStates,cStateObs,sThresh,'Gradient','PerturbationADNUM');

%BFGS iteration loop
[v, popt, history] = unitBoxBFGS(pinit, objh, 'objChangeTol', 1e-8, ...  %solving optimization step and return optimal parameter popt
    'gradTol', 1e-5, 'maxIt',15, 'lbfgsStrategy', 'dynamic', ...
    'lbfgsNum', 5, 'outputHessian',true, 'logPlot', true);
%%
%update schedule
cPredSchedFull  = upscaleSchedule(HMmodel, schedule_full);
cPredProblFull  = struct('model', HMmodel, 'schedule', cPredSchedFull, 'state0', cState0);


[HMmodel,misfitP,~,wellSolP, HMstates] = ...
    evaluateMatchFullAll(popt,mismatchFn,cPredProblFull,predPrms,statesFull,'Gradient','none');
%%
%Plume Compare
sThresh=0.1;
tObs=30;
comparePlume(cModelOri,cStateObs,HMstates,tObs,sThresh)
%%
%Compare Prediction and Reference Result
fh = plotWellSols({wellSolFull,wellSolP}, ...
    {schedule_full.step.val,cPredSchedFull.step.val}, ...
    'datasetnames', {'reference','CGnet'}, 'zoom', true, ...
    'field', 'qWs');
%%
fh = plotWellSols({wellSolObsRtype,wellSolObs}, ...
    {schedule_obs.step.val,schedule_obs.step.val}, ...
    'datasetnames', {'Rtype','Ori'}, 'zoom', true, ...
    'field', 'qWs');
%%
x=HMmodel.operators.T;
y=cModelOri.fluid.kW;

scatterhist(x,y,'Kernel','on')
xlabel('calibratedCG')
ylabel('reference')
%%
tiledlayout(2,3)

nexttile
barh([HMmodel.fluid.kW,cModelOri.fluid.kW])
xline(1,'LineWidth',1.5,'LineStyle','--')
xlabel('kW')
ylabel('Rocktype')

nexttile
barh([HMmodel.fluid.kG,cModelOri.fluid.kG])
xline(0.65,'LineWidth',1.5,'LineStyle','--')
xlabel('kG')
ylabel('Rocktype')

nexttile
barh([HMmodel.fluid.nW,cModelOri.fluid.nW])
xline(1,'LineWidth',1.5,'LineStyle','--')
xlabel('nW')
ylabel('Rocktype')

nexttile
barh([HMmodel.fluid.nG,cModelOri.fluid.nG])
xline(1,'LineWidth',1.5,'LineStyle','--')
xlabel('nG')
ylabel('Rocktype')

nexttile
barh([HMmodel.fluid.srW,cModelOri.fluid.srW])
xline(0.15,'LineWidth',1.5,'LineStyle','--')
xlabel('srW')
ylabel('Rocktype')
legend('CG','reference','location','best')
