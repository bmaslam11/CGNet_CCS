function CGsimEnsembleFwdChk(ens_id,m_prior,opt,outdir,tRun,cModel,W_fine,cState0)
mrstModule add ad-core ad-blackoil deckformat ...
    agglom upscaling coarsegrid book ...
    mrst-gui ad-props incomp optimization...
    network-models test-suite linearsolvers

if ~exist(outdir, 'dir')
    mkdir(outdir)
end

%fname='caseInput.txt';
%!!important note!! make sure array dimension for pv, perm consisten with
%active cell num in coarse grid


schedule_full.control    = struct('W', W_fine);
schedule_full.step.control = ones(tRun, 1);
schedule_full.step.val = repmat(year, tRun, 1); %simulate for 50 years

cPredSched  = upscaleSchedule(cModel, schedule_full);

W=cPredSched.control.W;

params=m_prior(ens_id,:);

%% #############------------ ADD PARAMETER ------------##################
%nParams = 9;


%% update of model parameter
%params order must follow ensemble init
%---(1) permIJ
if opt.permIJ.flag
    id=opt.permIJ.size(1);
    permIJ=params(1,1:id);
    permIJ=transpose(permIJ);

    assert(size(permIJ,1)==size(cModel.rock.perm(:,1),1),'permIJ dimension mismatch')

    cModel.rock.perm(:,1) =permIJ;
    cModel.rock.perm(:,2) =permIJ;

    params(1:id)=[]; %delete params for next

end

%---(2) permK
if opt.permK.flag
    id=opt.permK.size(1);
    permK=params(1,1:id);
    permK=transpose(permK);

    assert(size(permK,1)==size(cModel.rock.perm(:,3),1),'permK dimension mismatch')
    
    cModel.rock.perm(:,3) =permK;
    params(1:id)=[]; %delete params for next
end

%---(3) porevolume
if opt.pv.flag
    id=opt.pv.size(1);
    pv=params(1,1:id);
    pv=transpose(pv);

    assert(size(pv,1)==size(cModel.operators.pv,1),'pv dimension mismatch')

    cModel.operators.pv = pv;
    params(1:id)=[]; %delete params for next
end


% ---(4) wellIndex
if opt.wi.flag
    id=opt.wi.size(1);
    wi=params(1,1:id);
    wi=transpose(wi);
    for i=1:length(cPredSched.control.W)

        %ind = find(strcmp({cPredSched.control.W.name}, wellName(i,:))==1);
        cPredSched.control.W(i).WI=wi(i); %include all wells

    end
    params(1:id)=[]; %delete params for next
end

%---(5) relperm
%reservoir RTYPE


satnumC=ones(cModel.G.cells.num,1); 
cModel.rock.regions = struct('saturation', satnumC); 

rtype=numel(cModel.rock.regions);


%Brooks-Corey function update
if opt.nw.flag
    id=opt.nw.size(1);
    nw=params(1,1:id);
    nw=transpose(nw);
    cModel.fluid.nW=nw;

    params(1:id)=[]; %delete params for next
end

if opt.ng.flag
    id=opt.ng.size(1);
    ng=params(1,1:id);
    ng=transpose(ng);
    cModel.fluid.nG=ng;

    params(1:id)=[]; %delete params for next
end

if opt.kw.flag
    id=opt.kw.size(1);
    kw=params(1,1:id);
    kw=transpose(kw);
    cModel.fluid.kW=kw;

    params(1:id)=[]; %delete params for next
end

if opt.kg.flag
    id=opt.kg.size(1);
    kg=params(1,1:id);
    kg=transpose(kg);
    cModel.fluid.kG=kg;

    params(1:id)=[]; %delete params for next
end

if opt.srw.flag
    id=opt.srw.size(1);
    srw=params(1,1:id);
    srw=transpose(srw);
    cModel.fluid.srW=srw;
end



[kRw,kRg]=deal(cell([1,rtype]));

for i=1:rtype
    kRw{i}=coreyPhaseRelpermAD(cModel.fluid.nW(i), cModel.fluid.srW(i), cModel.fluid.kW(i), cModel.fluid.srW(i));
    kRg{i}=coreyPhaseRelpermAD(cModel.fluid.nG(i), cModel.fluid.srG(i), cModel.fluid.kG(i), cModel.fluid.srW(i));
end

cModel.fluid.krG=kRg;
cModel.fluid.krW=kRw;


%modify and update CG model class to generic BO model
cModel = GenericBlackOilModel(cModel.G, cModel.rock,...
    cModel.fluid, 'gas', true, 'water', true, 'oil', false);
%% Simulate Model
[wellSols,states] = simulateScheduleAD(cState0, cModel, cPredSched);
%% Plot Check
%Compare Prediction and Reference Result

fh = plotWellSols(wellSols, ...
    cPredSched.step.val, ...
    'datasetnames', 'CGnet', 'zoom', false, ...
    'field', 'qWs');

figure
plotToolbar(cModel.G,states)
plotWell(cModel.G,W,'color','blue','color2','black','colorM','red')

conv_md=0.9869233*10^-16;

figure; plotToolbar(cModel.G, cModel.rock.perm(:,1)/conv_md); c=colorbar; colormap jet% perm
clim([min(cModel.rock.perm(:,1)/conv_md) max(cModel.rock.perm(:,1)/conv_md)])
%% ############------------GENERATE OUTPUT---------------##################

% casename=['ens_',num2str(ens_id)];
% 
% % writer
% %
% % wellSols
% nWell=length(W); %include injector&monitor
% 
% wellResult =zeros(1,(nWell*3+1)); %3:qGs, qWs, bhp +1 timestep
% 
% for i=1:tRun
%     t=i;
% 
%     qWs=vertcat(states{i}.wellSol(:).qWs);
%     qWs=transpose(qWs);
% 
%     qGs=vertcat(states{i}.wellSol(:).qGs);
%     qGs=transpose(qGs);
% 
%     bhp=vertcat(states{i}.wellSol(:).bhp);
%     bhp=transpose(bhp);
% 
%     soltmp=[t qWs qGs bhp];
%     wellResult=[wellResult; soltmp];
% 
% end
% 
% wellResult(1,:)=[];
% 
% nMap=cModel.G.cells.num;
% satgas=zeros(nMap,tRun);
% 
% for i=1:tRun
%     satgas(:,i)=states{i}.s(:,2);
% end
% 
% %delimiter : tab
% writematrix(satgas,[pwd,'\',outdir,'\','Sgas_',casename,'.txt'],'Delimiter','tab')
% 
% %writeplume


%% naming results

% [qWsnm, qGsnm, bhpnm]=deal(cell(nWell,1));
% AllWname=vertcat(W.name);
% 
% for i=1:nWell
% 
%     qWsnm{i}=append('qwS_',AllWname(i,:));
%     qGsnm{i}=append('qgS_',AllWname(i,:));
%     bhpnm{i}=append('BHP_',AllWname(i,:));
% 
% end
% 
% varname=[{'time_yrs'} qWsnm' qGsnm' bhpnm'];
% 
% Tresult=array2table(wellResult,"VariableNames",varname);
% writetable(Tresult,[pwd,'\',outdir,'\','wellResult_',casename,'.txt']);


function currentDir=getcurrentdir()
        if isdeployed %stand-alone mode
            %[status, result] = system('path');
            [~, result] = system('path');
            currentDir = char(regexpi(result, 'Path=(.*?);','tokens','once'));
        else
            currentDir = pwd;
        end
    end

end
