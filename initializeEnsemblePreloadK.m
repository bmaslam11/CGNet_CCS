function [params,opt]=initializeEnsemblePreloadK(model,q,cModelOri,nEnsemble,opt,casename,...
    all_G,all_rock,PERMALL)
mrstModule add ad-core ad-blackoil deckformat ...
               agglom upscaling coarsegrid book ...
               mrst-gui ad-props incomp optimization...
               network-models test-suite linearsolvers


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
        perm_KL=PERMALL(:,:,:,i);
        perm_KL=perm_KL(:);
        perm_KL=perm_KL(all_G.cells.indexMap);

        perm_min = 100;      %md
        perm_max = 2000;    %md

        perm_KL=perm_min+((perm_max-perm_min)*(perm_KL-min(perm_KL))...
            /(max(perm_KL)-min(perm_KL)));

        conv_md=0.9869233*10^-16;

        perm_KL=perm_KL*conv_md;
        all_rock.perm=repmat(perm_KL,1,3);

        target_res = false(all_G.cartDims);
        target_res(:, :, 6) = true;
        ind_res = find(target_res(all_G.cells.indexMap));

        perm=all_rock.perm(ind_res,1);
        model.rock.perm(:,1:3)=repmat(perm,1,3);
        cModelOri = upscaleModelTPFA(model, q,'transFromRock',false);
        %         Kgauss=gaussianField(model.G.cartDims,...
        %             [opt.permIJ.min opt.permIJ.max],[15,15,15],15);
        %         Kgauss=Kgauss(:);
        %         Kgauss=Kgauss(model.G.cells.indexMap);
        %
        %         model.rock.perm(:,1:3)=repmat(Kgauss,1,3);
        %
        %         cModelOri = upscaleModelTPFA(model, q,'transFromRock',false);
        permIJ = cModelOri.rock.perm(:,1);
        
        fprintf(outID,'%s\n','permIJ');
        fprintf(outID,'%d\n',permIJ);
        fprintf(outID,'%s\n','\permIJ');
                    
        params{i}=horzcat(params{i},permIJ');
    end

    if opt.permK.flag

        permK=permIJ;
        
        fprintf(outID,'%s\n','permK');
        fprintf(outID,'%d\n',permK);
        fprintf(outID,'%s\n','\permK ');
        
        params{i}=horzcat(params{i},permK');
    end

    if opt.pv.flag 
        pvgauss=gaussianField([opt.NI opt.NJ opt.NK],...
            [opt.pv.min opt.pv.max],[15,15,15],5);
        pvgauss=pvgauss(:);
        pvgauss=pvgauss(1:cModelOri.G.cells.num);
        
        
        pv=pvgauss;
        pvTot_ens = sum(pv);
        
%         guarantee same total pv
        pv=pv.*(pvTot_ori/pvTot_ens);

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