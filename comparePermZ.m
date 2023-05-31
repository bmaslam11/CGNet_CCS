function comparePermZ(cModelOri,HMmodel)

tiledlayout('flow')

nbin=25;

permKref = cModelOri.rock.perm(:,3);
permKsim = HMmodel.rock.perm(:,3);

%log transform
permKrefLog=log10(permKref);
permKsimLog=log10(permKsim);

max_lim=max([permKrefLog;permKsimLog]);
min_lim=min([permKrefLog;permKsimLog]);

nexttile
plotToolbar(cModelOri.G, permKrefLog);
colorbar;
clim([min_lim max_lim])
title('Upscaled Reference Model (log Kz)')
%view(-55, 60);

nexttile
plotToolbar(HMmodel.G, permKsimLog);
colorbar;
clim([min_lim max_lim])
title('History Matched Model (log Kz)')
%view(-55, 60);

nexttile([1 2])
h1=histogram(permKref,nbin);
hold on
h2=histogram(permKsim,nbin);

h1.Normalization = 'probability';
%h1.BinWidth = 0.5*1e6;
h2.Normalization = 'probability';
%h2.BinWidth = 0.5*1e6;

xlabel('Permeability-Zdir, m^2')
ylabel('normalized frequency')
legend('reference','matched')
end

