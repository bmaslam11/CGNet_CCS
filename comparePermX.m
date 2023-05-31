function comparePermX(cModelOri,HMmodel)

tiledlayout('flow')

nbin=25;

permXref = cModelOri.rock.perm(:,1);
permXsim = HMmodel.rock.perm(:,1);

%log transform
permXrefLog=log10(permXref);
permXsimLog=log10(permXsim);

max_lim=max([permXrefLog;permXsimLog]);
min_lim=min([permXrefLog;permXsimLog]);

nexttile
plotToolbar(cModelOri.G, permXrefLog);
colorbar;
clim([min_lim max_lim])
title('Upscaled Reference Model (log Kx)')
%view(-55, 60);

nexttile
plotToolbar(HMmodel.G, permXsimLog);
colorbar;
clim([min_lim max_lim])
title('History Matched Model (log Kx)')
%view(-55, 60);

nexttile([1 2])
h1=histogram(permXref,nbin);
hold on
h2=histogram(permXsim,nbin);

h1.Normalization = 'probability';
%h1.BinWidth = 0.5*1e6;
h2.Normalization = 'probability';
%h2.BinWidth = 0.5*1e6;

xlabel('Permeability-Xdir, m^2')
ylabel('normalized frequency')
legend('reference','matched')
end

