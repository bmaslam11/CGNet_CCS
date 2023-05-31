function compareCellPV(cModelOri,HMmodel)

tiledlayout('flow')
title('cellPV comparison')
nbin=25;

max_lim=max([cModelOri.operators.pv;HMmodel.operators.pv]);
min_lim=min([cModelOri.operators.pv;HMmodel.operators.pv]);

nexttile
plotToolbar(cModelOri.G, cModelOri.operators.pv);
colorbar;
clim([min_lim max_lim])
title('Upscaled Reference Model (PV)')
%view(-55, 60);

nexttile
plotToolbar(HMmodel.G, HMmodel.operators.pv);
colorbar;
clim([min_lim max_lim])
title('History Matched Model (PV)')
%view(-55, 60);

nexttile([1 2])
h1=histogram(cModelOri.operators.pv,nbin);
hold on
h2=histogram(HMmodel.operators.pv,nbin);

h1.Normalization = 'probability';
h1.BinWidth = 0.5*1e6;
h2.Normalization = 'probability';
h2.BinWidth = 0.5*1e6;

xlabel('Cell PV, m^3')
ylabel('normalized frequency')
legend('reference','matched')
end

