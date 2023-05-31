function comparePlume(cModel,cStateRef,cStateSim,tObs,sThreshold)

sCO2max_obs=cStateRef{tObs}.s(:,2);
sCO2max_sim=cStateSim{tObs}.s(:,2);


sCO2bin_obs=imbinarize(sCO2max_obs,sThreshold);
sCO2bin_sim=imbinarize(sCO2max_sim,sThreshold);
sCO2bin_diff=abs(sCO2bin_sim-sCO2bin_obs);

[cStateObs,cStateSim,cStateDiff]=deal(cStateRef{1});

cStateObs.s(:,2)=sCO2bin_obs;
cStateSim.s(:,2)=sCO2bin_sim;
cStateDiff.s(:,2)=sCO2bin_diff;

figure
tiledlayout(1,3)

nexttile
plotToolbar(cModel.G,cStateObs.s(:,2))
title('Reference CG','FontWeight','bold')
axis('square')

nexttile
plotToolbar(cModel.G,cStateSim.s(:,2))
title('Simulated CG','FontWeight','bold')
axis('square')

nexttile
plotToolbar(cModel.G,cStateDiff.s(:,2))
title('PlumeDiff','FontWeight','bold')
axis('square')

RMSE=sum(sCO2bin_diff)/length(sCO2bin_diff)

[A,ref]=deal(zeros(length(sCO2bin_sim),1));
A(sCO2bin_sim)=1;
ref(sCO2bin_obs)=1;

SSIM=ssim(A,ref)

end

