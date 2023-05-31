function plotPlumeProb(hm,obs_dir,ens_dir, t_train_plume,nEnsemble,k)

%observation data
obs_data = fullfile(pwd,obs_dir,'plumeData.txt');
P_obs = readmatrix(obs_data);
P_obs = reshape(P_obs,[],2,hm.nWell,hm.plum.freq);
P_obs = P_obs(:,:,:,1:t_train_plume);

sz_obs=[size(P_obs) t_train_plume];

P_ens =zeros(sz_obs);
%ensemble data
for i=1:nEnsemble
    ens_data = fullfile(pwd,ens_dir,['plumeData_',num2str(i),'.txt']);
    P_tmp= readmatrix(ens_data);
    P_tmp= reshape(P_tmp,[],2,hm.nWell,t_train_plume);
    P_ens(:,:,:,:,i)=P_tmp;
end

%P_ens(theta,r_plume,well_id,plume_time,ensemble_id)

[r10,r50,r90]=deal(zeros(sz_obs(1),t_train_plume));

%get probablistic boundary
for n=1:t_train_plume
    r10(:,n)=prctile(squeeze(P_ens(:,2,1,n,:)),10,2);
    r50(:,n)=prctile(squeeze(P_ens(:,2,1,n,:)),50,2);
    r90(:,n)=prctile(squeeze(P_ens(:,2,1,n,:)),90,2);

end

tiledlayout('flow')

%k = 6; %temp counter
xq=rad2deg(P_obs(:,1,1,1));

coord_hi = [xq,r90(:,k)];
coord_lo = [xq,r10(:,k)];
coord_combine = [coord_hi;flipud(coord_lo)];

nexttile
%Inj-1
hold on
scatter(xq,P_obs(:,2,1,k),'black','filled')
%plot(xq,squeeze(P_ens(:,2,1,k,:)),'Color',[.2 .2 .2],'LineWidth',1.5)
plot(xq,r50(:,k),'Color',[.2 .2 .2],'LineWidth',1.5)

fill(coord_combine(:,1),coord_combine(:,2),'b','FaceAlpha',0.1,'EdgeAlpha',0)
ylim([0 max(r90(:,k))])
xlim([-180 180])
xticks(-180:30:180)
xline(0,'LineStyle','--','Color','red','LineWidth',2,'DisplayName','Azimuth')
%legend('Location','bestoutside')
xlabel('\theta (degree)','FontWeight','bold')
ylabel('Distance (m)','FontWeight','bold')
box on
title('CO_2 Plume Distribution CG (Polar) Inj-1')
subtitle(['at time : ',num2str(k),' years'])
hold off

% nexttile
% %Inj-2
% hold on
% scatter(xq,P_obs(:,2,2,k),'black','filled')
% plot(xq,squeeze(P_ens(:,2,2,k,:)),'Color',[.7 .7 .7])
% xlim([-180 180])
% xticks(-180:30:180)
% xline(0,'LineStyle','--','Color','red','LineWidth',2,'DisplayName','Azimuth')
% %legend('Location','bestoutside')
% xlabel('\theta (degree)','FontWeight','bold')
% ylabel('Distance (m)','FontWeight','bold')
% box on
% title('CO_2 Plume Distribution CG (Polar) Inj-2')
% subtitle(['at time : ',num2str(k),' years'])
% hold off

end