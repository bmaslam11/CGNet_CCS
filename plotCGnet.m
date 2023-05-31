function plotCGnet(model,cModel,schedule,cPredSched,q)
% Plot the CGNet model
%----------------------------------------------------%
%Inputs:
%model

fig1 = figure;
%original model
subplot(2,2,1)
G = model.G;
plotCellData(G, model.rock.poro,'EdgeColor','none');
plotWell(G, schedule.control(1).W,'color','k','FontSize',10);
view(85,65); 
axis tight off; set(gca,'Clipping',false); zoom(1.2);

%partitioned model
subplot(2,2,3), ax=gca;
nq = max(q);
G = model.G;
colormap(gca,tatarizeMap(nq));
explosionView(G,q);
set(ax.Children,'EdgeAlpha',.1); view(85,75); 
axis tight off; set(gca,'Clipping',false); zoom(1.2);

subplot(2,2,4), ax=gca;
CG = cModel.G;
N = getNeighbourship(CG);
A = getConnectivityMatrix(N,true,CG.cells.num);
network = graph(A-eye(CG.cells.num));
pg = plot(network,'LineWidth',2, ...
    'XData', CG.cells.centroids(:,1), ...
    'YData', CG.cells.centroids(:,2), ...
    'ZData', CG.cells.centroids(:,3));
labelnode(pg,1:nq,'');
cells=vertcat(cPredSched.control(1).W.cells);
names=rldecode({cPredSched.control(1).W.name},...
    vertcat(arrayfun(@(x) numel(x.cells), cPredSched.control(1).W)),2);
labelnode(pg,cells,names);
set(ax.Children,'NodeFontSize',10,'NodeFontWeight','bold');
plotGrid(G,'FaceColor','none','EdgeAlpha',.05);
view(85,75);
axis tight off; set(gca,'Clipping',false); zoom(1.2);

subplot(2,2,2), ax=gca;
pg = plot(network,'LineWidth',0.5,'Layout','circle');
labelnode(pg,1:nq,'');
labelnode(pg,cells,names);
set(ax.Children,'NodeFontSize',10,'NodeFontWeight','bold');
axis equal tight off; set(gca,'Clipping',false); zoom(1.2);

end 

