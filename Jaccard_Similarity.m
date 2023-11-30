%% Jaccard similirity of built models
% Name of the conditions   

% Reaction presence table, where columns are conditions, rows are reactions
model_keep

J = squareform(pdist(model_keep,'jaccard'));

altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0; 51 0 0]/255;

%get and save the default size
defaultPosition = get(0,'DefaultFigurePosition');
%get the current screen size
screensize = get( groot, 'Screensize' );
%screensize = get(0, 'Screensize'); %for earlier Matlab versions (e.g. Matlab 2010)
%set default figure position to full screen
screensize = [1,1,800,800];
set(0, 'DefaultFigurePosition', screensize);

cgo_J = clustergram(1-J,...
    'RowLabels', colnames,...
    'ColumnLabels', colnames,...
    'ColumnLabelsRotate',315, ...
    'Cluster', 'all', ...
    'symmetric','False',...
    'Colormap', altcolor)
addTitle(cgo_J,{"Model similarity using Jaccard distance"," based on the reconstructed models' reactions"})

%adding colorbar
cbButton = findall(gcf,'tag','HMInsertColorbar');
ccb = get(cbButton,'ClickedCallback')
set(cbButton,'State','on')
ccb{1}(cbButton,[],ccb{2})
%save plot
saveas(gcf,['Models_Jaccard_Similarity.png']);

% cgfig = findall(0,'type','figure', 'tag', 'Clustergram'); % Figure handle
% cgax = findobj(cgfig, 'type','axes','tag','HeatMapAxes'); % main axis handle
% % Change fontsize
% cgax.Fontsize = 8;

