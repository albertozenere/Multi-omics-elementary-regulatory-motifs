

function print_density(Y,xlab,varargin)

assert(size(Y,2)<=2)
flag_two_densities = size(Y,2)==2;
if ~isempty(varargin)
    filename = varargin{1};
else
    filename = 'temp';
end

%% Settings
BandWidth = 0.03; %bandwidth for ksdensity
nbin = 100;

font = 'Helvetica';
sizeFont = 20;

x_lim = [-1 1];

widthLine = 1.5;

alphaFace = 0.2;

double2rgb = @(d) uint8(d);
edgeColor1 = double2rgb( [100,0,0] ); %dark red
edgeColor2 = double2rgb( [0,0,100] ); %dark blue
faceColor1 = double2rgb( [230, 50, 0] ); %less dark red
faceColor2 = double2rgb( [0,0,150] ); %less dark blue

%% Plot
% [f,x] = ksdensity(Y{1},'BandWidth',BandWidth);
% if flag_two_densities
%     xbin2 = (max(Y{2})-min(Y{2}))/150;
%     [f_cl,x_cl] = ksdensity(Y{2},'BandWidth',BandWidth);
% end
% % figure
% hold on
% a = area(x,f); a.FaceAlpha = alphaFace; a.LineWidth = widthLine; 
% a.EdgeColor = edgeColor2;  a.FaceColor = faceColor2;
% if flag_two_densities
% 	a = area(x_cl,f_cl); a.FaceAlpha = alphaFace; a.LineWidth = widthLine;
% 	a.EdgeColor = edgeColor1; a.FaceColor = faceColor1;
% end
% xlim(x_lim)
% xlabel(xlab)
% set(gca,'Fontsize',sizeFont, 'FontName', font,'LineWidth',widthLine);
% 
% print(filename,'-depsc') %depsc


%% Plot
[Ncohe, Ecohe] = histcounts(Y{1},nbin,'Normalization','probability');
 Xcohe = Ecohe(1:end-1)+diff(Ecohe);
 
if flag_two_densities
	[Ninco, Einco] = histcounts(Y{2},nbin,'Normalization','probability');
    Xinco = Einco(1:end-1)+diff(Einco);
end
% figure
hold on
a = area(Xcohe,smooth(Ncohe)); a.FaceAlpha = alphaFace; a.LineWidth = widthLine; 
a.EdgeColor = edgeColor2;  a.FaceColor = faceColor2;
if flag_two_densities
	a = area(Xinco,smooth(Ninco)); a.FaceAlpha = alphaFace; a.LineWidth = widthLine;
	a.EdgeColor = edgeColor1; a.FaceColor = faceColor1;
end
xlim(x_lim)
xlabel(xlab)
set(gca,'Fontsize',sizeFont, 'FontName', font, 'ytick', [],'LineWidth',widthLine);

print(filename,'-depsc') %depsc
