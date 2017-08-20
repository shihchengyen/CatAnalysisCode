ps = ProcessDays(psparse,'IntraGroup','Cells',ndresp.SessionDirs,'AnalysisLevel','AllIntraGroup');

InspectGUI(ps)

figure; plot(ps)

% compare mean lifetime sparseness and median population sparseness
load sparsityobjs
% get indices grouped according to site
gind = groupDirs(sparseframe);
% find sites with more than 1 cell since we need more than 1 cell to
% compute population sparseness
si = find(sum(~isnan(gind))>1)
gind2 = gind(:,si);
% get the sparseness data grouped according to site
gdata = nanindex(sparseframe.data.sparsity,gind2)
% find mean sparseness for each site
gm = nanmean(gdata)
pm = prctile(ps.data.psparse,50);

% compute spree plot
% first compute variance for each cell
sfv = sparseframe.data.Values;
sv = nanstd(sfv).^2;
% now get variance for each site
vdata = nanindex(sv,gind2);
% find maximum variance
mv = nanmax(vdata);
% normalize by maximum variance
vdata2 = flipud(sort(vdata ./ repmat(mv,size(vdata,1),1)));

% data for 50 responsive cells
plot([nan nan 0 1],vdata2(:,[1 3 4 10 13]),'go-','MarkerSize',16)
hold on
plot([nan 0 0.5 1],vdata2(:,[2 5 6 7 9 11 12]),'bo-','MarkerSize',16)
plot([0 0.33 0.66 1],vdata2(:,8),'ro-','MarkerSize',16)
hold off

% data for 88 single units
subplot(1,6,1)
plot([nan nan nan nan nan 0 1],vdata2(:,[2 4 5 11 15 18 20]),'kx-','MarkerSize',16)
set(gca,'TickDir','out','Box','off')
% hold on
subplot(1,6,2)
plot([nan nan nan nan 0 0.5 1],vdata2(:,[3 10 12 14 16 21]),'kx-','MarkerSize',16)
set(gca,'YTickLabel','','TickDir','out','Box','off')
subplot(1,6,3)
plot([nan nan nan 0 0.33 0.66 1],vdata2(:,[1 17 19]),'kx-','MarkerSize',16)
set(gca,'YTickLabel','','TickDir','out','Box','off')
subplot(1,6,4)
plot([nan nan 0 0.25 0.50 0.75 1],vdata2(:,[6 13]),'kx-','MarkerSize',16)
set(gca,'YTickLabel','','TickDir','out','Box','off')
subplot(1,6,5)
plot([nan 0 0.2 0.4 0.6 0.8 1],vdata2(:,9),'kx-','MarkerSize',16)
set(gca,'YTickLabel','','TickDir','out','Box','off')
subplot(1,6,6)
plot([0 0.16 0.33 0.5 0.66 0.83 1],vdata2(:,[7 8]),'kx-','MarkerSize',16)
set(gca,'YTickLabel','','TickDir','out','Box','off')

% compute the area for 50 responsive cells
vd2 = vdata2(:,[1 3 4 10 13]);
vd3 = vdata2(:,[2 5 6 7 9 11 12]);
vd4 = vdata2(:,8);
a2 = 0.5 * (vd2(4,:) + 1);
a3 = 0.5 * repmat(0.5,1,2) * (vd3(2:3,:) + vd3(3:4,:));
a4 = 0.5 * repmat(1/3,1,3) * (vd4(1:3,:) + vd4(2:4,:));
screarea = [a2 a3 a4];
mean(screarea)
std(screarea)

% compute the area for 88 single units
vd2 = vdata2(:,[2 4 5 11 15 18 20]);
vd3 = vdata2(:,[3 10 12 14 16 21]);
vd4 = vdata2(:,[1 17 19]);
vd5 = vdata2(:,[6 13]);
vd6 = vdata2(:,9);
vd7 = vdata2(:,[7 8]);
a2 = 0.5 * (vd2(7,:) + 1);
a3 = 0.5 * repmat(0.5,1,2) * (vd3(5:6,:) + vd3(6:7,:));
a4 = 0.5 * repmat(1/3,1,3) * (vd4(4:6,:) + vd4(5:7,:));
a5 = 0.5 * repmat(1/4,1,4) * (vd5(3:6,:) + vd5(4:7,:));
a6 = 0.5 * repmat(1/5,1,5) * (vd6(2:6,:) + vd6(3:7,:));
a7 = 0.5 * repmat(1/6,1,6) * (vd7(1:6,:) + vd7(2:7,:));
screarea = [a2 a3 a4 a5 a6 a7];
mean(screarea)
std(screarea)

% compute scree plots by first normalizing each cell to its mean so that we
% can compare the variances across sites
% create normalized scree plot
sm = nanmean(sfv);
% normalize the spike counts by the mean so the variances 
% between sites are equivalent
sfv2 = sfv ./ repmat(sm,size(sfv,1),1);
sv2 = nanstd(sfv2).^2;
% now get variance for each site
vdata3 = nanindex(sv2,gind2);
% find maximum variance
mv2 = nanmax(vdata3);
% normalize by maximum variance
vdata4 = flipud(sort(vdata3 ./ repmat(mv2,size(vdata3,1),1)));

% plot for 50 responsive cells
plot([nan nan 0 1],vdata4(:,[1 3 4 10 13]),'k:')
hold on
plot([nan 0 0.5 1],vdata4(:,[2 5 6 7 9 11 12]),'k--')
plot([0 0.33 0.66 1],vdata4(:,8),'k-')
hold off

% plot for 88 single units
subplot(1,6,1)
plot([nan nan nan nan nan 0 1],vdata4(:,[2 4 5 11 15 18 20]),'kx-','MarkerSize',16)
set(gca,'TickDir','out','Box','off')
% hold on
subplot(1,6,2)
plot([nan nan nan nan 0 0.5 1],vdata4(:,[3 10 12 14 16 21]),'kx-','MarkerSize',16)
set(gca,'YTickLabel','','TickDir','out','Box','off')
subplot(1,6,3)
plot([nan nan nan 0 0.33 0.66 1],vdata4(:,[1 17 19]),'kx-','MarkerSize',16)
set(gca,'YTickLabel','','TickDir','out','Box','off')
subplot(1,6,4)
plot([nan nan 0 0.25 0.50 0.75 1],vdata4(:,[6 13]),'kx-','MarkerSize',16)
set(gca,'YTickLabel','','TickDir','out','Box','off')
subplot(1,6,5)
plot([nan 0 0.2 0.4 0.6 0.8 1],vdata4(:,9),'kx-','MarkerSize',16)
set(gca,'YTickLabel','','TickDir','out','Box','off')
subplot(1,6,6)
plot([0 0.16 0.33 0.5 0.66 0.83 1],vdata4(:,[7 8]),'kx-','MarkerSize',16)
set(gca,'YTickLabel','','TickDir','out','Box','off')

% compute the area for 50 responsive cells
vd2a = vdata4(:,[1 3 4 10 13]);
vd3a = vdata4(:,[2 5 6 7 9 11 12]);
vd4a = vdata4(:,8);
a2a = 0.5 * (vd2a(4,:) + 1);
a3a = 0.5 * repmat(0.5,1,2) * (vd3a(2:3,:) + vd3a(3:4,:))
a4a = 0.5 * repmat(1/3,1,3) * (vd4a(1:3,:) + vd4a(2:4,:))
screarea2 = [a2a a3a a4a];
mean(screarea2)
std(screarea2)

% compute the area for 88 single units
vd2a = vdata4(:,[2 4 5 11 15 18 20]);
vd3a = vdata4(:,[3 10 12 14 16 21]);
vd4a = vdata4(:,[1 17 19]);
vd5a = vdata4(:,[6 13]);
vd6a = vdata4(:,9);
vd7a = vdata4(:,[7 8]);
a2a = 0.5 * (vd2a(7,:) + 1);
a3a = 0.5 * repmat(0.5,1,2) * (vd3a(5:6,:) + vd3a(6:7,:));
a4a = 0.5 * repmat(1/3,1,3) * (vd4a(4:6,:) + vd4a(5:7,:));
a5a = 0.5 * repmat(1/4,1,4) * (vd5a(3:6,:) + vd5a(4:7,:));
a6a = 0.5 * repmat(1/5,1,5) * (vd6a(2:6,:) + vd6a(3:7,:));
a7a = 0.5 * repmat(1/6,1,6) * (vd7a(1:6,:) + vd7a(2:7,:));
screarea2 = [a2a a3a a4a a5a a6a a7a];
mean(screarea2)
std(screarea2)
