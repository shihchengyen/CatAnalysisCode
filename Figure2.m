% script for figure 2
if ispresent('fig2data.mat','file')
	load fig2data
else
	load jonscellsdata
	
	% set some constants
	ncells = 36;
	totalFrames = 27550;
	nocolor = 1;

	% column numbers
	cellC = 1;
	nframesC = 2;
	fanoC = 3;
	mscC = 4;
	varC = 5;
	trestd10C = 6;
	trestd5C = 7;
	trestd3C = 8;
	trestd2C = 9;
	trescore10C = 10;
	trescore5C = 11;
	trescore3C = 12;
	trescore2C = 13;

	% initialize array 
	framesdata = zeros(totalFrames,trescore2C);
	
	fstart = 1;
	for i = 1:ncells
		fend = fstart + jonscells(i).num_frames - 1;
		fs = fstart:fend;
		framesdata(fs,cellC) = jonscells(i).cell_num;
		framesdata(fs,nframesC) = jonscells(i).num_frames;
		framesdata(fs,fanoC) = jonscells(i).fanofactor';
		framesdata(fs,mscC) = jonscells(i).meanspikecount';
		framesdata(fs,varC) = jonscells(i).varspikecount';
		framesdata(fs,trestd10C) = jonscells(i).tres_std{1}';
		framesdata(fs,trestd5C) = jonscells(i).tres_std{2}';
		framesdata(fs,trestd3C) = jonscells(i).tres_std{3}';
		framesdata(fs,trestd2C) = jonscells(i).tres_std{4}';
		framesdata(fs,trescore10C) = jonscells(i).tresscores{1}';
		framesdata(fs,trescore5C) = jonscells(i).tresscores{2}';
		framesdata(fs,trescore3C) = jonscells(i).tresscores{3}';
		framesdata(fs,trescore2C) = jonscells(i).tresscores{4}';
		fstart = fend + 1;
	end
	
	% find rows that do not have nan in fanoC
	ff = find(isnan(framesdata(:,fanoC))==0);
	% find number of frames with valid fano
	nff = length(ff);
	
	% data without nan (there are still some nan's in trestd5C trestd3C trestd2C)
	fd = framesdata(ff,:);
		
	% find frames with msc >= 1
	msc = find(fd(:,mscC)>=1);
	% get number of frames with msc >= 1
	nmsc = length(msc);
	% get number of frames with msc >= 1 for each cell
	[cmsc,cnmsc] = unique(fd(msc,cellC));
	cnmsc =  [cnmsc(1); diff(cnmsc)];	
	
	% find frames with ff < 1
	ff1 = find(fd(:,fanoC)<1);
	% get number of frames with ff < 1
	nff1 = length(ff1);
	
	% find frames with msc >= 1 & ff < 1
	mscff = intersect(ff1,msc);
	nmscff = length(mscff);
	% find cells with msc >= 1 & ff < 1
	[cmscff,cnmscff] = unique(fd(mscff,cellC));
	% get number of cells
	ncmscff = length(cmscff);
	% get number of frames with msc >= 1 & ff < 1 for each cell
	cnmscff = [cnmscff(1); diff(cnmscff)];
	% find cells in cmsc that are not in cmscff
	diff1 = setdiff(cmsc,cmscff);
	% create matrix containing cell_num and number of frames
	diff2 = [cmscff cnmscff; diff1 zeros(size(diff1))];
	% sortrows so that we have a vector equivalent to cnmsc
	diff3 = sortrows(diff2,1);
	cnmscff = diff3(:,2);
	
	% find frames with trescore10 > 99
	tres = find(fd(:,trescore10C)>99);
	% find frames with msc >= 1 & trescore10 > 99
	msctres = intersect(msc,tres);
	% get number of frames with msc >= 1 & trescore10 > 99
	nmsctres = length(msctres);
	% find cells with msc >= 1 & trescore10 > 99
	[cmsctres,cnmsctres] = unique(fd(msctres,cellC));
	% get number of cells
	ncmsctres = length(cmsctres);
	% get number of frames with msc >= 1 & trescore10 > 99 for each cell
	cnmsctres = [cnmsctres(1); diff(cnmsctres)];
	% find cells in cmsc that are not in cmsctres
	diff1 = setdiff(cmsc,cmsctres);
	% create matrix containing cell_num and number of frames
	diff2 = [cmsctres cnmsctres; diff1 zeros(size(diff1))];
	% sortrows so that we have a vector equivalent to cnmsc
	diff3 = sortrows(diff2,1);
	cnmsctres = diff3(:,2);
	
	% find frames with msc >= 1 & ff < 1 & trescore10 > 99
	mscfftres = intersect(mscff,tres);
	% get number of frames
	nmscfftres = length(mscfftres);
	% find cells with msc >= 1 & ff < 1 & trescore10 > 99
	[cmscfftres,cnmscfftres] = unique(fd(mscfftres,cellC));
	% get number of cells
	ncmscfftres = length(cmscfftres);
	% get number of frames with msc >= 1 & ff < 1 & trescore10 > 99 for each cell
	cnmscfftres = [cnmscfftres(1); diff(cnmscfftres)];
	% find cells in cmsc that are not in cmscfftres
	diff1 = setdiff(cmsc,cmscfftres);
	% create matrix containing cell_num and number of frames
	diff2 = [cmscfftres cnmscfftres; diff1 zeros(size(diff1))];
	% sortrows so that we have a vector equivalent to cnmsc
	diff3 = sortrows(diff2,1);
	cnmscfftres = diff3(:,2);
	
	% put the various number of frames for each cell into a format that 
	% will be easier to plot
	dis = [cnmsc (cnmscff-cnmscfftres) cnmscfftres (cnmsctres-cnmscfftres)];
	[distrib,distribi] = sortrows(dis,1);	
	
	% compute histograms for FF
	ffbins = 0:0.5:6;
	histff = histcie(fd(msc,fanoC),ffbins);
	% compute histograms for FF & trescore10 > 99
	histfftres = histcie(fd(msctres,fanoC),ffbins);
	
	% compute histograms for trescore10
	trbins = 0:5:100;
	histtres = histcie(fd(:,trescore10C),trbins);
	% compute histograms for msc >= 1 & trescore10 > 99
	histmsctres = histcie(fd(msc,trescore10C),trbins);
end

subplot(3,2,1)
plot(fd(:,mscC),fd(:,fanoC),'Color',[0.5 0.5 0.5],'Marker','.','LineStyle','none')
ax = axis;
% draw line to indicate fano = 1
line(ax(1:2),[1 1],'Color','k')
% draw line to indicate msc = 1
line([1 1],ax(3:4),'Color','k')
xlabel('Mean spike count')
% ylabel('Fano Factor')
% create inset to show FF below 1
h = axes('Position',[0.23 0.8 0.22 0.115]);
plot(fd(:,mscC),fd(:,fanoC),'Color',[0.5 0.5 0.5],'Marker','.','LineStyle','none')
axis([0 5 0 1])
line([1 1],[0 1],'Color','k')
set(gca,'YTick',[0.5 1])

subplot(3,2,2)
plot(fd(:,mscC),fd(:,varC),'Color',[0.5 0.5 0.5],'Marker','.','LineStyle','none')
axis([0 3 0 1])
hold on
x = 0:0.1:1;
minvar = x .* (1 - x);
plot(x,minvar,'k')
plot(1+x,minvar,'k')
plot(2+x,minvar,'k')
line([0 1],[0 1],'Color','k')
xlabel('Mean spike count')
% ylabel('Variance')
hold off

subplot(3,2,3)
hb = bar(trbins,histtres,0.8,'histc');
set(hb,'FaceColor',[0.5 0.5 0.5]);
vert = get(hb,'vertices');
bv = 1;
vert(vert==0) = bv;
set(hb,'vertices',vert);
ticks = [1 10 100 1000 10000];
set(gca,'YScale','log','YTick',ticks,'YMinorTick','off');
axis([-1 101 1 10000]);
hold on
hb = bar(trbins,histmsctres,0.8,'histc');
set(hb,'FaceColor','k');
vert = get(hb,'vertices');
vert(vert==0) = bv;
set(hb,'vertices',vert);
hold off
xlabel('TRES')
% ylabel('# of frames')

subplot(3,2,4)
bv = 0.8;
hb1 = bar(distrib(:,1));
vert = get(hb1,'vertices');
vert(vert==0) = bv;
set(hb1,'vertices',vert);
set(gca,'YScale','log','YTick',[1 2 5 10 20 50 150],'YMinorTick','off');
axis([0 33 0 170])
hold on
hb2 = bar(distrib(:,2:4),'stacked');
if nocolor
	set(hb1,'FaceColor',[0.25 0.25 0.25]);
	set(hb2(1),'FaceColor',[0.75 0.75 0.75]);
	set(hb2(2),'FaceColor',[0 0 0]);
	set(hb2(3),'FaceColor',[1 1 1]);
else
	set(hb1,'FaceColor',[0.5 0.5 0.5]);
	set(hb2(2),'FaceColor','g')
end
vert = get(hb2,'vertices');
for i = 1:3
	v = vert{i};
	v(v==0) = bv;
	vert{i} = v;
	set(hb2(i),'vertices',vert{i});
end
xlabel('Cell number')
% ylabel('# of frames')
hold off

subplot(3,2,5)
hb = bar(ffbins,histff,0.8,'histc');
set(hb,'FaceColor',[0.5 0.5 0.5])
vert = get(hb,'vertices');
vert(vert==0) = bv;
set(hb,'vertices',vert);
set(gca,'YScale','log','YTick',[1 2 5 10 20 50 200],'YMinorTick','off');
hold on
hb = bar(ffbins,histfftres,0.8,'histc');
set(hb,'FaceColor','k');
vert = get(hb,'vertices');
vert(vert==0) = bv;
set(hb,'vertices',vert);
hold off
axis([-0.1 6.1 0.8 300])
xlabel('Fano Factor')
% ylabel('# of frames')
set(gca,'XTick',0:1:6)

subplot(3,2,6)
plot(fd(msc,trestd10C),fd(msc,fanoC),'Color',[0.5 0.5 0.5],'Marker','.','LineStyle','none')
% hold on
% [b,bint] = regress(fd(msc,fanoC),[ones(size(msc)) fd(msc,trestd10C)]);
% x = [-3 25];
% y = b(2) * x + b(1);
% line(x,y,'Color','r')
axis([-3 25 0 6])
xlabel('TREZ')
% ylabel('Fano Factor')
% hold off

% subplot(3,2,6)
% plot(fd(msc,mscC),fd(msc,trescore10C),'k.')
% ylim([-1 105])
% hold on
% [b,bint] = regress(fd(msc,trescore10C),[ones(size(msc)) fd(msc,mscC)]);
% x = [0 8];
% y = b(2) * x + b(1);
% line(x,y,'Color','r')
% xlabel('Mean Spike Count')
% ylabel('TRES')
% hold off

