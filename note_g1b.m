function testheo2()

% this function to show the spike train : 60 reps . for one cell on a 3 bins window.

load c:/bozeman/matlabfi/h1

h=figure('Color','w','PaperUnits','Centimeters','Units','centimeters','Position',[1 1 20 10]);

s_cellid=18;

Gss_cells=ss_cells;
Gv_thiscell_indice=s_cellid;
% what is the selected ini?
Gss_ini=ss_ini;
for s_i=1:length(Gss_ini)
   if (strmatch(Gss_ini(s_i).ini,Gss_cells(s_cellid).inifile,'exact'))
      % we find the indice in Gss_ini of the ini with the specified name, there should be only one
		Gs_thisini_indice=s_i;
   end;   
end; % for s_i     

Gs_xmin=350;
Gs_xmax=550;
%m_spike=shospik3(Gss_ini,Gss_cells,Gs_thisini_indice,Gv_thiscell_indice,Gs_xmin,Gs_xmax);

%save note_f5a m_spike;
load C:\MATLABR11\work\note_f5a

hold on;   
for s_i=1:60
   v_h=find(m_spike(s_i,:));
   for s_j=1:length(v_h)
   %h=stem(v_h,(s_i-1)*ones(size(v_h)),'b.');
   line([v_h(s_j) v_h(s_j)],[s_i-1 s_i],'Color','k','LineWidth',0.5);
	end;
end;
colormap(gray);

set(gca,'Units','Centimeters');
set(gca,'Position',[2 2 7 2]);
set(gca,'FontWeight','Bold','FontName','Arial','FontSize',8);






set(gca,'Units','Centimeters');
set(gca,'Position',[2 2 7 2]);
set(gca,'FontWeight','Bold','FontName','Arial','FontSize',8);

xlabel('time [ms]','FontWeight','Bold','FontName','Arial','FontSize',8);
%ylabel('repetitions','FontWeight','Bold','FontName','Arial','FontSize',8);

axis([0 6000 0 60])

line([0 6000],[0 0],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([0 6000],[60 60],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([6000 6000],[0 60],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([0 0],[0 60],'Color',[0.5 0.5 0.5],'LineWidth',2);

if (0)
   line([2350 2450],[-3 -3],'Color',[0.5 0.5 0.5],'LineWidth',2);
	line([2350 2450],[63 63],'Color',[0.5 0.5 0.5],'LineWidth',2);
	line([2450 2450],[-3 63],'Color',[0.5 0.5 0.5],'LineWidth',2);
	line([2350 2350],[-3 63],'Color',[0.5 0.5 0.5],'LineWidth',2);
end;

if (0)
h=figure('Color','w','PaperUnits','Centimeters','Units','centimeters','Position',[1 1 10 10]);
v_rate=(1000/size(m_spike,2))*sum(m_spike)/60;
bar(0.5*(1:length(v_rate)),v_rate,1,'Color','k');
set(gca,'Units','Centimeters');
set(gca,'Position',[2 2 6 2]);
set(gca,'FontWeight','Bold','FontName','Arial','FontSize',14);
end;
% the figure pointed by the arrow is : note_f4f




h=figure('Color','w','PaperUnits','Centimeters','Units','centimeters','Position',[1 1 20 10]);

hold on;   
for s_i=1:60
   v_h=find(m_spike(s_i,:));
   for s_j=1:length(v_h)
   %h=stem(v_h,(s_i-1)*ones(size(v_h)),'b.');
   line([v_h(s_j) v_h(s_j)],[s_i-1 s_i],'Color','k','LineWidth',1);
	end;
end;
colormap(gray);

m_f=Gss_ini(Gs_thisini_indice).vecframe/30;
m_f=floor(m_f);
v_f=m_f(1:850);
v_i=find((v_f>3900) & (v_f<4100));
v_f=v_f(v_i);

set(gca,'Units','Centimeters');
set(gca,'Position',[2 2 4 2]);
set(gca,'FontWeight','Bold','FontName','Arial','FontSize',8);

xlabel('time [ms]','FontWeight','Bold','FontName','Arial','FontSize',8);
%ylabel('repetitions','FontWeight','Bold','FontName','Arial','FontSize',8);

axis([min(v_f) max(v_f) 0 60])

if (0)
text(max(v_f)+2,60,'2','FontWeight','Bold','FontName','Arial','FontSize',8);
text(max(v_f)+2,30,'1','FontWeight','Bold','FontName','Arial','FontSize',8);
text(max(v_f)+2,0,'0','FontWeight','Bold','FontName','Arial','FontSize',8);
text(max(v_f)+20,22,'MC & FF','FontWeight','Bold','FontName','Arial','FontSize',8,'Rotation',90);
end;

line([min(v_f) max(v_f)],[0 0],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([min(v_f) max(v_f)],[60 60],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([max(v_f) max(v_f)],[0 60],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([min(v_f) min(v_f)],[0 60],'Color',[0.5 0.5 0.5],'LineWidth',2);

for s_i=1:length(v_f)-1
	line([v_f(s_i) v_f(s_i)],[0 60],'Color',[0.5 0.5 0.5],'LineWidth',1);
   m_help=m_spike(:,v_f(s_i):v_f(s_i+1));
   v_help=sum(m_help');
   v_help=v_help';
   v_meancount(s_i)=mean(v_help);
   v_varcount(s_i)=var(v_help);
   s_mean=v_meancount(s_i);
   s_var=v_varcount(s_i);
   s_fano=s_var/s_mean;
   v_xx=linspace(v_f(s_i),v_f(s_i+1),6);
	%line([v_xx(2) v_xx(3)],[60*s_mean/2 60*s_mean/2],'Color','k','LineWidth',2);
	%line([v_xx(4) v_xx(5)],[60*s_fano/2 60*s_fano/2],'Color','k','LineWidth',4);
	%line([v_xx(1) v_xx(6)],[30 30],'Color',[0.5 0.5 0.5],'LineWidth',1);
   
end;

v_fano=v_varcount./v_meancount;
% this is to plot the rate and the Fano
h=figure('Color','w','PaperUnits','Centimeters','Units','centimeters','Position',[1 1 20 10]);

set(gca,'Units','Centimeters');
set(gca,'Position',[2 2 4 2]);
set(gca,'FontWeight','Bold','FontName','Arial','FontSize',8);

line([min(v_f) max(v_f)],[0 0],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([min(v_f) max(v_f)],[2 2],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([max(v_f) max(v_f)],[0 2],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([min(v_f) min(v_f)],[0 2],'Color',[0.5 0.5 0.5],'LineWidth',2);

axis([min(v_f) max(v_f) 0 2])

for s_i=1:length(v_f)-1
	line([v_f(s_i) v_f(s_i)],[0 2],'Color',[0.5 0.5 0.5],'LineWidth',1);
   m_help=m_spike(:,v_f(s_i):v_f(s_i+1));
   v_help=sum(m_help');
   v_help=v_help';
   v_meancount(s_i)=mean(v_help);
   v_varcount(s_i)=var(v_help);
   s_mean=v_meancount(s_i);
   s_var=v_varcount(s_i);
   s_fano=s_var/s_mean;
   v_xx=linspace(v_f(s_i),v_f(s_i+1),10);
	patch([v_xx(2) v_xx(5) v_xx(5) v_xx(2)],[0 0 s_mean s_mean],[0.5 0.5 0.5]);
	patch([v_xx(5) v_xx(9) v_xx(9) v_xx(5)],[0 0 s_fano s_fano],'k');
end;


h=figure('Color','w','PaperUnits','Centimeters','Units','centimeters','Position',[1 1 20 10]);

hold on;   
for s_i=1:60
   v_h=find(m_spike(s_i,:));
   for s_j=1:length(v_h)
   %h=stem(v_h,(s_i-1)*ones(size(v_h)),'b.');
   line([v_h(s_j) v_h(s_j)],[s_i-1 s_i],'Color','k','LineWidth',1);
	end;
end;
colormap(gray);

m_f=Gss_ini(Gs_thisini_indice).vecframe/30;
m_f=floor(m_f);
v_f=m_f(1:850);
v_i=find((v_f>2550) & (v_f<2750));
v_f=v_f(v_i);

set(gca,'Units','Centimeters');
set(gca,'Position',[2 2 4 2]);
set(gca,'FontWeight','Bold','FontName','Arial','FontSize',8);

xlabel('time [ms]','FontWeight','Bold','FontName','Arial','FontSize',8);
%ylabel('repetitions','FontWeight','Bold','FontName','Arial','FontSize',8);

axis([min(v_f) max(v_f) 0 60])

line([min(v_f) max(v_f)],[0 0],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([min(v_f) max(v_f)],[60 60],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([max(v_f) max(v_f)],[0 60],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([min(v_f) min(v_f)],[0 60],'Color',[0.5 0.5 0.5],'LineWidth',2);

for s_i=1:length(v_f)-1
	line([v_f(s_i) v_f(s_i)],[0 60],'Color',[0.5 0.5 0.5],'LineWidth',1);
   m_help=m_spike(:,v_f(s_i):v_f(s_i+1));
   v_help=sum(m_help');
   v_help=v_help';
   v_meancount(s_i)=mean(v_help);
   v_varcount(s_i)=var(v_help);
   s_mean=v_meancount(s_i);
   s_var=v_varcount(s_i);
   s_fano=s_var/s_mean;
   v_xx=linspace(v_f(s_i),v_f(s_i+1),6);
%	line([v_xx(2) v_xx(3)],[60*s_mean/2 60*s_mean/2],'Color','k','LineWidth',2);
%	line([v_xx(4) v_xx(5)],[60*s_fano/2 60*s_fano/2],'Color','k','LineWidth',4);
%	line([v_xx(1) v_xx(6)],[30 30],'Color',[0.5 0.5 0.5],'LineWidth',1);
end;

v_fano=v_varcount./v_meancount;


set(gca,'Units','Centimeters');
set(gca,'Position',[2 2 4 2]);
set(gca,'FontWeight','Bold','FontName','Arial','FontSize',8);

xlabel('time [ms]','FontWeight','Bold','FontName','Arial','FontSize',8);
ylabel('repetitions','FontWeight','Bold','FontName','Arial','FontSize',8);

if (0)
text(max(v_f)+2,60,'2','FontWeight','Bold','FontName','Arial','FontSize',8);
text(max(v_f)+2,30,'1','FontWeight','Bold','FontName','Arial','FontSize',8);
text(max(v_f)+2,0,'0','FontWeight','Bold','FontName','Arial','FontSize',8);
text(max(v_f)+20,22,'MC & FF','FontWeight','Bold','FontName','Arial','FontSize',8,'Rotation',90);
end;

% this is to plot the rate and the Fano
h=figure('Color','w','PaperUnits','Centimeters','Units','centimeters','Position',[1 1 20 10]);

set(gca,'Units','Centimeters');
set(gca,'Position',[2 2 4 2]);
set(gca,'FontWeight','Bold','FontName','Arial','FontSize',8);

line([min(v_f) max(v_f)],[0 0],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([min(v_f) max(v_f)],[2 2],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([max(v_f) max(v_f)],[0 2],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([min(v_f) min(v_f)],[0 2],'Color',[0.5 0.5 0.5],'LineWidth',2);

axis([min(v_f) max(v_f) 0 2])

for s_i=1:length(v_f)-1
	line([v_f(s_i) v_f(s_i)],[0 2],'Color',[0.5 0.5 0.5],'LineWidth',1);
   m_help=m_spike(:,v_f(s_i):v_f(s_i+1));
   v_help=sum(m_help');
   v_help=v_help';
   v_meancount(s_i)=mean(v_help);
   v_varcount(s_i)=var(v_help);
   s_mean=v_meancount(s_i);
   s_var=v_varcount(s_i);
   s_fano=s_var/s_mean;
   v_xx=linspace(v_f(s_i),v_f(s_i+1),10);
	patch([v_xx(2) v_xx(5) v_xx(5) v_xx(2)],[0 0 s_mean s_mean],[0.5 0.5 0.5]);
	patch([v_xx(5) v_xx(9) v_xx(9) v_xx(5)],[0 0 s_fano s_fano],'k');
end;

save hola;