figure;

x11=coh_fmrib;x12=coh_sri24;x13=coh_icbm81;x14=coh_eve;x15=coh_iit2;x16=coh_iit3;

x21=cou_fmrib;x22=cou_sri24;x23=cou_icbm81;x24=cou_eve;x25=cou_iit2;x26=cou_iit3;

x31=ovl_fmrib;x32=ovl_sri24;x33=ovl_icbm81;x34=ovl_eve;x35=ovl_iit2;x36=ovl_iit3;

x41=dted_fmrib;x42=dted_sri24;x43=dted_icbm81;x44=dted_eve;x45=dted_iit2;x46=dted_iit3;

x51=dved_fmrib;x52=dved_sri24;x53=dved_icbm81;x54=dved_eve;x55=dved_iit2;x56=dved_iit3;

x61=tvdt_fmrib;x62=tvdt_sri24;x63=tvdt_icbm81;x64=tvdt_eve;x65=tvdt_iit2;x66=tvdt_iit3;

subplot(2,3,1)
h1=plot(coh_x,x11/sum(x11),coh_x,x12/sum(x12),coh_x,x13/sum(x13),coh_x,x14/sum(x14),coh_x,x15/sum(x15),coh_x,x16/sum(x16));
set(h1,'LineWidth',1,{'Color'},{'k';'c';'m';'b';'g';'r'});%,{'LineStyle'},{'-';'-';'-';'-';'-';'-'},{'Color'},{[0.5,0.5,0.5];[0.25,0.25,0.25];'k';[0.75,0.75,0.75]});
set(gca,'YGrid','on','Box','off','FontSize',14,'FontWeight','bold','YLim',[0 0.075]);
%title('coh','FontSize',20,'FontWeight','bold');%title('coh(mm^2/sec)','FontSize',20,'FontWeight','bold');
xlabel('COH','FontSize',20,'FontWeight','bold');
ylabel('Rel. # of WM voxels','FontSize',16,'FontWeight','bold');
%{
subplot(2,3,2)
h2=plot(cou_x,x21/sum(x21),cou_x,x22/sum(x22),cou_x,x23/sum(x23),cou_x,x24/sum(x24),cou_x,x25/sum(x25),cou_x,x26/sum(x26));
set(h2,'LineWidth',1,{'LineStyle'},{'-';'-';'-';'-';'-';'-'},{'Color'},{[0.5,0.5,0.5];[0.25,0.25,0.25];'k';[0.75,0.75,0.75]});
set(gca,'YGrid','on','Box','off','FontSize',14,'FontWeight','bold','XLim',[0 90],'YLim',[0 0.055],'XTick',[0 30 60 90],'XTickLabel',[0 30 60 90]);
%title('cou','FontSize',20,'FontWeight','bold');
xlabel('COU(\circ)','FontSize',20,'FontWeight','bold');

subplot(2,3,3)
h3=plot(ovl_x,x31/sum(x31),ovl_x,x32/sum(x32),ovl_x,x33/sum(x33),ovl_x,x34/sum(x34),ovl_x,x35/sum(x35),ovl_x,x36/sum(x36));
set(h3,'LineWidth',1,{'LineStyle'},{'-';'-';'-';'-';'-';'-'},{'Color'},{[0.5,0.5,0.5];[0.25,0.25,0.25];'k';[0.75,0.75,0.75]});
set(gca,'YGrid','on','Box','off','FontSize',14,'FontWeight','bold','XLim',[0.2 1],'YLim',[0 0.055]);%,'XLim',[0.3 1]);
%title('ovl','FontSize',20,'FontWeight','bold');
xlabel('OVL','FontSize',20,'FontWeight','bold');

subplot(2,3,4)
h4=plot(dted_x,x41/sum(x41),dted_x,x42/sum(x42),dted_x,x43/sum(x43),dted_x,x44/sum(x44),dted_x,x45/sum(x45),dted_x,x46/sum(x46));
set(h4,'LineWidth',1,{'LineStyle'},{'-';'-';'-';'-';'-';'-'},{'Color'},{[0.5,0.5,0.5];[0.25,0.25,0.25];'k';[0.75,0.75,0.75]});
set(gca,'YGrid','on','Box','off','FontSize',14,'FontWeight','bold','XLim',[0 8e-4],'YLim',[0 0.05]);
ylabel('Rel. # of WM voxels','FontSize',16,'FontWeight','bold');
xlabel('DTED(mm^2/s)','FontSize',20,'FontWeight','bold');
subplot(2,3,5)

h5=plot(dved_x,x51/sum(x51),dved_x,x52/sum(x52),dved_x,x53/sum(x53),dved_x,x54/sum(x54),dved_x,x55/sum(x55),dved_x,x56/sum(x56));
set(h5,'LineWidth',1,{'LineStyle'},{'-';'-';'-';'-';'-';'-'},{'Color'},{[0.5,0.5,0.5];[0.25,0.25,0.25];'k';[0.75,0.75,0.75]});
set(gca,'YGrid','on','Box','off','FontSize',14,'FontWeight','bold','XLim',[0 8e-4],'YLim',[0 0.06]);
xlabel('DVED(mm^2/s)','FontSize',20,'FontWeight','bold');

subplot(2,3,6)
h6=plot(tvdt_x,x61/sum(x61),tvdt_x,x62/sum(x62),tvdt_x,x63/sum(x63),tvdt_x,x64/sum(x64),tvdt_x,x65/sum(x65),tvdt_x,x66/sum(x66));
set(h6,'LineWidth',1,{'LineStyle'},{'-';'-';'-';'-';'-';'-'},{'Color'},{[0.5,0.5,0.5];[0.25,0.25,0.25];'k';[0.75,0.75,0.75]});
set(gca,'YGrid','on','Box','off','FontSize',14,'FontWeight','bold','XLim',[0 3e-7],'YLim',[0 0.16]);
xlabel('TVDT(mm^4/s^2)','FontSize',20,'FontWeight','bold');
%}
legend('FMRIB58','SRI24','ICBM81','Eve','IIT3mean','IIT2mean');
