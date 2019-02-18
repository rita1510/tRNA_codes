close all
clear
dicer_trna = importdata('metagene_stretched_100_SRR1230840_1_trimmedx6_aligned_all_hg38_F4_s_rmdup.txt');
pol3_trna = importdata('metagene_stretched_100_Pol3_Oler_hek293_hg38_F4_s_rmdup.txt');
pol2_trna = importdata('metagene_stretched_100_Pol2_ctrl_hg38_F4_s_rmdup.txt');
ctrl_trna = importdata('metagene_stretched_100_Pol2_chip_hg38_F4_s_rmdup.txt');
dicer_trna = movmean(dicer_trna ,25,2);
pol3_trna = movmean(pol3_trna ,25,2);
pol2_trna = movmean(pol2_trna ,25,2);
ctrl_trna = movmean(ctrl_trna ,25,2);

S_trna=sum(pol3_trna,2);
[Y_trna,i_trna]=sort(S_trna,'descend');


bluemap = makeColorMap([1 1 1], [0 0 1], 20);
redmap = makeColorMap([1 1 1], [1 0 0], 20);
greenmap = makeColorMap([1 1 1], [0 0.5 0], 20);
graymap = makeColorMap([1 1 1], [0 0 0], 100);
magmap = makeColorMap([1 1 1], [1 0 1], 20);
purplemap = makeColorMap([1 1 1], [0.5 0 0.5], 20);

figure;
subplot(1,4,1)
image(pol3_trna(i_trna,:))
colormap(graymap)
colorbar('southoutside')

title('PolIII')


subplot(1,4,2)
image(dicer_trna(i_trna,:))
colorbar('southoutside')

title('Dicer')

subplot(1,4,3)
image(pol2_trna(i_trna,:))
colorbar('southoutside')

title('PolII')


subplot(1,4,4)
image(ctrl_trna(i_trna,:))
title('Control')
colorbar('southoutside')
filename24=strcat('Figures/All_heatmaps.png');
print('-dpng',filename24)



figure;
image(pol3_trna_rmdup(i_trna_rmdup,:))
colormap(graymap)
title('Pol3 tRNA (Duplicates removed)')
colorbar('southoutside')
filename25=strcat('Figures/pol3_trna_heat_rmdup.png');
print('-dpng',filename25)

figure;
image(dicer_trna_rmdup(i_trna_rmdup,:))
colormap(redmap)
title('Dicer tRNA (Duplicates removed)')
colorbar('southoutside')
filename26=strcat('Figures/dicer_trna_heat_rmdup.png');
print('-dpng',filename26)

figure;
image(pol2_trna_rmdup(i_trna_rmdup,:))
colormap(greenmap)
title('Pol2 tRNA (Duplicates removed)')
colorbar('southoutside')
filename27=strcat('Figures/pol2_trna_heat_rmdup.png');
print('-dpng',filename27)


figure;

image(ctrl_trna_rmdup(i_trna_rmdup,:))
colormap(bluemap)
title('Control tRNA (Duplicates removed)')
colorbar('southoutside')
filename28=strcat('Figures/ctrl_trna_heat_rmdup.png');
print('-dpng',filename28)



%{
figure;
x=-200:200;
plot(x,trimmean(dicer_trna,10),'r','LineWidth',3)
hold on;
plot(x,trimmean(pol3_trna,10),'k','LineWidth',3)
plot(x,trimmean(pol2_trna,10),'LineWidth',3,'Color',[0 .7 .5])
plot(x,trimmean(ctrl_trna,10),'b','LineWidth',3)
legend('Dicer','Pol3','Pol2','Control');
xlabel('Distance from tRNA 5 prime end ');
ylabel('Average coverage (without top and bottom 10%),RPM');
title('tRNA metagene')
ylim([0 50]);
filename1=strcat('Figures/tRNA_metagene.png');
print('-dpng',filename1);

figure;
plot(x,trimmean(pol2_trna,10),'LineWidth',3,'Color',[0 .7 .5])
hold on;
plot(x,trimmean(ctrl_trna,10),'b','LineWidth',3)
legend('Pol2','Control');
xlabel('Distance from tRNA 5 prime end ');
ylabel('Average coverage (without top and bottom 10%), RPM');
title('tRNA metagene')
ylim([0 1]);
filename2=strcat('Figures/tRNA_metagene_pol2_ctrl.png');
print('-dpng',filename2);

figure;
x=-200:200;
plot(x,trimmean(dicer_trna_rmdup,10),'r','LineWidth',3)
hold on;
plot(x,trimmean(pol3_trna_rmdup,10),'k','LineWidth',3)
plot(x,trimmean(pol2_trna_rmdup,10),'LineWidth',3,'Color',[0 .7 .5])
plot(x,trimmean(ctrl_trna_rmdup,10),'b','LineWidth',3)
legend('Dicer','Pol3','Pol2','Control');
xlabel('Distance from tRNA 5 prime end - Duplicte reads removed');
ylabel('Average coverage (without top and bottom 10%),RPM');
title('tRNA metagene')
ylim([0 10]);
filename3=strcat('Figures/tRNA_metagene_rmdup.png');
print('-dpng',filename3);

figure;
plot(x,trimmean(pol2_trna_rmdup,10),'LineWidth',3,'Color',[0 .7 .5])
hold on;
plot(x,trimmean(ctrl_trna_rmdup,10),'b','LineWidth',3);
legend('Pol2','Control');
xlabel('Distance from tRNA 5 prime end ');
ylabel('Average coverage (without top and bottom 10%), RPM');
title('tRNA metagene')
ylim([0 0.9]);
filename4=strcat('Figures/tRNA_metagene_rmdup_pol2_ctrl.png');
print('-dpng',filename4)


figure;
plot(x,trimmean(dicer_rep,10),'r','LineWidth',3)
hold on;
plot(x,trimmean(pol3_rep,10),'k','LineWidth',3)
plot(x,trimmean(pol2_rep,10),'LineWidth',3,'Color',[0 .7 .5])
plot(x,trimmean(ctrl_rep,10),'b','LineWidth',3)
legend('Dicer','Pol3','Pol2','Control');
xlabel('Distance from Repeat 5 prime end ','FontSize',12);
ylabel('Average coverage (without top and bottom 10%),RPM','FontSize',12);
title('Repeat (tRNA like) metagene','FontSize',12,'FontWeight','bold')
ylim([0 50]);
filename5=strcat('Figures/rep_metagene.png');
print('-dpng',filename5)

figure;
plot(x,trimmean(dicer_rep,10),'r','LineWidth',3)
hold on;
plot(x,trimmean(pol3_rep,10),'k','LineWidth',3)
plot(x,trimmean(pol2_rep,10),'LineWidth',3,'Color',[0 .7 .5])
plot(x,trimmean(ctrl_rep,10),'b','LineWidth',3)
legend('Dicer','Pol3','Pol2','Control');
xlabel('Distance from Repeat 5 prime end ','FontSize',12);
ylabel('Average coverage (without top and bottom 10%),RPM','FontSize',12);
title('Repeat (tRNA like) metagene','FontSize',12,'FontWeight','bold')
filename6=strcat('Figures/rep_metagene_ZOOM.png');
print('-dpng',filename6)


figure;
plot(x,trimmean(pol2_rep,10),'LineWidth',3,'Color',[0 .7 .5])
hold on;
plot(x,trimmean(ctrl_rep,10),'b','LineWidth',3)
ylim([0 1]);
legend('Pol2','Control');
xlabel('Distance from Repeat 5 prime end ','FontSize',12);
ylabel('Average coverage (without top and bottom 10%), RPM','FontSize',12);
title('Repeat (tRNA like) metagene','FontSize',12,'FontWeight','bold')
filename7=strcat('Figures/rep_metagene_pol2_ctrl.png');
print('-dpng',filename7)

figure;
plot(x,trimmean(pol2_rep,10),'LineWidth',3,'Color',[0 .7 .5])
hold on;
plot(x,trimmean(ctrl_rep,10),'b','LineWidth',3)
legend('Pol2','Control');
xlabel('Distance from Repeat 5 prime end ','FontSize',12);
ylabel('Average coverage (without top and bottom 10%), RPM','FontSize',12);
title('Repeat (tRNA like) metagene','FontSize',12,'FontWeight','bold')
filename8=strcat('Figures/rep_metagene_pol2_ctrl_ZOOM.png');
print('-dpng',filename8)

figure;
plot(x,trimmean(dicer_rep_rmdup,10),'r','LineWidth',3)
hold on;
plot(x,trimmean(pol3_rep_rmdup,10),'k','LineWidth',3)
plot(x,trimmean(pol2_rep_rmdup,10),'LineWidth',3,'Color',[0 .7 .5])
plot(x,trimmean(ctrl_rep_rmdup,10),'b','LineWidth',3)
legend('Dicer','Pol3','Pol2','Control');
xlabel('Distance from Repeat 5 prime end - Duplicte reads removed','FontSize',12);
ylabel('Average coverage (without top and bottom 10%),RPM','FontSize',12);
title('Repeat (tRNA like) metagene','FontSize',12,'FontWeight','bold')
ylim([0 10]);

filename9=strcat('Figures/rep_metagene_rmdup.png');
print('-dpng',filename9)

figure;
plot(x,trimmean(dicer_rep_rmdup,10),'r','LineWidth',3)
hold on;
plot(x,trimmean(pol3_rep_rmdup,10),'k','LineWidth',3)
plot(x,trimmean(pol2_rep_rmdup,10),'LineWidth',3,'Color',[0 .7 .5])
plot(x,trimmean(ctrl_rep_rmdup,10),'b','LineWidth',3)
legend('Dicer','Pol3','Pol2','Control');
xlabel('Distance from Repeat 5 prime end - Duplicte reads removed','FontSize',12);
ylabel('Average coverage (without top and bottom 10%),RPM','FontSize',12);
title('Repeat (tRNA like) metagene','FontSize',12,'FontWeight','bold')
filename10=strcat('Figures/rep_metagene_rmdup_ZOOM.png');
print('-dpng',filename10)

figure;
plot(x,trimmean(pol2_rep_rmdup,10),'LineWidth',3,'Color',[0 .7 .5])
hold on;
plot(x,trimmean(ctrl_rep_rmdup,10),'b','LineWidth',3)
ylim([0 0.9]);
legend('Pol2','Control');
xlabel('Distance from Repeat 5 prime end ','FontSize',12);
ylabel('Average coverage (without top and bottom 10%), RPM','FontSize',12);
title('Repeat (tRNA like) metagene','FontSize',12,'FontWeight','bold')
filename11=strcat('Figures/rep_metagene_rmdup_pol2_ctrl.png');
print('-dpng',filename11)

figure;
plot(x,trimmean(pol2_rep_rmdup,10),'LineWidth',3,'Color',[0 .7 .5])
hold on;
plot(x,trimmean(ctrl_rep_rmdup,10),'b','LineWidth',3)
legend('Pol2','Control');
xlabel('Distance from Repeat 5 prime end ','FontSize',12);
ylabel('Average coverage (without top and bottom 10%), RPM','FontSize',12);
title('Repeat (tRNA like) metagene','FontSize',12,'FontWeight','bold')
filename12=strcat('Figures/rep_metagene_rmdup_pol2_ctrl_ZOOM.png');
print('-dpng',filename12)

figure;
plot(x,trimmean(dicer_pseudo,10),'r','LineWidth',3)
hold on;
plot(x,trimmean(pol3_pseudo,10),'k','LineWidth',3)
plot(x,trimmean(pol2_pseudo,10),'LineWidth',3,'Color',[0 .7 .5])
plot(x,trimmean(ctrl_pseudo,10),'b','LineWidth',3)
legend('Dicer','Pol3','Pol2','Control');
xlabel('Distance from pseudo tRNA 5 prime end ','FontSize',12);
ylabel('Average coverage (without top and bottom 10%),RPM','FontSize',12);
title('Pseudo tRNA metagene','FontSize',12,'FontWeight','bold')
ylim([0 50]);
filename13=strcat('Figures/pseudo_metagene.png');
print('-dpng',filename13)

figure;
plot(x,trimmean(dicer_pseudo,10),'r','LineWidth',3)
hold on;
plot(x,trimmean(pol3_pseudo,10),'k','LineWidth',3)
plot(x,trimmean(pol2_pseudo,10),'LineWidth',3,'Color',[0 .7 .5])
plot(x,trimmean(ctrl_pseudo,10),'b','LineWidth',3)
legend('Dicer','Pol3','Pol2','Control');
xlabel('Distance from pseudo tRNA 5 prime end ','FontSize',12);
ylabel('Average coverage (without top and bottom 10%),RPM','FontSize',12);
title('Pseudo tRNA metagene','FontSize',12,'FontWeight','bold')
filename14=strcat('Figures/pseudo_metagene_ZOOM.png');
print('-dpng',filename14)


figure;
plot(x,trimmean(pol2_pseudo,10),'LineWidth',3,'Color',[0 .7 .5])
hold on;
plot(x,trimmean(ctrl_pseudo,10),'b','LineWidth',3)
ylim([0 1]);
legend('Pol2','Control');
xlabel('Distance from pseudo tRNA 5 prime end ','FontSize',12);
ylabel('Average coverage (without top and bottom 10%), RPM','FontSize',12);
title('Pseudo tRNA metagene','FontSize',12,'FontWeight','bold')
filename15=strcat('Figures/pseudo_metagene_pol2_ctrl.png');
print('-dpng',filename15)

figure;
plot(x,trimmean(pol2_pseudo,10),'LineWidth',3,'Color',[0 .7 .5])
hold on;
plot(x,trimmean(ctrl_pseudo,10),'b','LineWidth',3);
legend('Pol2','Control');
xlabel('Distance from pseudo tRNA 5 prime end ','FontSize',12);
ylabel('Average coverage (without top and bottom 10%), RPM','FontSize',12);
title('Pseudo tRNA metagene','FontSize',12,'FontWeight','bold')
filename16=strcat('Figures/pseudo_metagene_pol2_ctrl_ZOOM.png');
print('-dpng',filename16)

figure;
plot(x,trimmean(dicer_pseudo_rmdup,10),'r','LineWidth',3)
hold on;
plot(x,trimmean(pol3_pseudo_rmdup,10),'k','LineWidth',3)
plot(x,trimmean(pol2_pseudo_rmdup,10),'LineWidth',3,'Color',[0 .7 .5])
plot(x,trimmean(ctrl_pseudo_rmdup,10),'b','LineWidth',3)
legend('Dicer','Pol3','Pol2','Control');
xlabel('Distance from pseudo tRNA 5 prime end - Duplicte reads removed','FontSize',12);
ylabel('Average coverage (without top and bottom 10%),RPM','FontSize',12);
title('Pseudo tRNA metagene','FontSize',12,'FontWeight','bold')
ylim([0 10]);
filename17=strcat('Figures/pseudo_metagene_rmdup.png');
print('-dpng',filename17)

figure;
x=-200:200;
plot(x,trimmean(dicer_pseudo_rmdup,10),'r','LineWidth',3)
hold on;
plot(x,trimmean(pol3_pseudo_rmdup,10),'k','LineWidth',3)
plot(x,trimmean(pol2_pseudo_rmdup,10),'LineWidth',3,'Color',[0 .7 .5])
plot(x,trimmean(ctrl_pseudo_rmdup,10),'b','LineWidth',3)
legend('Dicer','Pol3','Pol2','Control');
xlabel('Distance from pseudo tRNA 5 prime end - Duplicte reads removed','FontSize',12);
ylabel('Average coverage (without top and bottom 10%),RPM','FontSize',12);
title('Pseudo tRNA metagene','FontSize',12,'FontWeight','bold')
filename18=strcat('Figures/pseudo_metagene_rmdup_ZOOM.png');
print('-dpng',filename18)

figure;
plot(x,trimmean(pol2_pseudo_rmdup,10),'LineWidth',3,'Color',[0 .7 .5])
hold on;
plot(x,trimmean(ctrl_pseudo_rmdup,10),'b','LineWidth',3)
legend('Pol2','Control');
xlabel('Distance from pseudo tRNA 5 prime end ','FontSize',12);
ylabel('Average coverage (without top and bottom 10%), RPM','FontSize',12);
title('Pseudo tRNA metagene','FontSize',12,'FontWeight','bold')
ylim([0 0.9]);
filename19=strcat('Figures/pseudo_metagene_rmdup_pol2_ctrl.png');
print('-dpng',filename19)

figure;
plot(x,trimmean(pol2_pseudo_rmdup,10),'LineWidth',3,'Color',[0 .7 .5])
hold on;
plot(x,trimmean(ctrl_pseudo_rmdup,10),'b','LineWidth',3)
legend('Pol2','Control');
xlabel('Distance from pseudo tRNA 5 prime end ','FontSize',12);
ylabel('Average coverage (without top and bottom 10%), RPM','FontSize',12);
title('Pseudo tRNA metagene','FontSize',12,'FontWeight','bold')
filename20=strcat('Figures/pseudo_metagene_rmdup_pol2_ctrl_ZOOM.png');
print('-dpng',filename20)
%}
