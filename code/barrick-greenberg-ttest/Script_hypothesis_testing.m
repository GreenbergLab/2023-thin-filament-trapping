%This MATLAB code is associated with the following manuscript: Barrick,
%S.K., S.R. Clippinger, L. Greenberg, M.J. Greenberg. 2019. Computational
%tool to study perturbations in muscle regulation and its application to 
%heart disease.

%This is a script for hypothesis testing. The input is two matrices, each
%of which contain the best-fit parameters from the real data (best_fit_params)
%and from the accompanying bootstrapping simulations (bootstrap_params) for
%one of the two data sets to be compared (e.g., WT and mutant).

%%

%Change the names on the right side to the names of the variables in the workspace.
A_best = [5.5]%best_fit_params_A;
%A_boot =  readmatrix("D:\washu-thin-filament-regulation\data\isometric-force-clamp\Unregulated_bootstrap.csv")%^bootstrap_params_A;
A_boot =  readmatrix("D:\washu-thin-filament-regulation\data\standard-trap\unregulated_attachment-time-bootstrap.csv")%^bootstrap_params_A;
B_best = [3.9]%best_fit_params_B;
%B_boot = readmatrix("D:\washu-thin-filament-regulation\data\isometric-force-clamp\regulated_bootstrap.csv")%bootstrap_params_B;
B_boot = readmatrix("D:\washu-thin-filament-regulation\data\standard-trap\regulated_attachment-time-bootstrap.csv")%bootstrap_params_B;
tag='attachment-time-k-stats'; %When saving the data, this will be the name of the output.

ind=1;  %This is the index of the parameter to be compared. For example, 
%if best_fit_params is ordered as (KW,KTnocal,KTcal,KTmidcal,nH,A,B,C), 
%then ind=1 selects KW, ind=2 selects KTnocal, etc.

alpha=0.05; %This is the threshold value for significance.
null=0; %This is the value of the null hypothesis. In most cases, this is the hypothesis that A=B, implying A-B=0.
tail=2; %This defines whether the p-value is obtained for a 1- or 2-tailed
%test. Use 2 by default.

clearvars -except best_fit_params* bootstrap_params* *_best *_boot ind tag alpha null tail

close all

A = A_boot(:, ind); %This selects the parameters to be examined.
B = B_boot(:,ind);

A(isnan(A))=[]; %This gets rid of blanks.
B(isnan(B))=[];

%This extracts the true values of A and B and calculates the test statistic
%(i.e., the difference between A and B).
out_data = [A_best(ind) B_best(ind)]; %Values obtained from raw data
out_stat_data = A_best(ind)-B_best(ind); %This is the test statistic value.

%This calculates the cumulative distributions and plots them in Figure 1.
A_index = (1:length(A))./length(A);
B_index = (1:length(B))./length(B);
cumulative_A = sort(A);
cumulative_B = sort(B);

figure(1)
hold off
plot(cumulative_A,A_index,'k')
hold on
plot(cumulative_B,B_index,'r')

%This saves the plots. It can be commented out if desired.
saveas(gcf,strcat(tag,'_cumulative'),'epsc'); %Option for saving cumulative distributions
saveas(gcf,strcat(tag,'_cumulative'),'fig')

%The following code calculates the test statistic for each bootstrapping
%simulation to evaluate whether there is a statistically significant difference
%in the means of A and B.

%This creates a matrix of the parameter values obtained from the
%bootstrapping simulations.
clear out_sim

out_sim=[A B];

%This calculates the test statistic for the simulated data.
out_stat_sim = out_sim(:,1)-out_sim(:,2);
if median(out_stat_sim) <= 0 %This flips the curve for consistency.
    out_stat_sim=-out_stat_sim;
    out_stat_data=-out_stat_data;
end

%This calculates the confidence intervals for both the data and the test
%statistic.
CI_stat = prctile(out_stat_sim,[100*alpha/tail,100*(1-alpha/tail)]);
CI_data = prctile(out_sim,[100*alpha/tail,100*(1-alpha/tail)]);

%This makes the figure (Figure 2) showing the bootstrapped test statistic.
%Histogram shows the test statistic for all simulations. Red lines indicate
%the bounds of the 95% CI. The blue line indicates the null hypothesis. The
%yellow line marks the measured test statistic. If the blue line is outside
%of the boundaries defined by the red lines, then the null hypothesis is rejected.
figure(2)
hold off
histogram(out_stat_sim(:,1),'Normalization','probability')
hold on
plot(out_stat_data(1)*[1,1],ylim,'y-','LineWidth',2);
plot(CI_stat(1)*[1,1],ylim,'r-','LineWidth',2);
plot(CI_stat(2)*[1,1],ylim,'r-','LineWidth',2);
plot(null*[1,1],ylim,'b-','LineWidth',2);
xlabel('Difference between means');

set(gcf,'renderer','Painters')
print(strcat(tag,'_hist'),'-depsc') %Option for saving the histograms
saveas(gcf,strcat(tag,'_hist'),'fig')

%Next we get the p-value from the cumulative distribution of the test statistic (Figure 3).

out_cumulative_mean_ind = (1:length(out_stat_sim(:,1)))./length(out_stat_sim(:,1));
cumulative_mean = sort(out_stat_sim(:,1));

figure(3)
hold off
plot(cumulative_mean,out_cumulative_mean_ind,'k')
hold on
h1=plot(out_stat_data(1)*[1,1],ylim,'y-','LineWidth',2);
h2=plot(CI_stat(1)*[1,1],ylim,'r-','LineWidth',2);
plot(CI_stat(2)*[1,1],ylim,'r-','LineWidth',2);
h3=plot(null*[1,1],ylim,'b-','LineWidth',2);
xlabel('Difference between means');

saveas(gcf,strcat(tag,'_pval'),'epsc'); %Optional for saving cumulative distributions
saveas(gcf,strcat(tag,'_pval'),'fig')

%This displays a table containing the true fit value and the -/+ values for
%the confidence interval (CI) for A and B, as well as the p-value for the
%comparison of A and B.
p_ind = find(cumulative_mean>=0,1);
p_value = out_cumulative_mean_ind(p_ind)*tail;

rowNames = {'value','-','+'};
colNames = {'A','B'};
output = [out_data(1),out_data(2);out_data(1)-CI_data(1),out_data(2)-CI_data(3);
CI_data(2)-out_data(1),CI_data(4)-out_data(2)];
output_params = array2table(output,'RowNames',rowNames,'VariableNames',colNames)
p_value