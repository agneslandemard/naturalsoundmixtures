
function [fig1, fig2] = plot_in_CM_space(mods, fig_title)
if nargin < 2
    fig_title = '';
end
[vals, ~ ] = LoadCMparams;   
fig1 = figure('Position', [1000 1100 814 238]);
sz = [7 9 10];
mods = reshape(mods, sz);
spectrum = snm(mods,[1 2]);
nrj = snm(mods(:,:,:),3);

subplot(121);
plot(spectrum)
xticks(1:length(vals.freqs))
xticklabels(vals.freqs)
%caxis([-1 1])
%colormap('redblue')

subplot(122)
imagesc(nrj)
hold on;
plot(1+2*(snm(nrj,1)-min(snm(nrj,1)))./(max(snm(nrj,1))-min(snm(nrj,1))),'k')
plot(1+2*(snm(nrj,2)-min(snm(nrj,2)))./(max(snm(nrj,2))-min(snm(nrj,2))),1:sz(1),'k')
%clim([-1,1]*0.004)
colorbar
%colormap('redblue')
set(gca,'YDir','normal')
xticks(1:length(vals.rates));
yticks(1:length(vals.scales));
xticklabels(vals.rates);
yticklabels(vals.scales);
axis equal tight
ylabel('Scales (cyc/oct)')
xlabel('Rates (Hz)')
title(fig_title)


fig2 = figure('Position', [1000 1100 814 238]);
imagesc(mods(:,:))
vline(9.5:9:90.5,'k-')
hold on;
%clim([-1,1]*0.004)
colorbar
%colormap('redblue')
set(gca,'YDir','normal')
xticks(1:length(vals.rates));
yticks(1:length(vals.scales));
xticklabels(vals.rates);
yticklabels(vals.scales);
axis equal tight
title(fig_title)
ylabel('Scales (cyc/oct)')
xlabel('Rates (Hz)')

end