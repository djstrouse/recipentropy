function [] = prettyplot(fs)
% prettyplot makes the current figure pretty
% be sure to also: label axes, add a legend, fix linewidth/markers

% INPUTS
% fs - font size (default = 20)

% init
if ~exist('fs','var')||isempty(fs)
  fs = 20;
end

set(gca,'FontSize',fs,'LineWidth',2)
set(gcf,'PaperPositionMode','auto')
set(gca,'layer','top');
box off
legend('boxoff')

end

