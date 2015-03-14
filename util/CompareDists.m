function [output] = CompareDists(model1,model2,name1,name2,...
  file_name,saveon,showJS,showKL,meanflag,covarflag,varflag)

% Compares two (sets of) probability distributions over binary vectors

% model1 and model2 should be structures containing (at least) the
% following fields:
% 'probs' = binvec probabilities
% 'prob_ebs' = error bars on those probabilities
% 'samp_size' = number of samples in original (unssubsampled) data set
% used to estimate probs
% 'subsamp_size' = size of subsampled data sets used to estimate probs and
% 'subsamp_probs' = matrix of probabilities estimated on data subsamplings

% name1 and name2 should be strings indicating the names of each model (to
% be displayed on output plots)

% file_name is a string which serves as the root for figure file names
% saveon is a boolean indicating whether or not to save figs (default off)

% probs [=] num_patterns X 1
% subsamp_probs [=] num_patterns X num_subsamp

% defaults and pre-processing
if ~exist('name1','var')
  name1 = sprintf('model1');
end
if ~exist('name2','var')
  name2 = sprintf('model2');
end
if ~exist('saveon','var') % dont save figs by default
  saveon = false;
end
if ~exist('showJS','var') % show JS by default
  showJS = 1;
end
if ~exist('showKL','var') % show KL by default
  showKL = 1;
end
% if ~isfield('probs_10e6',model1(1))
%   if size(model1(1).probs_10e6,2)==1
%     model1(1).probs_10e6 = model1(1).subsamp_probs>10^-6;
%   end
% end
% if ~isfield('probs_10e6',model2(1))
%   if size(model2(1).probs_10e6,2)==1
%     model2(1).probs_10e6 = model2(1).subsamp_probs>10^-6;
%   end
% end
% if ~isfield('probs_10e2',model1(1))
%   if size(model1(1).probs_10e2,2)==1
%     model1(1).probs_10e2 = model1(1).subsamp_probs>10^-2;
%   end
% end
% if ~isfield('probs_10e2',model2(1))
%   if size(model2(1).probs_10e2,2)==1
%     model2(1).probs_10e2 = model2(1).subsamp_probs>10^-2;
%   end
% end
% if ~isfield('N10e6s',model1(1))
%   model1(1).N10e6s = sum(model1(1).probs_10e6,2);
% end
% if ~isfield('N10e6s',model2(1))
%   model2(1).N10e6s = sum(model2(1).probs_10e6,2);
% end
% if ~isfield('N10e2s',model1(1))
%   model1(1).N10e2s = sum(model1(1).probs_10e2,2);
% end
% if ~isfield('N10e2s',model2(1))
%   model2(1).N10e2s = sum(model2(1).probs_10e2,2);
% end
% if ~isfield('N10e6_mean',model1(1))
%   model1(1).N10e6_mean = mean(model1(1).N10e6s);
% end
% if ~isfield('N10e6_mean',model2(1))
%   model2(1).N10e6_mean = mean(model2(1).N10e6s);
% end
% if ~isfield('N10e2_mean',model1(1))
%   model1(1).N10e2_mean = mean(model1(1).N10e2s);
% end
% if ~isfield('N10e2_mean',model2(1))
%   model2(1).N10e2_mean = mean(model2(1).N10e2s);
% end
% if ~isfield('N10e6_std',model1(1))
%   model1(1).N10e6_std = std(model1(1).N10e6s);
% end
% if ~isfield('N10e6_std',model2(1))
%   model2(1).N10e6_std = std(model2(1).N10e6s);
% end
% if ~isfield('N10e2_std',model1(1))
%   model1(1).N10e2_std = std(model1(1).N10e2s);
% end
% if ~isfield('N10e2_std',model2(1))
%   model2(1).N10e2_std = std(model2(1).N10e2s);
% end
fs = 20; % font size
xfig = 1; % cm from left side of screen
yfig = 20; % cm from bottom of screen
wfig = 25; % width in cm of figures
hfig = 15; % height in cm of figures

% dim check
if not(size(model1(1).probs,1) == size(model2(1).probs,1))
  error('Models are not of equal size.')
end

num_patterns = size(model1(1).probs,1);
data_dim = round(log2(num_patterns));

% calculate mean and uncertainties for JS/KL divergences
num_subsamp1 = size(model1(1).subsamp_probs,2);
num_subsamp2 = size(model2(1).subsamp_probs,2);
js = zeros(num_subsamp1, num_subsamp2);
kl = zeros(num_subsamp1, num_subsamp2);

for i = 1:num_subsamp1
  for j = 1:num_subsamp2
    js(i,j) =...
      JSDiv(model1(1).subsamp_probs(:,i),model2(1).subsamp_probs(:,j));
    kl(i,j) =...
      KLDiv(model1(1).subsamp_probs(:,i),model2(1).subsamp_probs(:,j));
  end
end
clear i j;

js_vec = reshape(js,num_subsamp1*num_subsamp2,1);
kl_vec = reshape(kl,num_subsamp1*num_subsamp2,1);
mean_js = mean(js_vec);
mean_kl = mean(kl_vec);
std_js = std(js_vec);
std_kl = std(kl_vec);
perc_error_js = round(100*(std_js/mean_js)); % percentage error (rounded)
perc_error_kl = round(100*(std_kl/mean_kl));

% calculate common pattern overlap
shared_10e6 = model1(1).probs_10e6.*model2(1).probs_10e6;
Nshareds_10e6 = sum(shared_10e6,1)';
Nshared_10e6_mean = mean(Nshareds_10e6);
Nshared_10e6_std = std(Nshareds_10e6);
p1shareds_10e6 = Nshareds_10e6./model1(1).N10e6s;
p2shareds_10e6 = Nshareds_10e6./model2(1).N10e6s;
p1shared_10e6_mean = mean(p1shareds_10e6);
p2shared_10e6_mean = mean(p2shareds_10e6);
p1shared_10e6_std = std(p1shareds_10e6);
p2shared_10e6_std = std(p2shareds_10e6);
shared_10e2 = model1(1).probs_10e2.*model2(1).probs_10e2;
Nshareds_10e2 = sum(shared_10e2,1)';
Nshared_10e2_mean = mean(Nshareds_10e2);
Nshared_10e2_std = std(Nshareds_10e2);
p1shareds_10e2 = Nshareds_10e2./model1(1).N10e2s;
p2shareds_10e2 = Nshareds_10e2./model2(1).N10e2s;
p1shared_10e2_mean = mean(p1shareds_10e2);
p2shared_10e2_mean = mean(p2shareds_10e2);
p1shared_10e2_std = std(p1shareds_10e2);
p2shared_10e2_std = std(p2shareds_10e2);

%% compare probs

min_prob = min(cat(1,...
  model1(1).probs-model1(1).prob_ebs,...
  model2(1).probs-model2(1).prob_ebs));
max_prob = max(cat(1,...
  model1(1).probs-model1(1).prob_ebs,...
  model2(1).probs-model2(1).prob_ebs));
% can be used to set plotting bounds

h1 = figure(1);
set(h1,'units','centimeters',...
  'outerposition',[xfig yfig wfig 17+showJS*.5+showKL*.5])
set(gcf,'PaperPositionMode','auto');
hold off
errorbarxy(...
  model1(1).probs,model2(1).probs,...
  model1(1).prob_ebs,model2(1).prob_ebs,...
  [],[],... % ux=lx, uy=ly, so don't input
  10^-6,1);
prettyplot(fs)
xlabel([name1,' probabilities'])
ylabel([name2,' probabilities'])
% set(gca,'xscale','log');
% set(gca,'yscale','log');
axis([10^-6 1 10^-6 1])
if showJS
  [decplace_JS str_JS] = sigdig(std_js); % finds n: 10^n = 1st sigdig
  str(1) = {[
    sprintf(['JS = ',str_JS],mean_js),...
    char(177),...
    sprintf([str_JS,' (%i'],std_js,perc_error_js),...
    char(37),...
    sprintf(')')]};
end
if showKL
  ind = 1;
  if showJS
    ind = 2;
  end
  [decplace_KL str_KL] = sigdig(std_kl); % finds n: 10^n = 1st sigdig
  str(ind) = {[
    sprintf(['KL = ',str_KL],mean_kl),...
    char(177),...
    sprintf([str_KL,' (%i'],std_kl,perc_error_kl),...
    char(37),...
    sprintf(')')]};
end
if showJS||showKL
  title(str);
end
legend('probabilities','error bars','Location','NorthEastOutside');
legend('boxoff')
refline(1,0)
hold off
xs = [10^-6 10^-4 10^-2 10^0];
ys = [10^-6 10^-4 10^-2 10^0];
set(gca,'XTick',xs);
set(gca,'YTick',ys);
% ylab = get(gca,'YLabel');
% set(ylab,'Position',...
%   get(ylab,'Position')+[0 -.0005 0])
% draw grid lines by hand
miny = min(ys);
maxy = max(ys);
minx = min(xs);
maxx = max(xs);
for x = xs
  line([x x],[miny maxy],'LineStyle','--','LineWidth',1,'Color','k');
end
clear x;
for y = ys
  line([minx maxx],[y y],'LineStyle','--','LineWidth',1,'Color','k');
end
clear y xs ys;
% save
str1 = ...
  strcat(name1,...
  sprintf(...
  ' subsampled from original size (%i) to subsampled size (%i) %i times',...
  model1(1).samp_size, model1(1).subsamp_size, num_subsamp1));
disp(str1);
str2 = ...
  strcat(name2,...
  sprintf(...
  ' subsampled from original size (%i) to subsampled size (%i) %i times',...
  model2(1).samp_size, model2(1).subsamp_size, num_subsamp2));
disp(str2);
if saveon
  save_name = [file_name,'_probs'];
  if showJS
    save_name = [save_name,'JS'];
  end
  if showKL
    save_name = [save_name,'KL'];
  end
  saveas(h1,[save_name,'.fig']);
  export_fig(save_name,'-pdf','-eps','-transparent');
end

%% compare mean component frequencies

if meanflag
  % build matrix whose rows hold all possible patterns
  patterns = zeros(num_patterns, data_dim);
  for i = 1:num_patterns
    patterns(i,:) = int2bin(i-1, data_dim);
  end
  clear i;

  % calculate means
  model1means = model1(1).subsamp_probs' * patterns; % dim = num_sumbsamp1 X data_dim
  model2means = model2(1).subsamp_probs' * patterns; % dim = num_sumbsamp2 X data_dim
  mean_model1means = mean(model1means,1)'; % dim = data_dim X 1
  mean_model2means = mean(model2means,1)'; % dim = data_dim X 1
  std_model1means = std(model1means,0,1)'; % dim = data_dim X 1
  std_model2means = std(model2means,0,1)'; % dim = data_dim X 1

  min_mean = min(cat(1,mean_model1means',mean_model2means'));
  max_mean = max(cat(1,mean_model1means',mean_model2means'));

  % plot means
  h2 = figure(2);
  set(h2,'units','centimeters','outerposition',[xfig yfig 22 16])
  set(gcf, 'PaperPositionMode', 'auto');
  hold off
  errorbarxy(...
    mean_model1means*100,mean_model2means*100,...
    std_model1means*100,std_model2means*100,...
    [],[],... % ux=lx, uy=ly, so don't input
    0,100,... % lb/ub
    [],... % use default eb colors
    true); % linear
  prettyplot(fs)
  xlabel( [name1, ' firing rates (Hz)'] )
  ylabel( [name2, ' firing rates (Hz)'] )
  legend('means','error bars','Location','NorthEastOutside')
  legend('boxoff')
  axis([0 100 0 100])
  refline(1,0)
  % ylab = get(gca,'YLabel');
  % set(ylab,'Position',get(ylab,'Position')+[0 -.04 0])
  hold off
  xs = [0 20 40 60 80 100];
  ys = [0 20 40 60 80 100];
  set(gca,...
    'XTick',xs,...
    'XTickLabel',{'0','20','40','60','80','100'});
  set(gca,...
    'YTick',ys,...
    'YTickLabel',{'0','20','40','60','80','100'});
  % draw grid lines by hand
  miny = min(ys);
  maxy = max(ys);
  minx = min(xs);
  maxx = max(xs);
  for x = xs
    line([x x],[miny maxy],'LineStyle','--','LineWidth',1,'Color','k');
  end
  clear x;
  for y = ys
    line([minx maxx],[y y],'LineStyle','--','LineWidth',1,'Color','k');
  end
  clear y xs ys;
  % save
  if saveon
    saveas(h2,[file_name,'_means.fig']);
    export_fig([file_name,'_means'],'-pdf','-eps','-transparent');
  end
end

%% calculate covariances and variances

if covarflag||varflag
    sampled_cov1 = zeros((1/2)*data_dim*(data_dim-1),num_subsamp1);
  sampled_cov2 = zeros((1/2)*data_dim*(data_dim-1),num_subsamp2);
  sampled_var1 = zeros(data_dim,num_subsamp1);
  sampled_var2 = zeros(data_dim,num_subsamp2);

  % calculate off-diagonal covariance terms
  for l = 1:num_subsamp1
    count = 1;
    for i = 1:data_dim
      for j = i+1:data_dim
        for k = 1:num_patterns
          sampled_cov1(count,l) = sampled_cov1(count,l)...
            +model1(1).subsamp_probs(k,l)...
            *(patterns(k,i)-model1means(l,i))...
            *(patterns(k,j)-model1means(l,j));
        end
      count = count+1;
      end
    end
  end
  clear l i j k count;

  for l = 1:num_subsamp2
    count = 1;
    for i = 1:data_dim
      for j = i+1:data_dim
        for k = 1:num_patterns
          sampled_cov2(count,l) = sampled_cov2(count,l)...
            +model2(1).subsamp_probs(k,l)...
            *(patterns(k,i)-model2means(l,i))...
            *(patterns(k,j)-model2means(l,j));
        end
      count = count + 1;
      end
    end
  end
  clear l i j k count;

  % calculate variances

  for l = 1:num_subsamp1
    for i = 1:data_dim
      for k = 1:num_patterns
        sampled_var1(i,l) = sampled_var1(i,l)...
            +model1(1).subsamp_probs(k,l)...
            *(patterns(k,i)-model1means(l,i))^2;
      end
    end
  end
  clear l i k;

  for l = 1:num_subsamp2
    for i = 1:data_dim
      for k = 1:num_patterns
        sampled_var2(i,l) = sampled_var2(i,l)...
          +model2(1).subsamp_probs(k,l)...
          *(patterns(k,i)-model2means(l,i))^2;
      end
    end
  end
  clear l i k;

  % calculate mean and std for each covariance and variance term

  mean_cov1 = mean(sampled_cov1,2); % dim = (1/2) * data_dim * (data_dim-1) X 1
  mean_cov2 = mean(sampled_cov2,2); % dim = (1/2) * data_dim * (data_dim-1) X 1
  std_cov1 = std(sampled_cov1,0,2); % dim = (1/2) * data_dim * (data_dim-1) X 1
  std_cov2 = std(sampled_cov2,0,2); % dim = (1/2) * data_dim * (data_dim-1) X 1

  mean_var1 = mean(sampled_var1,2); % dim = data_dim X 1
  mean_var2 = mean(sampled_var2,2); % dim = data_dim X 1
  std_var1 = std(sampled_var1,0,2); % dim = data_dim X 1
  std_var2 = std(sampled_var2,0,2); % dim = data_dim X 1

  % find mean covariance and variance terms and mean stds
  mean_mean_cov1 = mean(mean_cov1);
  mean_mean_cov2 = mean(mean_cov2);
  mean_mean_var1 = mean(mean_var1);
  mean_mean_var2 = mean(mean_var2);

  std_mean_cov1 = std(mean_cov1);
  std_mean_cov2 = std(mean_cov2);
  std_mean_var1 = std(mean_var1);
  std_mean_var2 = std(mean_var2);

  % find bounds
  % min_cov = min(cat(1, mean_cov1vec-std_cov1vec, mean_cov2vec-std_cov2vec, mean_var1-std_var1, mean_var2-std_var2));
  % max_cov = max(cat(1, mean_cov1vec+std_cov1vec, mean_cov2vec+std_cov2vec, mean_var1+std_var1, mean_var2+std_var2));
  
end


%% plot covariances

if covarflag

  % plot covariances
  h3 = figure(3);
  set(h3,'units','centimeters','outerposition',[xfig yfig 23 16.5])
  set(gcf,'PaperPositionMode','auto');
  hold off
  errorbarxy (...
    mean_cov1,mean_cov2,...
    std_cov1,std_cov2,...
    [],[],... % ux=lx, uy=ly, so don't input
    -.25,.25,... % lb/ub
    [],... % use default eb colors
    true); % linear
  prettyplot(fs)
  % leg(1) = {strcat(...
  %     sprintf('mean1 = %f', mean_mean_cov1),...
  %     char(177),...
  %     sprintf('%f', std_mean_cov1))};
  % leg(2) = {strcat(...
  %     sprintf('mean2 = %f', mean_mean_cov2),...
  %     char(177),...
  %     sprintf('%f', std_mean_cov2))};
  % text(.02, -.2, leg, 'FontSize', 15);
  xlabel([name1,' covariances'])
  ylabel([name2,' covariances'])
  legend('covariances','error bars','Location','NorthEastOutside')
  legend('boxoff')
  axis([-.25 .25 -.25 .25])
  refline(1,0)
  % ylab = get(gca,'YLabel');
  % set(ylab,'Position',get(ylab,'Position')+[0 -.03 0])
  hold off
  xs = [-.2 -.1 0 .1 .2];
  ys = [-.2 -.1 0 .1 .2];
  set(gca,...
    'XTick',xs,...
    'XTickLabel',{'-0.2','-0.1','0','0.1','0.2'});
  set(gca,...
    'YTick',ys,...
    'YTickLabel',{'-0.2','-0.1','0','0.1','0.2'});
  % draw grid lines by hand
  miny = -.25;
  maxy = .25;
  minx = -.25;
  maxx = .25;
  for x = xs
    line([x x],[miny maxy],'LineStyle','--','LineWidth',1,'Color','k');
  end
  clear x;
  for y = ys
    line([minx maxx],[y y],'LineStyle','--','LineWidth',1,'Color','k');
  end
  clear y xs ys;
  % save
  if saveon
    saveas(h3,[file_name,'_covars.fig']);
    export_fig([file_name,'_covars'],'-pdf','-eps','-transparent');
  end
  
end

%% plot variances

if varflag
  h4 = figure(4);
  set(h4,'units','centimeters','outerposition',[xfig yfig 22.5 16.5])
  set(gcf,'PaperPositionMode','auto');
  hold off
  errorbarxy (...
    mean_var1,mean_var2,...
    std_var1,std_var2,...
    [],[],... % ux=lx, uy=ly, so don't input
    -.5,.5,... % lb/ub
    [],... % use default eb colors
    true); % linear
  prettyplot(fs)
  % leg(1) = {strcat(...
  %     sprintf('mean1 = %f', mean_mean_var1),...
  %     char(177),...
  %     sprintf('%f', std_mean_var1))};
  % leg(2) = {strcat(...
  %     sprintf('mean2 = %f', mean_mean_var2),...
  %     char(177),...
  %     sprintf('%f', std_mean_var2))};
  % text(.02, -.2, leg, 'FontSize', 15);
  xlabel([name1,' variances'])
  ylabel([name2,' variances'])
  legend('variances','error bars','Location','NorthEastOutside')
  legend('boxoff')
  axis([-.5 .5 -.5 .5])
  set(gca,...
    'XTick',[-.5 -.25 0 .25 .5],...
    'XTickLabel',{'-0.5','-0.25','0','0.25','0.5'});
  set(gca,...
    'YTick',[-.5 -.25 0 .25 .5],...
    'YTickLabel',{'-0.5','-0.25','0','0.25','0.5'});
  refline(1,0)
  % ylab = get(gca,'YLabel');
  % set(ylab,'Position',get(ylab,'Position')+[0 -.04 0])
  hold off
  % draw grid lines by hand
  xs = [-.5 -.25 0 .25 .5];
  ys = [-.5 -.25 0 .25 .5];
  miny = min(ys);
  maxy = max(ys);
  minx = min(xs);
  maxx = max(xs);
  for x = xs
    line([x x],[miny maxy],'LineStyle','--','LineWidth',1,'Color','k');
  end
  clear x;
  for y = ys
    line([minx maxx],[y y],'LineStyle','--','LineWidth',1,'Color','k');
  end
  clear y xs ys;
  % save
  if saveon
    saveas(h4,[file_name,'_vars.fig']);
    export_fig([file_name,'_vars'],'-pdf','-eps','-transparent');
  end
  
end

%% save important statistics

if (meanflag&&(covarflag||varflag))
  output = struct(...
    'name1',name1,...
    'name2',name2,...
    'js',js,...
    'kl',kl,...
    'js_vec',js_vec,...
    'kl_vec',kl_vec,...
    'mean_js',mean_js,...
    'mean_kl',mean_kl,...
    'std_js',std_js,...
    'std_kl',std_kl,...
    'perc_error_js',perc_error_js,...
    'perc_error_kl',perc_error_kl,...
    'model1means',model1means,...
    'model2means',model2means,...
    'mean_model1means',mean_model1means,...
    'mean_model2means',mean_model2means,...
    'std_model1means',std_model1means,...
    'std_model2means',std_model2means,...
    'sampled_cov1',sampled_cov1,...
    'sampled_cov2',sampled_cov2,...
    'sampled_var1',sampled_var1,...
    'sampled_var2',sampled_var2,...
    'mean_cov1',mean_cov1,...
    'mean_cov2',mean_cov2,...
    'std_cov1',std_cov1,...
    'std_cov2',std_cov2,...
    'mean_var1',mean_var1,...
    'mean_var2',mean_var2,...
    'std_var1',std_var1,...
    'std_var2',std_var2,...
    'mean_mean_cov1',mean_mean_cov1,...
    'mean_mean_cov2',mean_mean_cov2,...
    'mean_mean_var1',mean_mean_var1,...
    'mean_mean_var2',mean_mean_var2,...
    'std_mean_cov1',std_mean_cov1,...
    'std_mean_cov2',std_mean_cov2,...
    'std_mean_var1',std_mean_var1,...
    'std_mean_var2',std_mean_var2,...
    'shared_10e6',shared_10e6,...
    'Nshareds_10e6',Nshareds_10e6,...
    'Nshared_10e6_mean',Nshared_10e6_mean,...
    'Nshared_10e6_std',Nshared_10e6_std,...
    'p1shareds_10e6',p1shareds_10e6,...
    'p2shareds_10e6',p2shareds_10e6,...
    'p1shared_10e6_mean',p1shared_10e6_mean,...
    'p2shared_10e6_mean',p2shared_10e6_mean,...
    'p1shared_10e6_std',p1shared_10e6_std,...
    'p2shared_10e6_std',p2shared_10e6_std,...
    'shared_10e2',shared_10e2,...
    'Nshareds_10e2',Nshareds_10e2,...
    'Nshared_10e2_mean',Nshared_10e2_mean,...
    'Nshared_10e2_std',Nshared_10e2_std,...
    'p1shareds_10e2',p1shareds_10e2,...
    'p2shareds_10e2',p2shareds_10e2,...
    'p1shared_10e2_mean',p1shared_10e2_mean,...
    'p2shared_10e2_mean',p2shared_10e2_mean,...
    'p1shared_10e2_std',p1shared_10e2_std,...
    'p2shared_10e2_std',p2shared_10e2_std,...
    'hprob',h1,...
    'hmean',h2);
  
  if covarflag
    output.hcovar = h3;
  end
  if varflag
    output.hvar = h4;
  end
  
elseif ((~meanflag)&&(covarflag||varflag))
  output = struct(...
    'name1',name1,...
    'name2',name2,...
    'js',js,...
    'kl',kl,...
    'js_vec',js_vec,...
    'kl_vec',kl_vec,...
    'mean_js',mean_js,...
    'mean_kl',mean_kl,...
    'std_js',std_js,...
    'std_kl',std_kl,...
    'perc_error_js',perc_error_js,...
    'perc_error_kl',perc_error_kl,...
    'sampled_cov1',sampled_cov1,...
    'sampled_cov2',sampled_cov2,...
    'sampled_var1',sampled_var1,...
    'sampled_var2',sampled_var2,...
    'mean_cov1',mean_cov1,...
    'mean_cov2',mean_cov2,...
    'std_cov1',std_cov1,...
    'std_cov2',std_cov2,...
    'mean_var1',mean_var1,...
    'mean_var2',mean_var2,...
    'std_var1',std_var1,...
    'std_var2',std_var2,...
    'mean_mean_cov1',mean_mean_cov1,...
    'mean_mean_cov2',mean_mean_cov2,...
    'mean_mean_var1',mean_mean_var1,...
    'mean_mean_var2',mean_mean_var2,...
    'std_mean_cov1',std_mean_cov1,...
    'std_mean_cov2',std_mean_cov2,...
    'std_mean_var1',std_mean_var1,...
    'std_mean_var2',std_mean_var2,...
    'shared_10e6',shared_10e6,...
    'Nshareds_10e6',Nshareds_10e6,...
    'Nshared_10e6_mean',Nshared_10e6_mean,...
    'Nshared_10e6_std',Nshared_10e6_std,...
    'p1shareds_10e6',p1shareds_10e6,...
    'p2shareds_10e6',p2shareds_10e6,...
    'p1shared_10e6_mean',p1shared_10e6_mean,...
    'p2shared_10e6_mean',p2shared_10e6_mean,...
    'p1shared_10e6_std',p1shared_10e6_std,...
    'p2shared_10e6_std',p2shared_10e6_std,...
    'shared_10e2',shared_10e2,...
    'Nshareds_10e2',Nshareds_10e2,...
    'Nshared_10e2_mean',Nshared_10e2_mean,...
    'Nshared_10e2_std',Nshared_10e2_std,...
    'p1shareds_10e2',p1shareds_10e2,...
    'p2shareds_10e2',p2shareds_10e2,...
    'p1shared_10e2_mean',p1shared_10e2_mean,...
    'p2shared_10e2_mean',p2shared_10e2_mean,...
    'p1shared_10e2_std',p1shared_10e2_std,...
    'p2shared_10e2_std',p2shared_10e2_std,...
    'hprob',h1);
  
  if covarflag
    output.hcovar = h3;
  end
  if varflag
    output.hvar = h4;
  end
  
elseif (meanflag&&~(covarflag||varflag))
  output = struct(...
    'name1',name1,...
    'name2',name2,...
    'js',js,...
    'kl',kl,...
    'js_vec',js_vec,...
    'kl_vec',kl_vec,...
    'mean_js',mean_js,...
    'mean_kl',mean_kl,...
    'std_js',std_js,...
    'std_kl',std_kl,...
    'perc_error_js',perc_error_js,...
    'perc_error_kl',perc_error_kl,...
    'model1means',model1means,...
    'model2means',model2means,...
    'mean_model1means',mean_model1means,...
    'mean_model2means',mean_model2means,...
    'std_model1means',std_model1means,...
    'std_model2means',std_model2means,...
    'shared_10e6',shared_10e6,...
    'Nshareds_10e6',Nshareds_10e6,...
    'Nshared_10e6_mean',Nshared_10e6_mean,...
    'Nshared_10e6_std',Nshared_10e6_std,...
    'p1shareds_10e6',p1shareds_10e6,...
    'p2shareds_10e6',p2shareds_10e6,...
    'p1shared_10e6_mean',p1shared_10e6_mean,...
    'p2shared_10e6_mean',p2shared_10e6_mean,...
    'p1shared_10e6_std',p1shared_10e6_std,...
    'p2shared_10e6_std',p2shared_10e6_std,...
    'shared_10e2',shared_10e2,...
    'Nshareds_10e2',Nshareds_10e2,...
    'Nshared_10e2_mean',Nshared_10e2_mean,...
    'Nshared_10e2_std',Nshared_10e2_std,...
    'p1shareds_10e2',p1shareds_10e2,...
    'p2shareds_10e2',p2shareds_10e2,...
    'p1shared_10e2_mean',p1shared_10e2_mean,...
    'p2shared_10e2_mean',p2shared_10e2_mean,...
    'p1shared_10e2_std',p1shared_10e2_std,...
    'p2shared_10e2_std',p2shared_10e2_std,...
    'hprob',h1,...
    'hmean',h2);
  
else
  output = struct(...
    'name1',name1,...
    'name2',name2,...
    'js',js,...
    'kl',kl,...
    'js_vec',js_vec,...
    'kl_vec',kl_vec,...
    'mean_js',mean_js,...
    'mean_kl',mean_kl,...
    'std_js',std_js,...
    'std_kl',std_kl,...
    'perc_error_js',perc_error_js,...
    'perc_error_kl',perc_error_kl,...
    'shared_10e6',shared_10e6,...
    'Nshareds_10e6',Nshareds_10e6,...
    'Nshared_10e6_mean',Nshared_10e6_mean,...
    'Nshared_10e6_std',Nshared_10e6_std,...
    'p1shareds_10e6',p1shareds_10e6,...
    'p2shareds_10e6',p2shareds_10e6,...
    'p1shared_10e6_mean',p1shared_10e6_mean,...
    'p2shared_10e6_mean',p2shared_10e6_mean,...
    'p1shared_10e6_std',p1shared_10e6_std,...
    'p2shared_10e6_std',p2shared_10e6_std,...
    'shared_10e2',shared_10e2,...
    'Nshareds_10e2',Nshareds_10e2,...
    'Nshared_10e2_mean',Nshared_10e2_mean,...
    'Nshared_10e2_std',Nshared_10e2_std,...
    'p1shareds_10e2',p1shareds_10e2,...
    'p2shareds_10e2',p2shareds_10e2,...
    'p1shared_10e2_mean',p1shared_10e2_mean,...
    'p2shared_10e2_mean',p2shared_10e2_mean,...
    'p1shared_10e2_std',p1shared_10e2_std,...
    'p2shared_10e2_std',p2shared_10e2_std,...
    'hprob',h1);
end
  
end

