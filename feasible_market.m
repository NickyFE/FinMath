% FinMath Series Topic2 Feasible set of a market
% 
% this file serves as an example to construct all possible
% portfolios among multiple risky assets from a market, 
% and find out the minimum variance set and the tangent portfolio (market portfolio) 
% and also the efficient froniter (CML) with the existence of risk-free asset
%
% we use constituents from Dow Jones Industrial Average Index (DJIA) as an example
% the dataset of DJIA consists of 28 prominent companies listed in the US market
% we use monthly returns
% time range: July 2017 to June 2022
%
% Matlab built-in function used: fmincon
% try doc fmincon for more details on this function
%
% Copyright Nicky¿œ ¶ 2022
% nickycanhelp@qq.com
% Github: NickyFE/FinMath

clear;

%% read data from excel
[X,colname] = xlsread('DJIA_components_monthly_201707_202206.xlsx');
% in fact, DJIA has 30 components but in this dataset we have 28 
% because at the time when I prepare for this dataset, there are two stocks
% that have not been added into the index in the early days of our concerned time range

%% estimate sample mean and covariance from data
m=size(X,1); % number of observations (trading months)
n=size(X,2); % number of risky assets
all_one_m=ones(m,1); 
all_one_n=ones(n,1);

mu = (all_one_m'*X/m)'; % sample mean as a column vector
% or directly by:
% mu = mean(X,1)';
Sigma = 1/(m-1)*(X'*X-m*(mu*mu')); % sample covariance matrix
% or directly by:
% Sigma = cov(X);

% test existence of the inverse:
rank(Sigma) % should = n
% Sigma^(-1)

mu_F = 0.005; % set an artificial risk-free return
% mu_F < mu_MVP

%% other basic settings
% set right limit of x-axis in the figure:
xmax=0.15; % round(sqrt(max(diag(Sigma)))*1.1,2);

var_fn = @(w) w'*Sigma*w; % variance function of the portolio
mu_fn = @(w) w'*mu; % mean of the portfolio

w_0 = 1/n*ones(n,1);%[1;zeros(n-1,1)];
% initial guess on the weight when using fmincon below

shorting_yes = 1; % =1 if shorting is allowed, =0 if forbidden

% RGB codes for colors used:
color_blue = [0 0.4470 0.7410];
color_orange = [0.8500 0.3250 0.0980];
color_yellow = [0.9290 0.6940 0.1250];
color_purple = [0.4940 0.1840 0.5560];
color_green = [0.4660 0.6740 0.1880];
color_gray = [.85 .85 .85];
color_gray_dark = [.65 .65 .65];

%% random compute portfolios in order to plot later
num_P = 10000;
w_boundary = 10;
mu_rand = zeros(num_P,1);
sd_rand = zeros(num_P,1);
for i = 1:num_P
    if shorting_yes
        w_temp = zeros(n,1);
        while sum(w_temp)==0
            w_temp = -w_boundary+2*w_boundary*rand(n,1);
        end
    else
        w_temp = rand(n,1); 
    end
    w_rand = w_temp/sum(w_temp);

    mu_rand(i) = mu_fn(w_rand);
    sd_rand(i) = sqrt(var_fn(w_rand));
end

%% find MVP
options = optimoptions('fmincon','MaxFunEvals',5000);
Aeq_MVP = ones(1,n);
beq_MVP = 1;
w_MVP = fmincon(var_fn,w_0,[],[],Aeq_MVP,beq_MVP,[],[],[],options);
mu_MVP = w_MVP'*mu;
sd_MVP = sqrt(var_fn(w_MVP));

% check with analytical solution for MVP:
v = Sigma^(-1)*all_one_n;
w_MVP_analytical = v/(all_one_n'*v);
mu_MVP_analytical = w_MVP_analytical'*mu;
sd_MVP_analytical = sqrt(var_fn(w_MVP_analytical));

%% find min variance for each level of return z
zmin = -0.01;
zmax = 0.08;
z = zmin:0.0005:zmax;
num_level = length(z);
Aeq_z = [mu';ones(1,n)];
w_min_var = zeros(num_level,n);
mu_min_var = zeros(num_level,1);
sd_min_var = zeros(num_level,1);
options = optimoptions('fmincon','Display','off','MaxFunEvals',5000);
for i = 1:num_level
    beq_z = [z(i);1];
    w_temp = fmincon(var_fn,w_0,[],[],Aeq_z,beq_z,[],[],[],options);
    w_min_var(i,:) = w_temp';
    mu_min_var(i) = mu_fn(w_temp);
    sd_min_var(i) = sqrt(var_fn(w_temp));
end

%% find tangent portfolio
% in this code, we find tangent portfolio through fmincon
% instead of directly invoke the analytical solution
obj = @(w) -(mu_fn(w)-mu_F)/sqrt(var_fn(w)); 
% construct Sharpe ratio but do not forget minus sign above
% since we use fmincon later
w_tangent = fmincon(obj,w_0,[],[],Aeq_MVP,beq_MVP);
mu_tangent = w_tangent'*mu;
sd_tangent = sqrt(var_fn(w_tangent));

% check with the analytical solution of tangent portfolio
v = Sigma^(-1)*(mu-mu_F*all_one_n);
w_tangent_analytical = v/(all_one_n'*v);
mu_tangent_analytical = w_tangent_analytical'*mu;
sd_tangent_analytical = sqrt(var_fn(w_tangent_analytical));

%% efficient frontier: tangent line

eff = @(x) (mu_tangent-mu_F)/sd_tangent*x+mu_F;
sd_eff_P = xmax;

%% to plot feasible region and efficient frontier
% adjust figure size and location based on your own PC screen
f=figure;
floc = f.Position(1:2);
f.Position(1:2) = [floc(1)/1.1,floc(2)/2];
fsize = f.Position(3:4);
f.Position(3:4) = 1.3*fsize;

% feasible region plot:
for i = 1:num_level
    plot([sd_min_var(i),xmax],[mu_min_var(i),mu_min_var(i)],...
        'LineWidth',3,'color',color_gray);
    hold on;
end

% simulated portfolios plot:
scatter_color=color_purple;%'r';%color_gray
scatter_pt_size = 6;
plot(sd_rand,mu_rand,'.','markersize',scatter_pt_size,...
    'MarkerFaceColor',scatter_color,'MarkerEdgeColor',scatter_color);
hold on;
p(2)=plot(sd_rand(7),mu_rand(7),'.','markersize',scatter_pt_size,...
    'MarkerFaceColor',scatter_color,'MarkerEdgeColor',scatter_color);
%p(5)=scatter(sd_rand(5),mu_rand(5),scatter_pt_size,scatter_color,'filled');
mylgd{2}=['possible portfolios by simulations'];

% constituents plot
for i = 1:n
    plot(sqrt(Sigma(i,i)),mu(i),'.','color','black','markersize',10);
    hold on;
    if i==n
        p(1)=plot(sqrt(Sigma(i,i)),mu(i),'.','color','black','markersize',10);
        %mylgd{1}=['30 contituents of DJIA (estimated from' char(10) 'monthly data, 2020.07~2022.06)'];
        mylgd{1}=['28 contituents of DJIA (estimated from' char(10) 'monthly data 2017.07~2022.06)'];
        hold on;
    end    
end

% minimum variance set plot:
p(3)=plot(sd_min_var,mu_min_var,'LineWidth',2,'color',color_blue);
mylgd{3}=['minimum variance set'];
hold on;

% eff frontier plot:
p(4)=plot([0,sd_eff_P],[mu_F,eff(sd_eff_P)],'LineWidth',2,'color','r');
mylgd{4}=['efficient frontier (with inclusion of' char(10) 'risk-free asset, borrowing allowed)'];%' char(10) '
hold on;

plot(sd_MVP,mu_MVP,'.','color',color_blue,'markersize',15);
%mylgd{4}=['MVP'];
hold on;
plot(0,mu_F,'.','color',color_green,'markersize',15);
%mylgd{5}=['riskless asset'];
hold on;
plot(sd_tangent,mu_tangent,'.','color','r','markersize',15);
%mylgd{6}=['Market portfolio (tangent portfolio)'];
hold on;

%text(sigma_A-0.0005,mu_A+0.002,'A', 'FontSize',12);
%text(sigma_B+0.0001,mu_B+0.002,'B', 'FontSize',12);
text(sd_MVP-0.011,mu_MVP,'MVP', 'FontSize',12);
str={'feasible region', ' (unbounded)'};
text(0.06,0.045,str, 'FontSize',10);
str2={'Market portfolio  ', '(tangent portfolio)'};
annotation('textarrow',[0.365,0.435],[0.77,0.72],'String',...
    str2);
str2={'risk-free asset'};
annotation('textarrow',[0.165,0.135],[0.22,0.24],'String',...
    str2);
%set(gca,'XTick',0.01:0.005:0.025,'xlim',[0.01,0.025],...
%    'ytick',0:0.01:0.06,'ylim',[0,0.06]);
xlabel('\sigma');
ylabel('\mu');
set(gca,'xlim',[0,xmax],...
    'ylim',[zmin,zmax]);
lg=legend(p,mylgd,'location','northeast',...
    'box','off');%'Orientation','Horizontal',

hold off;
