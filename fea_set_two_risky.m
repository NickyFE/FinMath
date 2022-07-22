% FinMath Series Topic1 Feasible set between two assets
% 
% this file serves as a simple example to construct all possible
% portfolios between two given risky assets A and B,
% namely, the feasible set between them

% Copyright Nicky¿œ ¶ 2022
% nickycanhelp@qq.com
% Github: NickyFE/FinMath

clear;

%% assets settings
mu_A = 0.05;
sigma_A = 0.02;
mu_B = 0.02;
sigma_B = 0.015;
rho = 0.15;

%% compute weight of A in (global) MVP

MVP_w_A = (sigma_B^2-rho*sigma_A*sigma_B)/...
    (sigma_A^2+sigma_B^2-2*rho*sigma_A*sigma_B);

%% define the variance of the portfolio as the function of w_A
var_P = @(w_A) (sigma_A^2+sigma_B^2-2*rho*sigma_A*sigma_B)*...
    ((w_A*mu_A+(1-w_A)*mu_B-mu_B)/(mu_A-mu_B))^2+...
    2*(rho*sigma_A*sigma_B-sigma_B^2)*...
    (w_A*mu_A+(1-w_A)*mu_B-mu_B)/(mu_A-mu_B)+sigma_B^2;

%% compute mu and sigma of MVP

mu_MVP = MVP_w_A*mu_A+(1-MVP_w_A)*mu_B;
sigma_MVP = sqrt(var_P(MVP_w_A));

%% compute mu and sigma of any portfolio between A and B
% by varying w_A
% suppose shorting is allowed

w_A_list = -0.3:0.0001:1.2;
sigma_P_list = zeros(1,length(w_A_list));
mu_P_list = zeros(1,length(w_A_list));
for i=1:length(w_A_list)
    sigma_P_list(i)=sqrt(var_P(w_A_list(i)));
    mu_P_list(i)=w_A_list(i)*mu_A+(1-w_A_list(i))*mu_B;
end

%% plot the feasible curve 
figure;
plot(sigma_P_list,mu_P_list,'LineWidth',2);
hold on;
plot(sigma_A,mu_A,'.','color','black','markersize',15);
hold on;
plot(sigma_B,mu_B,'.','color','black','markersize',15);
hold on;
plot(sigma_MVP,mu_MVP,'.','color','black','markersize',15);
text(sigma_A-0.0005,mu_A+0.002,'A', 'FontSize',12);
text(sigma_B+0.0001,mu_B+0.002,'B', 'FontSize',12);
text(sigma_MVP-0.0014,mu_MVP+0.001,'MVP', 'FontSize',12);
set(gca,'XTick',0.01:0.005:0.025,'xlim',[0.01,0.025],...
    'ytick',0:0.01:0.06,'ylim',[0,0.06]);
xlabel('\sigma');
ylabel('\mu');

hold off;
