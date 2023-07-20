%% COMPARISON OF COMMUNICATION COST VS EQUILIBRIUM ACCURACY
% Juan Parras, GAPS-UPM, March 2017
clear all; clc; close all;

%% Initial parameters

np=2; %Number of players
na1=2; %Number of actions per player
na2=2;

n_avg=1e2; %Number of averages per equilibrium
n_iter_rm=1e3; %Iterations in RM

l_grid=50;
np_corr=round((6*l_grid^2)^(1/3));  %Grid parameters for correlated
N_com_v=10:10:200;
N_com_v=[5 N_com_v];

delta=0.9; %Test delta

% Plot parameter
pl=1;
% Save parameter
sa=1;

% SOO parameters
kdev=10;

%% Define games to be tested
ng=1; %Number of games to be tested

U=cell(ng,2); %To store player's payoff
static_nash_str=cell(ng,2); %To store player's nash strategies
inc_eq=cell(ng,2);
inc_cor=cell(ng,2);

% % Game 5: Network game
% U{1,1}=[-1 0; 1 -1]; %Payoff to player 1
% U{1,2}=[1 0;-1 0]; %Payoff to player 2
% static_nash_str{1,1}=1/2;
% static_nash_str{1,2}=1/3;
% 
% pareto_v1=linspace(-1/3,0,500);
% pareto_v2=-pareto_v1;
% pareto_m=[pareto_v1' pareto_v2'];

% Game 2: Prisonners dilemma
U{1,1}=[2 -1; 3 0];
U{1,2}=[2 3; -1 0];
static_nash_str{1,1}=0;
static_nash_str{1,2}=0;
% Pareto region
pareto_v1a=linspace(0,2,100);
pareto_v2a=(-pareto_v1a+8)/3;
pareto_v1b=linspace(2,8/3,100);
pareto_v2b=-3*pareto_v1b+8;
pareto_v1=[pareto_v1a pareto_v1b];
pareto_v2=[pareto_v2a pareto_v2b];
pareto_m=[pareto_v1' pareto_v2'];
%% Test games
% Cases
% 1 - Regret matching
% 2 - Nash, equisp,
% 3 - Nash, unif,
% 4 - Nash, soo, lambda=0.5
% 5 - Nash, soo, lambda =0.75
% 6 - Nash, soo, lambda = 1
% 7 - Corr, equisp,
% 8 - Corr, unif,
% 9 - Corr, soo, lambda=0.5
% 10 - Corr, soo, lambda =0.75
% 11 - Corr, soo, lambda = 1

dist=zeros(length(N_com_v),n_avg,11,ng);

for game=1:ng
    for i=1:length(N_com_v)
        N_com=N_com_v(i);

        display(['Game ' num2str(game) ' of ' num2str(ng) ', N_com number = ' num2str(i) ' of ' num2str(length(N_com_v))]);

        u_1=U{game,1}; %Payoff to player 1
        u_2=U{game,2}; %Payoff to player 2

        g1=@(y,z) y.*z.*u_1(1,1)+y.*(1-z).*u_1(1,2)+(1-y).*z.*u_1(2,1)+(1-y).*(1-z).*u_1(2,2); % Payoff function for player 1
        g2=@(y,z) y.*z.*u_2(1,1)+y.*(1-z).*u_2(1,2)+(1-y).*z.*u_2(2,1)+(1-y).*(1-z).*u_2(2,2); % Payoff function for player 2

        yn=static_nash_str{game,1};
        zn=static_nash_str{game,2};

        v1n=g1(yn,zn);
        v2n=g2(yn,zn);
        for j=1:n_avg
            % Obtain RM strategy
            [~,~,a1rm,a2rm]=regret_min(na1,na2,u_1,u_2,n_iter_rm);
            a1rm=sum(a1rm(:,1))/n_iter_rm;
            a2rm=sum(a2rm(:,1))/n_iter_rm;
            % In the incomings, use RM as threat
            lambda=0.5;
            % Case 1: sample using equidistant
            [a1eq,a2eq,~]=obtain_CA_nash(N_com,l_grid,delta,g1,g2,a1rm, a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,1);
            % Case 2: sample using uniform distribution
            [a1eq_in,a2eq_in,~]=obtain_CA_nash(N_com,l_grid,delta,g1,g2,a1rm, a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,2);
            % Case 3: sample using SOO (lambda=0)
            [a1eq_soo1,a2eq_soo1,~]=obtain_CA_nash(N_com,l_grid,delta,g1,g2,a1rm, a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,3);
            % Case 3: sample using SOO
            lambda=0.75;
            [a1eq_soo2,a2eq_soo2,~]=obtain_CA_nash(N_com,l_grid,delta,g1,g2,a1rm, a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,3);
            % Case 3: sample using SOO
            lambda=1;
            [a1eq_soo3,a2eq_soo3,~]=obtain_CA_nash(N_com,l_grid,delta,g1,g2,a1rm, a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,3);

            % Case 1 corr: sample using equidistant
            lambda=0.5;
            [phi,v_corr,~]=obtain_CA_corr(N_com,np_corr,delta,u_1,u_2,a1rm,a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,1);
            % Case 2 corr: sample using random
            [phi_in,v_corr_in,~]=obtain_CA_corr(N_com,np_corr,delta,u_1,u_2,a1rm,a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,2);
            % Case 3 corr: sample using SOO
            [phi_soo,v_corr_soo1,~]=obtain_CA_corr(N_com,np_corr,delta,u_1,u_2,a1rm,a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,3);
            % Case 3 corr: sample using SOO
            lambda=0.75;
            [phi_soo,v_corr_soo2,~]=obtain_CA_corr(N_com,np_corr,delta,u_1,u_2,a1rm,a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,3);
            % Case 3 corr: sample using SOO
            lambda=1;
            [phi_soo,v_corr_soo3,~]=obtain_CA_corr(N_com,np_corr,delta,u_1,u_2,a1rm,a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,3);

            % Save distances
            dist(i,j,1,game)=min(pdist2([g1(a1rm,a2rm) g2(a1rm,a2rm)],pareto_m));
            dist(i,j,2,game)=min(pdist2([g1(a1eq,a2eq) g2(a1eq,a2eq)],pareto_m));
            dist(i,j,3,game)=min(pdist2([g1(a1eq_in,a2eq_in) g2(a1eq_in,a2eq_in)],pareto_m));
            dist(i,j,4,game)=min(pdist2([g1(a1eq_soo1,a2eq_soo1) g2(a1eq_soo1,a2eq_soo1)],pareto_m));
            dist(i,j,5,game)=min(pdist2([g1(a1eq_soo2,a2eq_soo2) g2(a1eq_soo2,a2eq_soo2)],pareto_m));
            dist(i,j,6,game)=min(pdist2([g1(a1eq_soo3,a2eq_soo3) g2(a1eq_soo3,a2eq_soo3)],pareto_m));
            dist(i,j,7,game)=min(pdist2([v_corr(1) v_corr(2)],pareto_m));
            dist(i,j,8,game)=min(pdist2([v_corr_in(1) v_corr_in(2)],pareto_m));
            dist(i,j,9,game)=min(pdist2([v_corr_soo1(1) v_corr_soo1(2)],pareto_m));
            dist(i,j,10,game)=min(pdist2([v_corr_soo2(1) v_corr_soo2(2)],pareto_m));
            dist(i,j,11,game)=min(pdist2([v_corr_soo3(1) v_corr_soo3(2)],pareto_m));

        end
    end
end


%% Plotting section
if pl
    figure();
    semilogy(N_com_v,mean(dist(:,:,1)'),'c*-');
    hold on;
    grid on;
    plot(N_com_v,mean(dist(:,:,2)'),'r-','LineWidth',2);
    plot(N_com_v,mean(dist(:,:,3)'),'k--','LineWidth',2);
    plot(N_com_v,mean(dist(:,:,4)'),'b-.','LineWidth',2);
    %plot(N_com_v,mean(dist(:,:,5)'),'gs-');
    plot(N_com_v,mean(dist(:,:,6)'),'m:','LineWidth',2);
    plot(N_com_v,mean(dist(:,:,7)'),'r-','LineWidth',3);
    plot(N_com_v,mean(dist(:,:,8)'),'k--','LineWidth',3);
    plot(N_com_v,mean(dist(:,:,9)'),'b-.','LineWidth',3);
    %plot(N_com_v,mean(dist(:,:,10)'),'gs:');
    plot(N_com_v,mean(dist(:,:,11)'),'m:','LineWidth',3);
    xlabel('N_{c}');
    ylabel('\xi');
    legend('RM','Nash, eq','Nash, unif','Nash, SOO, \lambda=0.5','Nash, SOO, \lambda=1','Corr, eq','Corr, unif','Corr, SOO, \lambda=0.5','Corr, SOO, \lambda=1');

    figure();
    semilogy(N_com_v,max(dist(:,:,1)'),'c*-');
    hold on;
    grid on;
    plot(N_com_v,max(dist(:,:,2)'),'r-','LineWidth',2);
    plot(N_com_v,max(dist(:,:,3)'),'k--','LineWidth',2);
    plot(N_com_v,max(dist(:,:,4)'),'b-.','LineWidth',2);
   % plot(N_com_v,max(dist(:,:,5)'),'gs-');
    plot(N_com_v,max(dist(:,:,6)'),'m:','LineWidth',2);
    plot(N_com_v,max(dist(:,:,7)'),'r-','LineWidth',3);
    plot(N_com_v,max(dist(:,:,8)'),'k--','LineWidth',3);
    plot(N_com_v,max(dist(:,:,9)'),'b-.','LineWidth',3);
    %plot(N_com_v,max(dist(:,:,10)'),'gs:');
    plot(N_com_v,max(dist(:,:,11)'),'m:','LineWidth',3);
    xlabel('N_{c}');
    ylabel('\xi');
    legend('RM','Nash, eq','Nash, unif','Nash, SOO, \lambda=0.5','Nash, SOO, \lambda=1','Corr, eq','Corr, unif','Corr, SOO, \lambda=0.5','Corr, SOO, \lambda=1');

end
%% Saving section
if sa
    save('Data_nash_learning_communication_testbench_cost_vs_accuracy');
end