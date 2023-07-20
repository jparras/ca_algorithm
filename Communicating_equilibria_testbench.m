%% COMMUNICATING EQUILIBRIA TESTBENCH
% Juan Parras, GAPS-UPM, February 2017
clear all; clc; close all;

%% Initial parameters

np=2; %Number of players
na1=2; %Number of actions per player
na2=2;

n_avg=100; %Number of averages per equilibrium
n_iter_rm=1e3; %Iterations in RM

l_grid=50; %Grid parameters for nash spe
np_corr=round((6*l_grid^2)^(1/3));  %Grid parameters for correlated
n_com=100; %Number of communications

delta_v=[0.1 0.5 0.9]; %For testing

% Plot parameter
pl=0;
% Save parameter
sa=1;
%% Prepare storing values
ng=5; %Number of games to be tested

U=cell(ng,2); %To store player's payoff
static_nash_str=cell(ng,2); %To store player's nash strategies
Nc=9; %Number of cases
v_out=cell(ng,length(delta_v),2,Nc);
inc=cell(ng,length(delta_v),2,Nc); %9 cases as:
% 1 - Nash equispaced
% 2 - Nash random
% 3 - Nash SOO, lambda=0.5
% 4 - Nash SOO, lambda = 1
% 5 - Corr equispaced
% 6 - Corr random
% 7 - Corr SOO, lambda=0.5
% 8 - Corr SOO, lambda = 1
% 9 - RM

reg_plot=cell(ng,length(delta_v),10,Nc-1); %(v1region,v2region,v1feasible,v2feasible, v1n,v2n,v1com,v2com,v1rm,v2rm)
% All cases except for RM

%% Prepare game values
% Game 1: Matching pennies
U{1,1}=[1 -1; -1 1]; %Payoff to player 1
U{1,2}=[-1 1; 1 -1]; %Payoff to player 2
static_nash_str{1,1}=1/2;
static_nash_str{1,2}=1/2;

% Game 2: Prisonners dilemma
U{2,1}=[2 -1; 3 0];
U{2,2}=[2 3; -1 0];
static_nash_str{2,1}=0;
static_nash_str{2,2}=0;

% Game 3: Battle of sexes
U{3,1}=[2 0; 0 1]; %Payoff to player 1
U{3,2}=[1 0; 0 2]; %Payoff to player 2
static_nash_str{3,1}=2/3;
static_nash_str{3,2}=1/3;

% Game 4: Chicken game
U{4,1}=[-10 1; -1 0]; %Payoff to player 1
U{4,2}=[-10 -1; 1 0];
static_nash_str{4,1}=1/10;
static_nash_str{4,2}=1/10;

% Game 5: Network game
U{5,1}=[-1 0; 1 -1]; %Payoff to player 1
U{5,2}=[1 0;-1 0]; %Payoff to player 2
static_nash_str{5,1}=1/2;
static_nash_str{5,2}=1/3;

%% Test games

for game=1:ng;
    inc_aux=zeros(length(delta_v),n_avg,2,Nc); % (kd,j, player,case)
    v_aux=zeros(length(delta_v),n_avg,2,Nc); % (kd,j, player,case)
    r_aux=cell(8,1);
    
    % Load game parameters
    u_1=U{game,1}; %Payoff to player 1
    u_2=U{game,2}; %Payoff to player 2

    g1=@(y,z) y.*z.*u_1(1,1)+y.*(1-z).*u_1(1,2)+(1-y).*z.*u_1(2,1)+(1-y).*(1-z).*u_1(2,2); % Payoff function for player 1
    g2=@(y,z) y.*z.*u_2(1,1)+y.*(1-z).*u_2(1,2)+(1-y).*z.*u_2(2,1)+(1-y).*(1-z).*u_2(2,2); % Payoff function for player 2

    yn=static_nash_str{game,1};
    zn=static_nash_str{game,2};

    v1n=g1(yn,zn);
    v2n=g2(yn,zn);
    
    for kd=1:length(delta_v)
        delta=delta_v(kd);
        display(['Game ' num2str(game) ' of ' num2str(ng) ', delta = ' num2str(delta)]);
        for j=1:n_avg
            % Case 9: RM
            [~,~,a1rm,a2rm]=regret_min(na1,na2,u_1,u_2,n_iter_rm);
            a1rm=sum(a1rm(:,1))/n_iter_rm;
            a2rm=sum(a2rm(:,1))/n_iter_rm;
            inc_aux(kd,j,1,9)=g1(a1rm,a2rm)-g1(a1rm,a2rm); %Always 0
            inc_aux(kd,j,2,9)=g2(a1rm,a2rm)-g2(a1rm,a2rm);
            v_aux(kd,j,1,9)=g1(a1rm,a2rm);
            v_aux(kd,j,2,9)=g2(a1rm,a2rm);
            % In the incomings, use RM as threat
            kdev=10;
            lambda=0.5;
            % Case 1: Nash equispaced
            [a1_1,a2_1,r_aux{1}]=obtain_CA_nash(n_com,l_grid,delta,g1,g2,a1rm, a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,1);
            inc_aux(kd,j,1,1)=g1(a1_1,a2_1)-g1(a1rm,a2rm);
            inc_aux(kd,j,2,1)=g2(a1_1,a2_1)-g2(a1rm,a2rm);
            v_aux(kd,j,1,1)=g1(a1_1,a2_1);
            v_aux(kd,j,2,1)=g2(a1_1,a2_1);
            % Case 2: Nash uniform
            [a1_2,a2_2,r_aux{2}]=obtain_CA_nash(n_com,l_grid,delta,g1,g2,a1rm, a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,2);
            inc_aux(kd,j,1,2)=g1(a1_2,a2_2)-g1(a1rm,a2rm);
            inc_aux(kd,j,2,2)=g2(a1_2,a2_2)-g2(a1rm,a2rm);
            v_aux(kd,j,1,2)=g1(a1_2,a2_2);
            v_aux(kd,j,2,2)=g2(a1_2,a2_2);
            % Case 3: Nash SOO lambda=0.5
            lambda=0.5;
            [a1_3,a2_3,r_aux{3}]=obtain_CA_nash(n_com,l_grid,delta,g1,g2,a1rm, a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,3);
            inc_aux(kd,j,1,3)=g1(a1_3,a2_3)-g1(a1rm,a2rm);
            inc_aux(kd,j,2,3)=g2(a1_3,a2_3)-g2(a1rm,a2rm);
            v_aux(kd,j,1,3)=g1(a1_3,a2_3);
            v_aux(kd,j,2,3)=g2(a1_3,a2_3);
            % Case 4: Nash SOO, lambda = 1
            lambda=1;
            [a1_4,a2_4,r_aux{4}]=obtain_CA_nash(n_com,l_grid,delta,g1,g2,a1rm, a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,3);
            inc_aux(kd,j,1,4)=g1(a1_4,a2_4)-g1(a1rm,a2rm);
            inc_aux(kd,j,2,4)=g2(a1_4,a2_4)-g2(a1rm,a2rm);
            v_aux(kd,j,1,4)=g1(a1_4,a2_4);
            v_aux(kd,j,2,4)=g2(a1_4,a2_4);
            % Case 5: Corr equispaced
            [phi_5,v_corr_5,r_aux{5}]=obtain_CA_corr(n_com,np_corr,delta,u_1,u_2,a1rm,a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,1);
            inc_aux(kd,j,1,5)=v_corr_5(1)-g1(a1rm,a2rm);
            inc_aux(kd,j,2,5)=v_corr_5(2)-g2(a1rm,a2rm);
            v_aux(kd,j,1,5)=v_corr_5(1);
            v_aux(kd,j,2,5)=v_corr_5(2);
            % Case 6: Corr uniform
            [phi_6,v_corr_6,r_aux{6}]=obtain_CA_corr(n_com,np_corr,delta,u_1,u_2,a1rm,a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,2);
            inc_aux(kd,j,1,6)=v_corr_6(1)-g1(a1rm,a2rm);
            inc_aux(kd,j,2,6)=v_corr_6(2)-g2(a1rm,a2rm);
            v_aux(kd,j,1,6)=v_corr_6(1);
            v_aux(kd,j,2,6)=v_corr_6(2);
            % Case 7: Corr SOO lambda = 0.5
            lambda=0.5;
            [phi_7,v_corr_7,r_aux{7}]=obtain_CA_corr(n_com,np_corr,delta,u_1,u_2,a1rm,a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,3);
            inc_aux(kd,j,1,7)=v_corr_7(1)-g1(a1rm,a2rm);
            inc_aux(kd,j,2,7)=v_corr_7(2)-g2(a1rm,a2rm);
            v_aux(kd,j,1,7)=v_corr_7(1);
            v_aux(kd,j,2,7)=v_corr_7(2);
            % Case 8: Corr SOO lambda = 1
            lambda=1;
            [phi_8,v_corr_8,r_aux{8}]=obtain_CA_corr(n_com,np_corr,delta,u_1,u_2,a1rm,a2rm,g1(a1rm,a2rm),g2(a1rm,a2rm),lambda,3);
            inc_aux(kd,j,1,8)=v_corr_8(1)-g1(a1rm,a2rm);
            inc_aux(kd,j,2,8)=v_corr_8(2)-g2(a1rm,a2rm);
            v_aux(kd,j,1,8)=v_corr_8(1);
            v_aux(kd,j,2,8)=v_corr_8(2);
        end
        % Save region plots
        for c=1:Nc-1
           reg_plot{game,kd,1,c}=r_aux{c}{1};
           reg_plot{game,kd,2,c}=r_aux{c}{2};
           reg_plot{game,kd,3,c}=r_aux{c}{3};
           reg_plot{game,kd,4,c}=r_aux{c}{4};
           reg_plot{game,kd,5,c}=r_aux{c}{5};
           reg_plot{game,kd,6,c}=r_aux{c}{6};
           reg_plot{game,kd,7,c}=v_aux(kd,:,1,c);
           reg_plot{game,kd,8,c}=v_aux(kd,:,2,c);
           reg_plot{game,kd,9,c}=v_aux(kd,:,1,9);
           reg_plot{game,kd,10,c}=v_aux(kd,:,2,9);
        end
        % Save increments
        for c=1:Nc
            inc{game,kd,1,c}=inc_aux(kd,:,1,c);
            inc{game,kd,2,c}=inc_aux(kd,:,2,c);
            v_out{game,kd,1,c}=v_aux(kd,:,1,c);
            v_out{game,kd,2,c}=v_aux(kd,:,2,c);
        end
    end
end


%% Plotting section
if pl
    
    for i=1:ng
        y1=zeros(Nc-1,length(delta_v));
        y2=zeros(Nc-1,length(delta_v));
        for j=1:length(delta_v)
            for k=1:Nc-1
                y1(k,j)=mean(inc{i,j,1,k});
                y2(k,j)=mean(inc{i,j,2,k});
            end
        end
%         figure()
%         bar([y1;y2]);
%         grid on;
%         legend(['\delta=' num2str(delta_v(1))],['\delta=' num2str(delta_v(2))],['\delta=' num2str(delta_v(3))])
        figure();
        plot(1:Nc-1,y1(:,1),'bx-',1:Nc-1,y1(:,2),'rx-',1:Nc-1,y1(:,3),'kx-',1:Nc-1,y2(:,1),'b--o',1:Nc-1,y2(:,2),'r--o',1:Nc-1,y2(:,3),'k--o');
        grid on;
        legend(['P1, \delta=' num2str(delta_v(1))],['P1, \delta=' num2str(delta_v(2))],['P1, \delta=' num2str(delta_v(3))],['P2, \delta=' num2str(delta_v(1))],['P2, \delta=' num2str(delta_v(2))],['P2, \delta=' num2str(delta_v(3))])

    end

    for pic=1:Nc-1
        figure();
        id=0;
        for i=1:ng
            id=i;
            for j=1:length(delta_v)
                subplot(length(delta_v),ng,id);
                plot(reg_plot{i,j,1,pic},reg_plot{i,j,2,pic},'.c'); %Possible region
                grid on;
                hold on;
                plot(reg_plot{i,j,3,pic},reg_plot{i,j,4,pic},'.k'); %Feasible points
                plot(reg_plot{i,j,9,pic},reg_plot{i,j,10,pic},'sg'); %RM eq
                plot(reg_plot{i,j,7,pic},reg_plot{i,j,8,pic},'xb'); %Comm eq nash
                plot(reg_plot{i,j,5,pic},reg_plot{i,j,6,pic},'or'); %Static nash eq
                title(['Game ' num2str(i) ', delta = ' num2str(delta_v(j)) ', case = ' num2str(pic)]);
                id=id+ng;
            end
        end
    end
    
    
end
%% Error computation for Nash
max_1=zeros(1,ng);
max_2=zeros(1,ng);
e_1=zeros(ng,length(delta_v));
e_2=zeros(ng,length(delta_v));
for i=1:ng
    u_1=U{i,1}; %Payoff to player 1
    u_2=U{i,2}; %Payoff to player 2
    grad_1=@(x) -sqrt(((u_1(1,1)-u_1(1,2)-u_1(2,1)+u_1(2,2)).*x(2)+(u_1(1,2)-u_1(2,2))).^2+((u_1(1,1)-u_1(1,2)-u_1(2,1)+u_1(2,2)).*x(1)+(u_1(2,1)-u_1(2,2))).^2);
    grad_2=@(x) -sqrt(((u_2(1,1)-u_2(1,2)-u_2(2,1)+u_2(2,2)).*x(2)+(u_2(1,2)-u_2(2,2))).^2+((u_2(1,1)-u_2(1,2)-u_2(2,1)+u_2(2,2)).*x(1)+(u_2(2,1)-u_2(2,2))).^2);
    [p1, max_1(i)]=fmincon(grad_1,[0 0],[],[],[],[],[0,0],[1 1]);
    [p2, max_2(i)]=fmincon(grad_2,[0 0],[],[],[],[],[0,0],[1 1]);
    for j=1:length(delta_v)
        e_1(i,j)=-(1-delta_v(j))*max_1(i)/l_grid;
        e_2(i,j)=-(1-delta_v(j))*max_2(i)/l_grid;
    end
end
e_1
e_2
%% Saving section
if sa
    save('Data_nash_learning_communication_testbench');
end