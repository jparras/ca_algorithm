function[phi,v_corr,v_out]=obtain_CA_corr(N_com,n_points,delta,u_1,u_2,threat_str_1, threat_str_2,threat_v1,threat_v2,lambda,sampling_method)

%Define payoff functions
g1=@(x) u_1(1,1)*x(1)+u_1(1,2)*x(2)+u_1(2,1)*x(3)+u_1(2,2)*x(4);
g2=@(x) u_2(1,1)*x(1)+u_2(1,2)*x(2)+u_2(2,1)*x(3)+u_2(2,2)*x(4);

% phi(1,1)=>1
% phi(1,2)=>2
% phi(2,1)=>3
% phi(2,2)=>4

% Declare variables

Nmax_iter=1e3; %Max number of iterations in each point search - to avoid stalling
valid_grid=[];
v1c=[];
v2c=[];
v_out=[];
global param; % Declared as global
param=struct('lambda',lambda,'valid_grid',valid_grid,'g1',g1,'g2',g2,'v1c',v1c,'v2c',v2c,'threat_v1',threat_v1,'threat_v2',threat_v2,'delta',delta,'u_1',u_1,'u_2',u_2);

% Sampling method: 1 equidistant, 2 random uniform, 3 SOO
switch sampling_method
    case 1
        [valid_grid,v_out]=sample_equidist_corr(N_com,n_points,delta,g1,g2,u_1,u_2,threat_str_1, threat_str_2,threat_v1,threat_v2);
    case 2
        [valid_grid,v_out]=sample_random_corr(N_com,Nmax_iter,delta,g1,g2,u_1,u_2,threat_str_1, threat_str_2,threat_v1,threat_v2);
    case 3
        % Use SOO to sample

        settings.dim = 3; % Sample in a R^3: simplex is dim(3) (although not all points are valid!)
        settings.type = 'det';
        [jitter]=oo(@function_cost_soo_corr_p1,N_com,settings);
        [jitter]=oo(@function_cost_soo_corr_p2,N_com,settings);

        % Update global variables
        par=getGlobalx;
        v1c=par.v1c;
        v2c=par.v2c;
        valid_grid=par.valid_grid;

        % Update region

        v1f=zeros(size(valid_grid,1),1);
        v2f=zeros(size(valid_grid,1),1);
        for i=1:size(valid_grid,1)
            v1f(i)=g1(valid_grid(i,:));
            v2f(i)=g2(valid_grid(i,:));
        end
        phit=[threat_str_1*threat_str_2,threat_str_1*(1-threat_str_2),(1-threat_str_1)*threat_str_2,(1-threat_str_1)*(1-threat_str_2)];
        v1t=g1(phit);
        v2t=g2(phit);
        v_out=cell(6,1);
        v_out{1}=v1c;
        v_out{2}=v2c;
        v_out{3}=v1f;
        v_out{4}=v2f;
        v_out{5}=v1t;
        v_out{6}=v2t;
end

if ~isempty(valid_grid)
    % Each player sorts its payoffs
    % Player 1
    y_grid=valid_grid;
    for i=1:size(y_grid,1);
        y_grid(i,5)=g1(y_grid(i,1:4));
    end
    y_grid=sortrows(y_grid,-5); %Sort by payoff value
    % Player 2
    z_grid=valid_grid;
    for i=1:size(z_grid,1);
        z_grid(i,5)=g2(z_grid(i,1:4));
    end
    z_grid=sortrows(z_grid,-5); %Sort by payoff value


    % Find a pareto value
    pareto_search=1;
    phi=[threat_str_1*threat_str_2,threat_str_1*(1-threat_str_2),(1-threat_str_1)*threat_str_2,(1-threat_str_1)*(1-threat_str_2)];

    while pareto_search
        %Jointly controlled lottery
        w1=rand(1);
        w2=rand(1);
        w=w1+w2;
        if w>1
            w=w-1;
        end
        id=round(1+(size(valid_grid,1)-1)*w);
        phi=valid_grid(id,:); % Starting point for search
        %Player 1 erases dominated strategies
        id1=find(y_grid(:,5)==g1(phi));
        y_grid=y_grid(1:id1,:); %Update y_grid
        %Player 2 erases dominated strategies
        id2=find(z_grid(:,5)==g2(phi));
        z_grid=z_grid(1:id2,:); %Update z_grid
        % Intersect and reorder
        valid_grid=intersect(y_grid(:,1:4),z_grid(:,1:4),'rows');
        if isempty(valid_grid) || size(valid_grid,1)==1 % Dominating strategy found
            pareto_search=0;
        else
            % Each player sorts its payoffs
            % Player 1
            y_grid=valid_grid;
            for i=1:size(y_grid,1);
                y_grid(i,5)=g1(y_grid(i,1:4));
            end
            y_grid=sortrows(y_grid,-5); %Sort by payoff value
            % Player 2
            z_grid=valid_grid;
            for i=1:size(z_grid,1);
                z_grid(i,5)=g2(z_grid(i,1:4));
            end
            z_grid=sortrows(z_grid,-5); %Sort by payoff value
        end

    end
else %No valid equilibrium for both players
    phi=[threat_str_1*threat_str_2,threat_str_1*(1-threat_str_2),(1-threat_str_1)*threat_str_2,(1-threat_str_1)*(1-threat_str_2)];
end
v_corr=[g1(phi),g2(phi)];