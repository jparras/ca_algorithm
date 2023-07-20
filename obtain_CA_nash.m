function[a1,a2,v_out]=obtain_CA_nash(N_com,l_grid,delta,g1,g2,threat_str_1, threat_str_2,threat_v1,threat_v2,lambda,sampling_method)

% Declare variables

Nmax_iter=1e3; %Max number of tryouts in each point search - to avoid stalling
valid_grid=[];
ac_dev=linspace(0,1,l_grid);
v1c=[];
v2c=[];
v_out=[];
global param; % Declared as global
param=struct('lambda',lambda,'valid_grid',valid_grid,'g1',g1,'g2',g2,'ac_dev',ac_dev,'v1c',v1c,'v2c',v2c,'threat_v1',threat_v1,'threat_v2',threat_v2,'delta',delta);

% Sampling method: 1 equidistant, 2 random uniform, 3 SOO
switch sampling_method
    case 1
        [valid_grid,v_out]=sample_equidist_nash(N_com,l_grid,delta,g1,g2,threat_v1,threat_v2,threat_str_1, threat_str_2);
    case 2
        [valid_grid,v_out]=sample_random_nash(N_com,Nmax_iter,delta,g1,g2,threat_v1,threat_v2,ac_dev,threat_str_1, threat_str_2);
    case 3
        % Use SOO to sample

        settings.dim = 2;
        settings.type = 'det';
        [jitter]=oo(@function_cost_soo_p1,N_com,settings);
        [jitter]=oo(@function_cost_soo_p2,N_com,settings);

        % Update global variables
        par=getGlobalx;
        v1c=par.v1c;
        v2c=par.v2c;
        valid_grid=par.valid_grid;
        %Region obtention

        v1f=zeros(size(valid_grid,1),1);
        v2f=zeros(size(valid_grid,1),1);
        for i=1:size(valid_grid,1)
            v1f(i)=g1(valid_grid(i,1),valid_grid(i,2));
            v2f(i)=g2(valid_grid(i,1),valid_grid(i,2));
        end
        v1t=g1(threat_str_1,threat_str_2);
        v2t=g2(threat_str_1,threat_str_2);
        v_out=cell(6,1);
        v_out{1}=v1c;
        v_out{2}=v2c;
        v_out{3}=v1f;
        v_out{4}=v2f;
        v_out{5}=v1t;
        v_out{6}=v2t;
end

if ~isempty(valid_grid)
    %plot(g1(valid_grid(:,1),valid_grid(:,2)),g2(valid_grid(:,1),valid_grid(:,2)),'o')
    % Each player sorts its payoffs
    % Player 1
    y_grid=valid_grid;
    for i=1:size(y_grid,1);
        y_grid(i,3)=g1(y_grid(i,1),y_grid(i,2));
    end
    y_grid=sortrows(y_grid,-3); %Sort by payoff value
    % Player 2
    z_grid=valid_grid;
    for i=1:size(z_grid,1);
        z_grid(i,3)=g2(z_grid(i,1),z_grid(i,2));
    end
    z_grid=sortrows(z_grid,-3); %Sort by payoff value


    % Find a pareto value
    pareto_search=1;
    a1=threat_str_1; %Security intialization: update if possible
    a2=threat_str_2;

    while pareto_search
        %Jointly controlled lottery
        w1=rand(1);
        w2=rand(1);
        w=w1+w2;
        if w>1
            w=w-1;
        end
        id=round(1+(size(valid_grid,1)-1)*w);
        st=valid_grid(id,:); % Starting point for search
        a1=st(1);
        a2=st(2);
        %Player 1 erases dominated strategies
        id1=find(y_grid(:,3)==g1(a1,a2),1);
        y_grid=y_grid(1:id1,:); %Update y_grid
        %Player 2 erases dominated strategies
        id2=find(z_grid(:,3)==g2(a1,a2),1);
        z_grid=z_grid(1:id2,:); %Update z_grid
        % Intersect and reorder
        valid_grid=intersect(y_grid(:,1:2),z_grid(:,1:2),'rows');
        if isempty(valid_grid) || size(valid_grid,1)==1 % Dominating strategy found
            pareto_search=0;
        else
            % Each player sorts its payoffs
            % Player 1
            y_grid=valid_grid;
            for i=1:size(y_grid,1);
                y_grid(i,3)=g1(y_grid(i,1),y_grid(i,2));
            end
            y_grid=sortrows(y_grid,-3); %Sort by payoff value
            % Player 2
            z_grid=valid_grid;
            for i=1:size(z_grid,1);
                z_grid(i,3)=g2(z_grid(i,1),z_grid(i,2));
            end
            z_grid=sortrows(z_grid,-3); %Sort by payoff value
        end

    end
else %No valid equilibrium for both players
    a1=threat_str_1;
    a2=threat_str_2;
end