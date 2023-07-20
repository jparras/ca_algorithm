function[valid_grid,v_out]=sample_equidist_corr(N_com,n_points,delta,g1,g2,u_1,u_2,threat_str_1, threat_str_2,threat_v1,threat_v2);

% Obtain sample grid of valid equilibria: 
phi_grid=[];
% for i=1:n_points %This idea can be found in Donald B. Rubin, The Bayesian bootstrap Ann. Statist. 9, 1981, 130-134
%     r_v=rand(1,3);
%     aux=[0 sort(r_v) 1];
%     phi_grid(i,:)=diff(aux);
% end
aux=linspace(0,1,n_points);
for i1=1:n_points
    for i2=1:n_points
        for i3=1:n_points
            va=aux(i1)+aux(i2)+aux(i3);
            if va<=1
                phi_grid(end+1,:)=[aux(i1) aux(i2) aux(i3) 1-va];
            end
        end
    end
end

% Player 1
p1_grid=[];
for i=1:size(phi_grid,1)
    phi_test=phi_grid(i,:);
    v1o=g1(phi_test);
    restr_1=(1-delta)*(phi_test(1)*(u_1(2,1)-u_1(1,1))+phi_test(2)*(u_1(2,2)-u_1(1,2)))+delta*(threat_v1-v1o)*(phi_test(1)+phi_test(2));
    restr_2=(1-delta)*(phi_test(3)*(u_1(1,1)-u_1(2,1))+phi_test(4)*(u_1(1,2)-u_1(2,2)))+delta*(threat_v1-v1o)*(phi_test(3)+phi_test(4));
    if restr_1<=0 && restr_2<=0 %No profitable deviation
        if v1o>=threat_v1 % Better eq than threat
            p1_grid(end+1,:)=phi_test;
        end
    end
end

% Player 1
p2_grid=[];
for i=1:size(phi_grid,1)
    phi_test=phi_grid(i,:);
    v2o=g2(phi_test);
    restr_1=(1-delta)*(phi_test(1)*(u_2(1,2)-u_2(1,1))+phi_test(3)*(u_2(2,2)-u_2(2,1)))+delta*(threat_v2-v2o)*(phi_test(1)+phi_test(3));
    restr_2=(1-delta)*(phi_test(2)*(u_2(1,1)-u_2(1,2))+phi_test(4)*(u_2(2,1)-u_2(2,2)))+delta*(threat_v2-v2o)*(phi_test(2)+phi_test(4));
    if restr_1<=0 && restr_2<=0 %No profitable deviation
        if v2o>=threat_v2 % Better eq than threat
            p2_grid(end+1,:)=phi_test;
        end
    end
end

% Randomly take N_com values from y_grid, z_grid only
if size(p1_grid,1)>N_com
    p1_grid=p1_grid(randperm(size(p1_grid,1)),:);
    p1_grid=p1_grid(1:N_com,:);
end
if size(p2_grid,1)>N_com
    p2_grid=p2_grid(randperm(size(p2_grid,1)),:);
    p2_grid=p2_grid(1:N_com,:);
end

% Refine valid points: only those in which both players gain
valid_grid=intersect(p1_grid,p2_grid,'rows');

% Obtain region points
v1c=zeros(size(phi_grid,1),1);
v2c=zeros(size(phi_grid,1),1);
for i=1:size(phi_grid,1)
    v1c(i)=g1(phi_grid(i,:));
    v2c(i)=g2(phi_grid(i,:));
end
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