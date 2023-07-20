function[valid_grid,v_out]=sample_equidist_nash(N_com,l_grid,delta,g1,g2,threat_v1,threat_v2,threat_str_1, threat_str_2);

% Obtain grid of valid equilibria
ac_grid=linspace(0,1,l_grid);
% Player 1
y_grid=[];
for i1=1:length(ac_grid)
    for i2=1:length(ac_grid)
        u1=g1(ac_grid(i1),ac_grid(i2)); %Point to test
        ac_dev=ac_grid;
        ac_dev(i1)=[]; %All actions but current
        dev=max(g1(ac_dev,ac_grid(i2)));
        idx=(u1-(1-delta)*dev-delta*threat_v1<0);
        if sum(idx)==0 && u1>=threat_v1 %No profitable deviation
            y_grid(end+1,:)=[ac_grid(i1),ac_grid(i2),u1];
        end
    end
end

% Player 2
z_grid=[];
for i1=1:length(ac_grid)
    for i2=1:length(ac_grid)
        u2=g2(ac_grid(i1),ac_grid(i2)); %Point to test
        ac_dev=ac_grid;
        ac_dev(i2)=[]; %All actions but current
        dev=max(g2(ac_grid(i1),ac_dev));
        idx=(u2-(1-delta)*dev-delta*threat_v2<0);
        if sum(idx)==0 && u2>=threat_v2 %No profitable deviation
            z_grid(end+1,:)=[ac_grid(i1),ac_grid(i2),u2];
        end
    end
end

% Randomly take N_com values from y_grid, z_grid only
if size(y_grid,1)>N_com
    y_grid=y_grid(randperm(size(y_grid,1)),:);
    y_grid=y_grid(1:N_com,:);
end
if size(z_grid,1)>N_com
    z_grid=z_grid(randperm(size(z_grid,1)),:);
    z_grid=z_grid(1:N_com,:);
end
% Refine valid points: only those in which both players gain
valid_grid=intersect(y_grid(:,1:2),z_grid(:,1:2),'rows');

% Obtain region points
v1c=zeros(l_grid^2,1);
v2c=zeros(l_grid^2,1);
for i=1:l_grid
    v1c((i-1)*l_grid+1:i*l_grid)=g1(ac_grid(i),ac_grid)';
    v2c((i-1)*l_grid+1:i*l_grid)=g2(ac_grid(i),ac_grid)';
end
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