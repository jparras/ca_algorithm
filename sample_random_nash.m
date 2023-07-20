function[valid_grid,v_out]=sample_random_nash(N_com,Nmax_iter,delta,g1,g2,threat_v1,threat_v2,ac_dev,threat_str_1, threat_str_2);

stalled=0; % Set to 1 if algorithm stalls
valid_grid=[];
v1c=[];
v2c=[];
i=1;
while i<=N_com && stalled==0 % Test for a determined number of communications
    % Player 1
    searching_point=1;
    ac_test=[];
    j=0;
    while searching_point==1 && stalled==0
        j=j+1;
        % Generate random action value
        ac_test=rand(1,2);
        u1=g1(ac_test(1),ac_test(2));
        % Check if satisfies p1 restrictions
        ac_dev(ac_dev==ac_test(1))=[]; %All actions but current
        dev=max(g1(ac_dev,ac_test(2)));
        idx=(u1-(1-delta)*dev-delta*threat_v1<0);
        if sum(idx)==0 && u1>=threat_v1 %No profitable deviation
            searching_point=0;
        end
        % Add to region tested
        v1c(end+1)=g1(ac_test(1),ac_test(2));
        v2c(end+1)=g2(ac_test(1),ac_test(2));
        if j>=Nmax_iter
            stalled=1;
        end
    end
    if stalled==0
        % Test if it is valid for player 2 (send it the point and check)
        u2=g2(ac_test(1),ac_test(2));
        ac_dev(ac_dev==ac_test(2))=[]; %All actions but current
        dev=max(g2(ac_test(1),ac_dev));
        idx=(u2-(1-delta)*dev-delta*threat_v2<0);        
        if sum(idx)==0 && u2>=threat_v2 %Valid for player 2
            valid_grid(end+1,:)=ac_test; %Add to valid grid
        end
    end
    
    % Player 2
    searching_point=1;
    j=0;
    while searching_point==1 && stalled==0
        j=j+1;
        % Generate random action value
        ac_test=rand(1,2);
        u2=g2(ac_test(1),ac_test(2));
        % Check if satisfies p2 restrictions
        ac_dev(ac_dev==ac_test(2))=[]; %All actions but current
        dev=max(g2(ac_test(1),ac_dev));
        idx=(u2-(1-delta)*dev-delta*threat_v2<0);
        if sum(idx)==0 && u2>=threat_v2 %No profitable deviation
            searching_point=0;
        end
        % Add to region tested
        v1c(end+1)=g1(ac_test(1),ac_test(2));
        v2c(end+1)=g2(ac_test(1),ac_test(2));
        if j>=Nmax_iter
            stalled=1;
        end
    end
    if stalled==0
        % Test if it is valid for player 1 (send it the point and check)
        u1=g1(ac_test(1),ac_test(2));
        ac_dev(ac_dev==ac_test(1))=[]; %All actions but current
        dev=max(g1(ac_dev,ac_test(2)));
        idx=(u1-(1-delta)*dev-delta*threat_v1<0);        
        if sum(idx)==0 && u1>=threat_v1 %Valid for player 1
            valid_grid(end+1,:)=ac_test; %Add to valid grid
        end
    end
    i=i+1;
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