function[valid_grid,v_out]=sample_random_corr(N_com,Nmax_iter,delta,g1,g2,u_1,u_2,threat_str_1, threat_str_2,threat_v1,threat_v2);

stalled=0; % Set to 1 if algorithm stalls
valid_grid=[];
v1c=[];
v2c=[];
i=1;

while i<=N_com && stalled==0 % Test for a determined number of communications
    % Player 1
    searching_point=1;
    phi_test=[];
    j=0;
    while searching_point==1 && stalled==0
        j=j+1;
        % Generate random phi value
        r_v=rand(1,3);
        aux=[0 sort(r_v) 1];
        phi_test=diff(aux);
        % Check if satisfies p1 restrictions
        v1o=g1(phi_test);
        restr_1=(1-delta)*(phi_test(1)*(u_1(2,1)-u_1(1,1))+phi_test(2)*(u_1(2,2)-u_1(1,2)))+delta*(threat_v1-v1o)*(phi_test(1)+phi_test(2));
        restr_2=(1-delta)*(phi_test(3)*(u_1(1,1)-u_1(2,1))+phi_test(4)*(u_1(1,2)-u_1(2,2)))+delta*(threat_v1-v1o)*(phi_test(3)+phi_test(4));
        if restr_1<=0 && restr_2<=0 %No profitable deviation
            if v1o>=threat_v1 % Better eq than threat
                searching_point=0; %Exit loop: valid point found
            end
        end
        % Add to region tested
        v1c(end+1)=g1(phi_test);
        v2c(end+1)=g2(phi_test);
        if j>=Nmax_iter
            stalled=1;
        end
    end
    if stalled==0
        % Test if it is valid for player 2 (send it the point and check)
        v2o=g2(phi_test);
        restr_1=(1-delta)*(phi_test(1)*(u_2(1,2)-u_2(1,1))+phi_test(3)*(u_2(2,2)-u_2(2,1)))+delta*(threat_v2-v2o)*(phi_test(1)+phi_test(3));
        restr_2=(1-delta)*(phi_test(2)*(u_2(1,1)-u_2(1,2))+phi_test(4)*(u_2(2,1)-u_2(2,2)))+delta*(threat_v2-v2o)*(phi_test(2)+phi_test(4));
        if restr_1<=0 && restr_2<=0 %Valid for player 2
            if v2o>=threat_v2 % Better eq than threat
                valid_grid(end+1,:)=phi_test; %Add to valid grid
            end
        end
    end
    
    % Player 2
    searching_point=1;
    phi_test=[];
    j=0;
    while searching_point==1 && stalled==0
        j=j+1;
        % Generate random phi value
        r_v=rand(1,3);
        aux=[0 sort(r_v) 1];
        phi_test=diff(aux);
        % Check if satisfies p2 restrictions
        v2o=g2(phi_test);
        restr_1=(1-delta)*(phi_test(1)*(u_2(1,2)-u_2(1,1))+phi_test(3)*(u_2(2,2)-u_2(2,1)))+delta*(threat_v2-v2o)*(phi_test(1)+phi_test(3));
        restr_2=(1-delta)*(phi_test(2)*(u_2(1,1)-u_2(1,2))+phi_test(4)*(u_2(2,1)-u_2(2,2)))+delta*(threat_v2-v2o)*(phi_test(2)+phi_test(4));
        if restr_1<=0 && restr_2<=0 %No profitable deviation
            if v2o>=threat_v2 % Better eq than threat
                searching_point=0; %Exit loop: valid point found
            end
        end
        % Add to region tested
        v1c(end+1)=g1(phi_test);
        v2c(end+1)=g2(phi_test);
        if j>=Nmax_iter
            stalled=1;
        end
    end
    if stalled==0
        % Test if it is valid for player 1 (send it the point and check)
        v1o=g1(phi_test);
        restr_1=(1-delta)*(phi_test(1)*(u_1(2,1)-u_1(1,1))+phi_test(2)*(u_1(2,2)-u_1(1,2)))+delta*(threat_v1-v1o)*(phi_test(1)+phi_test(2));
        restr_2=(1-delta)*(phi_test(3)*(u_1(1,1)-u_1(2,1))+phi_test(4)*(u_1(1,2)-u_1(2,2)))+delta*(threat_v1-v1o)*(phi_test(3)+phi_test(4));
        if restr_1<=0 && restr_2<=0 %Valid for player 1
            if v1o>=threat_v1 % Better eq than threat
                valid_grid(end+1,:)=phi_test; %Add to valid grid
            end
        end
    end
    i=i+1;
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