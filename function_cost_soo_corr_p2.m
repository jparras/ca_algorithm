function [value_of_function]=function_cost_soo_corr_p2(s)

param=getGlobalx;
lambda=param.lambda;
valid_grid=param.valid_grid;
g1=param.g1;
g2=param.g2;
v1c=param.v1c;
v2c=param.v2c;
delta=param.delta;
threat_v1=param.threat_v1;
threat_v2=param.threat_v2;
u_1=param.u_1;
u_2=param.u_2;

gamma_i=[];
gamma_rest=[];
value_of_function=[];
valid_i=0;

% Check that s is valid!
if sum(s)<=1
    phi_test=[s(1) s(2) s(3) 1-sum(s)];
    % Check action for player 2
    v2o=g2(phi_test);
    restr_1=(1-delta)*(phi_test(1)*(u_2(1,2)-u_2(1,1))+phi_test(3)*(u_2(2,2)-u_2(2,1)))+delta*(threat_v2-v2o)*(phi_test(1)+phi_test(3));
    restr_2=(1-delta)*(phi_test(2)*(u_2(1,1)-u_2(1,2))+phi_test(4)*(u_2(2,1)-u_2(2,2)))+delta*(threat_v2-v2o)*(phi_test(2)+phi_test(4));
    if restr_1<=0 && restr_2<=0 %No profitable deviation
        if v2o>=threat_v2 % Better eq than threat
            gamma_i=norm(v2o-threat_v2);
            valid_i=1;
        end
    end
    if isempty(gamma_i)
        gamma_i=-norm(v2o-threat_v2);
    end
    % Add to region tested
    v1c(end+1)=g1(phi_test);
    v2c(end+1)=g2(phi_test);

    % Test if it is valid for player 1 (send it the point and check)
    v1o=g1(phi_test);
    restr_1=(1-delta)*(phi_test(1)*(u_1(2,1)-u_1(1,1))+phi_test(2)*(u_1(2,2)-u_1(1,2)))+delta*(threat_v1-v1o)*(phi_test(1)+phi_test(2));
    restr_2=(1-delta)*(phi_test(3)*(u_1(1,1)-u_1(2,1))+phi_test(4)*(u_1(1,2)-u_1(2,2)))+delta*(threat_v1-v1o)*(phi_test(3)+phi_test(4));
    if restr_1<=0 && restr_2<=0 %Valid for player 1
        if v1o>=threat_v1 % Better eq than threat
            if valid_i==1
                valid_grid(end+1,:)=phi_test; %Add to valid grid: valid for all
            end
            gamma_rest=norm(v1o-threat_v1);
        end
    end
    if isempty(gamma_rest)
        gamma_rest=-norm(v1o-threat_v1);
    end

    %Update global
    setGlobalgrid(valid_grid,v1c,v2c)

    %output value

    value_of_function=lambda*gamma_i+(1-lambda)*gamma_rest;
else
    value_of_function=-10; %High
end

return;