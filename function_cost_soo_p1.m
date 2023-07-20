function [value_of_function]=function_cost_soo_p1(s)

param=getGlobalx;
lambda=param.lambda;
valid_grid=param.valid_grid;
g1=param.g1;
g2=param.g2;
ac_dev=param.ac_dev;
v1c=param.v1c;
v2c=param.v2c;
delta=param.delta;
threat_v1=param.threat_v1;
threat_v2=param.threat_v2;

gamma_i=[];
gamma_rest=[];
valid_i=0;

% Check action for player 1
ac_test=s;
u1=g1(ac_test(1),ac_test(2));
% Check if satisfies p1 restrictions
ac_dev(ac_dev==ac_test(1))=[]; %All actions but current
dev=max(g1(ac_dev,ac_test(2)));
idx=(u1-(1-delta)*dev-delta*threat_v1<0);
if sum(idx)==0 && u1>=threat_v1 %No profitable deviation
    gamma_i=norm(u1-threat_v1);
    valid_i=1;
else
    gamma_i=-norm(u1-threat_v1);
end
% Add to region tested
v1c(end+1)=g1(ac_test(1),ac_test(2));
v2c(end+1)=g2(ac_test(1),ac_test(2));

% Test if it is valid for player 2 (send it the point and check)
u2=g2(ac_test(1),ac_test(2));
ac_dev(ac_dev==ac_test(2))=[]; %All actions but current
dev=max(g2(ac_test(1),ac_dev));
idx=(u2-(1-delta)*dev-delta*threat_v2<0);

if sum(idx)==0 && u2>=threat_v2 %Valid for player 2
    if valid_i==1
        valid_grid(end+1,:)=ac_test; %Add to valid grid: valid for all
    end
    gamma_rest=norm(u2-threat_v2);
else
    gamma_rest=-norm(u2-threat_v2);
end

%Update global
setGlobalgrid(valid_grid,v1c,v2c)

%output value

value_of_function=lambda*gamma_i+(1-lambda)*gamma_rest;

return;