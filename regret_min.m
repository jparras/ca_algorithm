function[u1v,u2v,a1v,a2v]=regret_min(na1,na2,u_1,u_2,n_iter);
W1=zeros(na1,1);
W2=zeros(na2,1);
a1v=zeros(n_iter,na1);
a2v=zeros(n_iter,na2);
u1v=zeros(n_iter,1);
u2v=zeros(n_iter,1);
for t=1:n_iter
    % Obtain actions realizations for the current iteration
    if max(W1)<=0 && max(W2)<=0 %No regret: use random action
        a1 = round(1 + (na1-1).*rand(1)); %Random distrubution in [1, na1]
        a2 = round(1 + (na2-1).*rand(1)); %Random distribution in [1, na2]
    elseif max(W1)>0 && max(W2)>0 %Both have positive regrets
        aux1=W1;
        aux1(aux1<=0)=0;
        p1_aux=aux1/sum(aux1);
        a1_aux=mnrnd(1,p1_aux);
        a1=find(a1_aux==1);
        
        aux2=W2;
        aux2(aux2<=0)=0;
        p2_aux=aux2/sum(aux2);
        a2_aux=mnrnd(1,p2_aux);
        a2=find(a2_aux==1);
        
    elseif max(W1)<=0 && max(W2)>0
        
        a1 = round(1 + (na1-1).*rand(1));
        aux2=W2;
        aux2(aux2<=0)=0;
        p2_aux=aux2/sum(aux2);
        a2_aux=mnrnd(1,p2_aux);
        a2=find(a2_aux==1);
        
    elseif max(W1)>0 && max(W2)<=0
        
        aux1=W1;
        aux1(aux1<=0)=0;
        p1_aux=aux1/sum(aux1);
        a1_aux=mnrnd(1,p1_aux);
        a1=find(a1_aux==1);
        a2 = round(1 + (na2-1).*rand(1));
    end
    % Update actions vector
    a1v(t,a1)=1;
    a2v(t,a2)=1;
    % Update payoffs
    u1v(t)=u_1(a1,a2);
    u2v(t)=u_2(a1,a2);
    % Update regrets
    for i=1:na1
        W1(i)=W1(i)+u_1(i,a2)-u_1(a1,a2);
    end
    for i=1:na2
        W2(i)=W2(i)+u_2(a1,i)-u_2(a1,a2);
    end
end