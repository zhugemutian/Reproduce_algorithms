clc;
clear;
%% 论文相关信息
% Ervural, Bilal, and Huseyin Hakli. A binary reptile search algorithm based on ...
% transfer functions with a new stochastic repair method for 0–1 knapsack problems. 
% Computers & Industrial Engineering, 2023, 178: 109080.

Ztt = 0;
ZZ = [];
ZZ_con_1 = [];
tic
for kkk = 1:30
%% parameter setting
load('kp_sc_1000');
E1 = [P./W;1:n];
E2 = sortrows(E1',1,'descend')';
E = E2(2,:);
P = P(E);
W = W(E);

MaxIter = 200;
pop = 30;
UB = 100;
LB = -100;

alpha = 0.1;
beta = 0.005;

%SRate 
SRate = 200/n^1.45;
%SRate = 0.5;
%% 初始化
X = rand(pop,n)*(UB-LB)+LB;

%对应的binary 解 Tf1 transfer function
Y = (rand(pop,n)<=sqrt(abs(X/UB)))*1;
%对解进行修复
pw = C-Y*W';
for i = 1:pop
    if pw(i)<0
        while pw(i)<0
            T11 = find(Y(i,:)>0);
            if rand>= SRate
                at = max(T11);
            else
                t11 = randi(length(T11));
                at = T11(t11);
            end
            Y(i, at) = 0;
            pw(i) = C-Y(i,:)*W';
        end
    end
    if pw(i)>0
        terminate = 0;
        while terminate == 0
            dif = pw(i);
            candidates_1 = find(Y(i,:)==0);
            candidates_2 = find(W<=dif);
            candidate = intersect(candidates_1,candidates_2);
            if ~isempty(candidate)
                if rand>= SRate
                    bt = min(candidate);
                else
                    t12 = randi(length(candidate));
                    bt = candidate(t12);
                end
                Y(i, bt) = 1;
                pw(i) = C-Y(i,:)*W';
            else 
                terminate = 1;
            end
        end
    end
end

Record = [];
t = 1;
while t <= MaxIter
    %更新最优解
    [a,b] = max(Y*P');
    Best = a;
    Record = [Record,a];
    Best_ind = X(b,:);
    Best_Ind = Y(b,:);
    %更新参数
    r3 = randi([-1,1]);
    ES = 2*r3*(1-t/MaxIter);
    for i = 1:pop
        Mxi = sum(X(i,:))/n;
        for j = 1:n
            Rij = (Best_ind(j)-X(randi(pop),j))/(Best_ind(j)+eps);
            Pij =alpha+(X(i,j)-Mxi)/(Best_ind(j)*(UB-LB)+eps);
            etaij = Best_ind(j)*Pij;
            if t<MaxIter/4
                Xt1(i,j) = Best_ind(j)*(-etaij)*beta-Rij*rand;
            end
            if t<2*MaxIter/4 && t>MaxIter/4
                Xt1(i,j) = Best_ind(j)*X(randi(pop),j)*ES*rand;
            end
            if t<3*MaxIter/4 && t>2*MaxIter/4
                Xt1(i,j) = Best_ind(j)*Pij*rand;
            end
            if t>3*MaxIter/4
                Xt1(i,j) = Best_ind(j)-etaij*eps-Rij*rand;
            end
        end
        %映射得到binary solution
        Yt1(i,:) = (rand(1,n)<=sqrt(abs(Xt1(i,:)/UB)))*1;
        pwi = C-Yt1(i,:)*W';
        if pwi<0
            while pwi<0
                T11 = find(Yt1(i,:)>0);
                if rand>= SRate
                    at = max(T11);
                else
                    t11 = randi(length(T11));
                    at = T11(t11);
                end
                Yt1(i, at) = 0;
                pwi = C-Yt1(i,:)*W';
            end
        end
        if pwi>0
            terminate = 0;
            while terminate == 0
                dif = pwi;
                candidates_1 = find(Yt1(i,:)==0);
                candidates_2 = find(W<=dif);
                candidate = intersect(candidates_1,candidates_2);
                if ~isempty(candidate)
                    if rand>= SRate
                        bt = min(candidate);
                    else
                        t12 = randi(length(candidate));
                        bt = candidate(t12);
                    end
                    Yt1(i, bt) = 1;
                    pwi = C-Yt1(i,:)*W';
                else 
                    terminate = 1;
                end
            end
        end
    end
    %更新
    for i = 1:pop
        if Yt1(i,:)*P'>Y(i,:)*P'
            Y(i,:) = Yt1(i,:);
            X(i,:) = Xt1(i,:);
        end
    end
    t = t+1;
end
    
%plot(1:MaxIter,Record)

ZZ = [ZZ,max(Record)];
ZZ_con_1 = [ZZ_con_1;Record];
end
Ztt = toc; %计算时长
ZZ_con = sum(ZZ_con_1/kkk); %收敛曲线
max(ZZ) %多次独立计算的best solution
    
    
    
    
    