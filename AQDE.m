clc;
clear;
%% 论文相关信息
%Ashish Ranjan Hota, Pat Ankit. An adaptive quantum-inspired differential evolution...
%algorithm for 0–1 knapsack problem. 2010 second world congress on nature and...
%biologically inspired computing (NaBIC). IEEE, 2010.
% 论文部分地方根据下面的文献进行了适当的修改
% Hu Peng Zhijian Wu, Peng Shao, Changshou Deng. Dichotomous Binary Differential...
%Evolution for Knapsack Problems [J]. Mathematical Problems in Engineering, 2016, 5732489.

Ztt = 0;
ZZ = [];
ZZ_con_1 = [];
tic
for kkk = 1:30
%% 背包问题参数
load('kp_sc_1000');
%the Instances should include
% profits of all items $P = [p_1,p_2,\dots,p_n]$
% weights of all items $W = [w_1,w_2,\dots,w_n]$
% knapsack capacity $C$
% number of items $n$

E1 = [P./W;1:n];
E2 = sortrows(E1',1,'descend')';
E = E2(2,:);
P = P(E);
W = W(E);
%% AQDE
%算法参数
%最大迭代次数
MaxIter = 200;
%种群个体数量
pop = 30;
% 论文中的theta , 范围为0到2pi
X = rand(pop,n)*2*pi;
% 对应的
Y = (rand(pop,n)<sin(X).*sin(X))*1;


Record = [];
%修复不可行解
for i = 1:pop
    j = n;
    while Y(i,:)*W'>=C
        Y(i,j) = 0;
        j = j-1;
    end
%     while Y(i,:)*W'>=C 
%         T11 = find(Y(i,:)>0);
%         t11 = length(T11);
%         t12 = randperm(t11,1);
%         j = T11(t12);
%         Y(i,j) = 0;
%     end
    while Y(i,:)*W'< C
        T12 = find(Y(i,:)==0);
        t13 = length(T12);
        t14 = randperm(t13,1);
        j = T12(t14);
        T13 = Y(i,:);
        T13(j) = 1;
        if T13*W'< C
            Y(i,j) = 1;
        else
            break
        end
    end
end

t = 1;
while t <= MaxIter
    for i = 1:pop
        %变异
        T21 = randperm(pop,3);
        t21 = T21(1);
        t22 = T21(2);
        t23 = T21(3);
        F = rand*rand*0.1;
        Theta(i,:) = X(t21,:)+F*(X(t22,:)-X(t23,:));
        % 交叉
        %cr = normrnd(0.5,0.0375);%原始文献参数
        cr = 0.5+0.0375*randn;% 参数来源 Dichotomous Binary Differential Evolution for Knapsack Problems
        Irand = randi(n);
        for j = 1:n
            if rand < cr || j == Irand
                Theta_c(i,j) = Theta(i,j);
            else
                Theta_c(i,j) = X(i,j);
            end
        end
    end
    %选择
    % 对应的
    Y_c = (rand(pop,n)<sin(Theta_c).*sin(Theta_c))*1;

    %修复不可行解
    for i = 1:pop
        j = n;
        while Y_c(i,:)*W'>=C
            Y_c(i,j) = 0;
            j = j-1;
        end
%         while Y_c(i,:)*W'>=C
%             T11 = find(Y_c(i,:)>0);
%             t11 = length(T11);
%             t12 = randperm(t11,1);
%             j = T11(t12);
%             Y_c(i,j) = 0;
%         end
        while Y_c(i,:)*W'< C
            T12 = find(Y_c(i,:)==0);
            t13 = length(T12);
            t14 = randperm(t13,1);
            j = T12(t14);
            T13 = Y_c(i,:);
            T13(j) = 1;
            if T13*W'< C
                Y_c(i,j) = 1;
            else
                break
            end
        end
    end
    for i = 1:pop
        if Y_c(i,:)*P'> Y(i,:)*P'
            Y(i,:) = Y_c(i,:);
            X(i,:) = Theta_c(i,:);
        end
    end
    [a,b] = max(Y*P');
    Record = [Record,a];
    t = t+1;
end

%plot(1:MaxIter,Record)

ZZ = [ZZ,max(Record)];
ZZ_con_1 = [ZZ_con_1;Record];
end
Ztt = toc; %计算时长
ZZ_con = sum(ZZ_con_1/kkk); %收敛曲线
max(ZZ) %多次独立计算的best solution