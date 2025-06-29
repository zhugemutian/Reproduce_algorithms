clc;
clear;
% Kennedy, J.; Eberhart, R.C. Discrete binary version of the particle swarm...
% algorithm. In Proceedings of the IEEE International Conference on Systems, ...
% Man and Cybernetics, 1997, 5: 4104�C4108.

% ���Ĳ��ֵط�������������׽������ʵ����޸�
% Hu Peng Zhijian Wu, Peng Shao, Changshou Deng. Dichotomous Binary Differential...
%Evolution for Knapsack Problems [J]. Mathematical Problems in Engineering, 2016, 5732489.

Ztt = 0;
ZZ = [];
ZZ_con_1 = [];
tic
for kkk = 1:30
%% �����������
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
%% BPSO 
%�㷨����
%����������
MaxIter = 200;
%ѧϰ����
c1 = 2;
c2 = 2;
%��Ⱥ��������
pop = 30;
%����ϵ��
%omega = 0.5;%Ҳ�������ŵ��������𽥼�С
%�㷨��ʼ��
%�ٶ��Ͻ�
Vmax = 6; %�����ٶȹ�����С��������ֲ����Ž��޷�����
%�ٶȳ�ʼ��
%�ٶȷֲ�Ϊ[-Vmax, Vmax]
V = rand(pop,n)*Vmax-1.6*Vmax; %���һ��, ���ٱ�ѡ�����Ʒ, ������ִ������Ϊ�����н�����...
%1.4��ֵ���趨��C��sum(W)�ı�ֵ���
%ѡ�����
SV = 1./(1+exp(-V));
%�ٶ�����Ӧ��λ��
X = [rand(pop,n)<SV]*1;

for i1 = 1:pop
    if X(i1,:)*W' <= C %�޳������н�
        pBest_value(i1) = X(i1,:)*P';
        pBest(i1,:) = X(i1,:);
    else %���ڲ����н���ȫ0���滻
        pBest_value(i1) = 0;
        pBest(i1,:) = -1*ones(1,n);
    end
end
[a,b] = max(pBest_value);
gBest_value = a;
gBest = X(b,:);
%��¼�㷨����
Record = gBest_value;
%�㷨ѭ��
for Iter = 2:MaxIter
    %�ٶȸ���
    for i1 = 1:pop %��һ���Ե�����Ϊ��ÿ��Ԫ�ص��������ͬ
        for i2 = 1:n
            V(i1,i2) = V(i1,i2)+c1*rand*(X(i1,i2)-pBest(i1,i2))+c2*rand*(X(i1,i2)-gBest(i2));
            %�����ٶȽ��б߽��޸�
            if V(i1,i2)>Vmax
                V(i1,i2) = Vmax;
            end
            if V(i1,i2)<-Vmax
                V(i1,i2) = -Vmax;
            end
        end
    end
	%����ѡ�����
    SV = 1./(1+exp(-V));
    %����λ��
    X = [rand(pop,n)<SV]*1;
    for i = 1:pop
        j = n;
        while X(i,:)*W'>C
            if X(i,j) == 1;
                X(i,j) = 0;
            end
            j = j-1;
        end
    end
    %���¸������ź�ȫ������
    for i1 = 1:pop
        if X(i1,:)*W' <= C %�޳������н�
            if X(i1,:)*P' > pBest_value(i1)
                pBest_value(i1) = X(i1,:)*P';
                pBest(i1,:) = X(i1,:);
            else %���ڲ����н���ȫ0���滻
                pBest_value(i1) = 0;
                pBest(i1,:) = -1*ones(1,n);
            end
        end
    end
    [a,b] = max(pBest_value);
    gBest_value = max(Record(Iter-1),a);
    gBest = X(b,:);
    Record = [Record,gBest_value];
end
%plot(1:MaxIter,Record)

ZZ = [ZZ,max(Record)];
ZZ_con_1 = [ZZ_con_1;Record];
end
Ztt = toc; %����ʱ��
ZZ_con = sum(ZZ_con_1/kkk); %��������
max(ZZ) %��ζ��������best solution