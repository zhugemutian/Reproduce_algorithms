clc;
clear;
%% ���������Ϣ
% Ismail M. Ali, Daryl Essam, Kathryn Kasmarik. Novel binary differential...
%evolution algorithm for knapsack problems [J]. Information Sciences, 2021, 542: 177-194.

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

%% Nbin-DE
%�㷨����
%����������
MaxIter = 200;
%��Ⱥ��������
pop = 30;


%������Դ������
% scale factor (F), crossover rate (cr)
F = 1;
cr = 0.5;
% Similarity threshold, Maximum value of threshold, Minimum value of threshold
delta = 0.2;
delta_max = 0.9*n;
delta_min = 0.1*n;
%flip����δ����, ���ж���
flip = 0.5;

%BMV parameter
Ind = 0.9;
Gen = 0.8;

%Section 3.1
con_pop = rand(pop,n);

bin_pop = zeros(pop,n);
for i = 1:pop
    t11 = sum(con_pop(i,:))/n;
    bin_pop(i,:) = (con_pop(i,:)>t11).*1;
end

%Section 3.2

%�洢��j-th��Ʒ������ȵ���Ʒ��� SameWeight
SameWeight = zeros(n,round(n/20));
for j = 1:n
    T11 = find(W==W(j));
    t12 = length(T11);
    SameWeight(j,1:t12) = T11;
end
T12 = sum(SameWeight);
T13 = find(T12<1);
if ~isempty(T13)
    SameWeight = SameWeight(:,1:T13(1)-1); 
end
Num_SameWeight = sum(sign(SameWeight)'); %ͳ��ÿ��sameweight������
vio = (bin_pop*W')-C;
for i = 1:pop
    if vio(i) > 0 %(Alg. 2, lines 6-13)
        j = n;
        while vio(i)>0
            %������Ʒ˳��ͼ�ֵ�ܶ�˳����ͬ, ���԰�����Ʒ��������޳���Ʒ
            %���intancesû�а�����Ʒ��ֵ�ܶȽ�������, ������Ҫ��������
            if bin_pop(i,j)==1
                bin_pop(i,j) = 0;
                con_pop(i,j) = 0;
                vio(i) = bin_pop(i,:)*W'-C;
            end
            j = j-1;
        end
    end
    if vio(i) < 0 %(Alg. 2, lines 14-20)
        for j = 1:n
            if bin_pop(i,j) == 0 && W(j)+vio(i)<0
                bin_pop(i,j) = 1;
                con_pop(i,j) = 0;
                vio(i) = bin_pop(i,:)*W'-C;
            end
        end
    end
    for j = 1:n
        if bin_pop(i,j)==1
            if ~isempty(Num_SameWeight)
                for k = 1:Num_SameWeight(j)
                    t21 = SameWeight(j,k);
                    if t21 < j && bin_pop(i,t21) == 0
                        %���ཻ��
                        bin_pop(i,j) = 0; 
                        bin_pop(i,t21) = 1;
                        %������ֵ
                        t22 = con_pop(i,j);
                        con_pop(i,j) = con_pop(i,t21);
                        con_pop(i,t21) = t22;
                        break
                    end
                end
            end
        end
    end
end
% BMV operator 
%�ο�����
%Ismail M. Ali, Daryl Essam, Kathryn Kasmarik. A novel differential...
%evolution mapping technique for generic combinatorial optimization
%problems...Applied Soft Computing, 2019, 80: 297-309.
[a,b] = max(bin_pop*P');
Record = a;
Best_ind_1 = [con_pop(b,:);1:n];
Best_ind_2 = [sortrows(Best_ind_1',1)';1:n];
Best_ind_3 = sortrows(Best_ind_2',2)';
Best_ind = Best_ind_3(3,:);

t = 2;
while t <= MaxIter 
    for i = 1:pop
        %���������������ȵ�����
        T21 = randperm(pop,3);
        Vig = con_pop(T21(1),:)+F*(con_pop(T21(2),:)-con_pop(T21(3),:)); %(eq. 2)
        %���н������
        aj = randi(n);
        for j = 1:n
            if rand <= cr || j == aj
                Uig(j) = Vig(j);
            else
                Uig(j) = con_pop(i,j);
            end
        end
        %�ٴν�con_popӳ��Ϊbin_pop
        T22 = [Uig;1:n];
        T23 = sortrows(T22',1)';
        T24 = [T23;1:n];
        T25 = sortrows(T24',2)';
        Re_mapping = T25(3,:);
        if rand<= Ind
            for j = 1:n
                if rand<=Gen
                    t31 = Best_ind(j); %t21=5 %j=3
                    t32 = find(Re_mapping==t31); %t22 = 6-th λ��
                    t33 = Re_mapping(t32); %t23 = 6 ��ֵ
                    Re_mapping(j) = Best_ind(j);
                    Re_mapping(t32) = t33;
                end
            end
        end
        avg_Re_mapping = sum(Re_mapping)/n;
        New_Re_mapping = (Re_mapping>avg_Re_mapping)*1;
        vio(i) = New_Re_mapping*W'-C;
        %�ٴ�̰���޸�
        if vio(i) > 0 %(Alg. 2, lines 6-13)
            j = n;
            while vio(i)>0
                if New_Re_mapping(j) == 1
                    New_Re_mapping(j) = 0;
                    vio(i) = New_Re_mapping*W'-C;
                end
                j = j-1;
            end
        end
        if vio(i) < 0 %(Alg. 2, lines 14-20)
            for j = 1:n
                if New_Re_mapping(j) == 0 && W(j)+vio(i)<0
                    New_Re_mapping(j) = 1;
                    vio(i) = New_Re_mapping*W'-C;
                end
            end
        end
        for j = 1:n
            if New_Re_mapping(j)==1
                for k = 1:Num_SameWeight(j)
                    t21 = SameWeight(j,k);
                    if t21 < j && New_Re_mapping(t21) == 0
                        %���ཻ��
                        New_Re_mapping(j) = 0; 
                        New_Re_mapping(t21) = 1;
                        break
                    end
                end
            end
        end
        %selected operator
        if New_Re_mapping*P'>=bin_pop(i,:)*P'
            con_pop_2(i,:) = Uig;
            bin_pop_2(i,:) = New_Re_mapping;
        else
            bin_pop_2(i,:) = bin_pop(i,:);
            con_pop_2(i,:) = con_pop(i,:);
            vio(i) = bin_pop_2(i,:)*W'-C;
            if vio(i) > 0 %(Alg. 2, lines 6-13)
                j = n;
                while vio(i)>0
                    %������Ʒ˳��ͼ�ֵ�ܶ�˳����ͬ, ���԰�����Ʒ��������޳���Ʒ
                    %���intancesû�а�����Ʒ��ֵ�ܶȽ�������, ������Ҫ��������
                    if bin_pop_2(i,j)==1
                        bin_pop_2(i,j) = 0;
                        con_pop_2(i,j) = 0;
                        vio(i) = bin_pop_2(i,:)*W'-C;
                    end
                    j = j-1;
                end
            end
            if vio(i) < 0 %(Alg. 2, lines 14-20)
                for j = 1:n
                    if bin_pop_2(i,j) == 0 && W(j)+vio(i)<0
                        bin_pop_2(i,j) = 1;
                        con_pop_2(i,j) = 0;
                        vio(i) = bin_pop_2(i,:)*W'-C;
                    end
                end
            end
            for j = 1:n
                if bin_pop_2(i,j)==1
                    for k = 1:Num_SameWeight(j)
                        t21 = SameWeight(j,k);
                        if t21 < j && bin_pop_2(i,t21) == 0
                            %���ཻ��
                            bin_pop_2(i,j) = 0; 
                            bin_pop_2(i,t21) = 1;
                            %������ֵ
                            t22 = con_pop_2(i,j);
                            con_pop_2(i,j) = con_pop_2(i,t21);
                            con_pop_2(i,t21) = t22;
                            break
                        end
                    end
                end
            end
        end
    end
    bin_pop = bin_pop_2;
    con_pop = con_pop_2;
    [a,b] = max(bin_pop*P');
    Record = [Record,a];
    %��֤��Ⱥ������
    %�Ը�����밴��Ŀ�꺯�����н�������
    fitness_bin_pop_1 = bin_pop*P';
    fitness_bin_pop_2 = [bin_pop,fitness_bin_pop_1];
    fitness_bin_pop_3 = sortrows(fitness_bin_pop_2, n+1,'descend');
    fitness_con_pop_2 = [con_pop,fitness_bin_pop_1];
    fitness_con_pop_3 = sortrows(fitness_con_pop_2, n+1,'descend');
    bin_pop = fitness_bin_pop_3(:,1:n);
    con_pop = fitness_con_pop_3(:,1:n);
    delta = min(delta_min+(delta_min+delta_max)*t/MaxIter,delta_max);
    for i = 2:pop
        diff(i,:) = abs(bin_pop(i,:)-bin_pop(i-1,:));
        if sum(diff(i,:)) < delta*n
            for j = 1:n
                if rand<flip && diff(i,j)==0
                    bin_pop(i,j) = 1-bin_pop(i,j);
                    if bin_pop(i,j) == 1
                        con_pop(i,j) = 0.6+(1-0.6)*rand;
                    else
                        con_pop(i,j) = 0.4*rand;
                    end
                end
            end
        end
    end
    t = t+1;
end
%plot(1:MaxIter,Record)

ZZ = [ZZ,max(Record)];
ZZ_con_1 = [ZZ_con_1;Record];

end
Ztt = toc; %����ʱ��
ZZ_con = sum(ZZ_con_1/kkk); %��������
max(ZZ) %��ζ��������best solution