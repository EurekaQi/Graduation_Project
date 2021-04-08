clear
clc
tic
V=5;  %������
N=13; %�ܽڵ�����
C=6; %��������
demands=[0,1.2,1.7,1.5,1.4,1.7,1.4,1.2,1.9,1.8,1.6,1.7,1.1]; %������
x=[81.5,87,75,85,89,77,76,87,73,77,73,91,92];
y=[41.5,37,53,52,41,58,45,53,38,38,31,47,44];
axis=[x' y']; %��������
for i=1:N
    for j=1:N
        Dij(i,j)= sqrt((axis(i,1)-axis(j,1))^2+(axis(i,2)-axis(j,2))^2);
    end
end
%% ���߱���
xijk=binvar(N,N,V,'full');%i��j�ڵ�֮���Ƿ��ɵ�k������������
yik=binvar(N,V,'full'); %��k�����Ƿ񾭹�i�ڵ�
zik=sdpvar(N,V,'full'); %����k�ڷ���i�ڵ�󣬳��ӵ�ʣ������
%% Ŀ�꺯��
obj=0;
for i=1:N
    for j=1:N
        for k=1:V
            obj=obj+Dij(i,j)*xijk(i,j,k);
        end
    end
end
%% Լ������
Con=[ ];
for i = 2:N 
    Con = [Con; sum(yik(i,:)) == 1];
end

Con = [Con; sum(yik(1,:)) <= V];

for k = 1:V
    Con = [Con; sum(demands(:).*yik(:,k)) <= C];
end

for i = 1:N
    for j = 1:N
        for k = 1:V
            if i == j
                Con = [Con; xijk(i,j,k) == 0];
            end
        end
    end
end

for i = 1:N
    for k = 1:V
        Con = [Con; sum(xijk(i,:,k)) == sum(xijk(:,i,k))];
    end
end

for j = 1:N
    for k = 1:V
        Con = [Con; sum(xijk(:,j,k)) == yik(j,k)];
    end
end

for i = 1:N
    for k = 1:V
        Con = [Con; sum(xijk(i,:,k)) == yik(i,k)];
    end
end

for i = 2:N
    for j = 2:N
        for k = 1:V
            if i ~= j
                if demands(i)+demands(j) <= C
                    Con = [Con; zik(i,k)-zik(j,k)+C*xijk(i,j,k) <= C-demands(i)];
                end
            end
        end
    end
end

for i = 2:N
    for k = 1:V
        Con = [Con; zik(i,k) <= C];
        Con = [Con; zik(i,k) >= demands(i)];
    end
end
%% ���
ops = sdpsettings( 'solver','gurobi');
sol = solvesdp(Con,obj,ops);
obj = double(obj);
xijk = double(xijk);
yik = double(yik);
zik = value(zik);