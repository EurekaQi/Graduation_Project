%%
clc,clear
%%
AllGens=cell(3,1);
AllGens{1}=struct('capacity',20,'cost',16);
AllGens{2}=struct('capacity',10,'cost',19);
AllGens{3}=struct('capacity',25,'cost',15);

AllLoads=cell(3,1);
AllLoads{1}=struct('demand',5,'cost',18);
AllLoads{2}=struct('demand',20,'cost',20);
AllLoads{3}=struct('demand',15,'cost',21);

B=[0 100 125;
   100 0 150;
   125 150 0];

F=[0 5 10;
   5 0 10;
   10 10 0];
%%
gen=sdpvar(3,1,'full');
dem=sdpvar(3,1,'full');
theta=sdpvar(3,1,'full');
Offer_Price=sdpvar(1,1,'full');
LMP=sdpvar(3,1,'full');
phi_gmin=sdpvar(3,1,'full');
phi_gmax=sdpvar(3,1,'full');
phi_dmin=sdpvar(3,1,'full');
phi_dmax=sdpvar(3,1,'full');
v=sdpvar(3,3,'full');
r=sdpvar(1,1,'full');
c_gen_min=binvar(3,1,'full');
c_gen_max=binvar(3,1,'full');
c_dem_min=binvar(3,1,'full');
c_dem_max=binvar(3,1,'full');
c_plmax=binvar(3,3,'full');
M=500000;
%%
const=[];

const=[const,Offer_Price-LMP(1)+phi_gmax(1)-phi_gmin(1)==0];
const=[const,AllGens{2}.cost-LMP(2)+phi_gmax(2)-phi_gmin(2)==0];
const=[const,AllGens{3}.cost-LMP(3)+phi_gmax(3)-phi_gmin(3)==0];

const=[const,LMP(1)-AllLoads{1}.cost+phi_dmax(1)-phi_dmin(1)==0];
const=[const,LMP(2)-AllLoads{2}.cost+phi_dmax(2)-phi_dmin(2)==0];
const=[const,LMP(3)-AllLoads{3}.cost+phi_dmax(3)-phi_dmin(3)==0];

const=[const,B(1,2)*(LMP(1)-LMP(2))+B(1,3)*(LMP(1)-LMP(3))+B(1,2)*(v(1,2)-v(2,1))+B(1,3)*(v(1,3)-v(3,1))==0];
const=[const,B(2,1)*(LMP(2)-LMP(1))+B(2,3)*(LMP(2)-LMP(3))+B(2,1)*(v(2,1)-v(1,2))+B(2,3)*(v(2,3)-v(3,2))==0];
const=[const,B(3,1)*(LMP(3)-LMP(1))+B(3,2)*(LMP(3)-LMP(2))+B(3,1)*(v(3,1)-v(1,3))+B(3,2)*(v(3,2)-v(2,3))+r==0];

for i=1:3
    const=[const,0<=gen(i)<=M*c_gen_min(i)];
    const=[const,0<=phi_gmin(i)<=M*(1-c_gen_min(i))];
    
    const=[const,0<=AllGens{i}.capacity-gen(i)<=M*c_gen_max(i)];
    const=[const,0<=phi_gmax(i)<=M*(1-c_gen_max(i))];
end

for i=1:3
    const=[const,0<=dem(i)<=M*c_dem_min(i)];
    const=[const,0<=phi_dmin(i)<=M*(1-c_dem_min(i))];
    
    const=[const,0<=AllLoads{i}.demand-dem(i)<=M*c_dem_max(i)];
    const=[const,0<=phi_dmax(i)<=M*(1-c_dem_max(i))];
end

for i=1:3
    for j=1:3
        if i~=j
            const=[const,0<=F(i,j)-B(i,j)*(theta(i)-theta(j))<=M*c_plmax(i,j)];
            const=[const,0<=v(i,j)<=M*(1-c_plmax(i,j))];
        end
    end
end

const=[const,gen(1)-dem(1)-B(1,2)*(theta(1)-theta(2))-B(1,3)*(theta(1)-theta(3))==0];
const=[const,gen(2)-dem(2)-B(2,1)*(theta(2)-theta(1))-B(2,3)*(theta(2)-theta(3))==0];
const=[const,gen(3)-dem(3)-B(3,1)*(theta(3)-theta(1))-B(3,2)*(theta(3)-theta(2))==0];

const=[const,theta(3)==0];

%const=[const,Offer_Price==AllGens{1}.cost];


Obj=(AllGens{1}.cost-LMP(1))*gen(1);
Ops=sdpsettings('verbose',1,'solver','gurobi');
Result=solvesdp(const,Obj,Ops);















