function [node,elem,Db,Nb_L,Nb_R,Cb]=generation_elem_node(h)
% Omega初始边界及单元
% 3----6----9
% |\ 4 |\ 8 |
% | 3\ | 7\ |
% 2----5----8
% |\ 2 |\ 6 |
% | 1\ | 5\ |
% 1----4----7  
if 1/h ~= round(1/h)
    error('输入的划分不好,请重新修改\n');
end
N=1/h+1; % 每条边的节点数
nt=N-1;  % 每条边的单元数
elem=zeros(nt^2*2,3);
for i=1:1:nt
    start=(i-1)*N;numt=(i-1)*2*nt;
    % [2 1 4];      % 生成正向集
    elem(numt+1:2:numt+2*nt,:)=[(start+2:start+nt+1)',(start+1:start+nt)',(start+1+N:start+N+nt)'];
    % [4 5 2]
    elem(numt+2:2:numt+2*nt,:)=[(start+1+N:start+N+nt)',(start+2+N:start+N+nt+1)',(start+2:start+nt+1)'];
end

node=zeros(N^2,2);
a=[0:h:1]';
b=ones(N,1);
node(:,1)=kron(a,b);
node(:,2)=kron(b,a);

% 生成边界点
Db=[];Nb_L=[];Nb_R=[];Cb=[];
left=[1:N-1;2:N]';
right=[N*(N-1)+1:N*N-1;N*(N-1)+2:N*N]';
up=[N:N:N*(N-1);2*N:N:N*N]';
down=[1:N:N*(N-2)+1;1+N:N:N*(N-1)+1]';
% 组装各类边界
Db=up; % 平衡边界
Nb_L=left;Nb_R=right; %应力边界
Cb=down;
end
 