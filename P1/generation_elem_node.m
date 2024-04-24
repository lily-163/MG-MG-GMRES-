function [node,elem,Db,Nb_L,Nb_R,Cb]=generation_elem_node(h)
% Omega��ʼ�߽缰��Ԫ
% 3----6----9
% |\ 4 |\ 8 |
% | 3\ | 7\ |
% 2----5----8
% |\ 2 |\ 6 |
% | 1\ | 5\ |
% 1----4----7  
if 1/h ~= round(1/h)
    error('����Ļ��ֲ���,�������޸�\n');
end
N=1/h+1; % ÿ���ߵĽڵ���
nt=N-1;  % ÿ���ߵĵ�Ԫ��
elem=zeros(nt^2*2,3);
for i=1:1:nt
    start=(i-1)*N;numt=(i-1)*2*nt;
    % [2 1 4];      % ��������
    elem(numt+1:2:numt+2*nt,:)=[(start+2:start+nt+1)',(start+1:start+nt)',(start+1+N:start+N+nt)'];
    % [4 5 2]
    elem(numt+2:2:numt+2*nt,:)=[(start+1+N:start+N+nt)',(start+2+N:start+N+nt+1)',(start+2:start+nt+1)'];
end

node=zeros(N^2,2);
a=[0:h:1]';
b=ones(N,1);
node(:,1)=kron(a,b);
node(:,2)=kron(b,a);

% ���ɱ߽��
Db=[];Nb_L=[];Nb_R=[];Cb=[];
left=[1:N-1;2:N]';
right=[N*(N-1)+1:N*N-1;N*(N-1)+2:N*N]';
up=[N:N:N*(N-1);2*N:N:N*N]';
down=[1:N:N*(N-2)+1;1+N:N:N*(N-1)+1]';
% ��װ����߽�
Db=up; % ƽ��߽�
Nb_L=left;Nb_R=right; %Ӧ���߽�
Cb=down;
end
 