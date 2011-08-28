clear
close all
M  = [80.81210000000002,2.32343142623549;2.32343142623549,0.30126570934256];
K0 = [-794.119500000000,-25.739089291258;-25.739089291258,-8.139414705882];
K2 = [0,76.40620875965657;0,2.67560553633218];
C1 = [0,33.77386947593010;-0.84823447825693,1.70696539792387];
eigval = zeros(4,1000);
n=1000;
vmax = 20;
for i=1:n  
    v(i) = (i-1)/n*vmax;
    C = C1.*v(i);
    K = K0 + K2.*v(i)^2;
    stab = zeros(4);
    stab(1,3)=1;
    stab(2,4)=1;
    stab(3,1)=(K(2,1)/M(2,2)-K(1,1)/M(1,2))/(M(1,1)/M(1,2)-M(2,1)/M(2,2));
    stab(3,2)=(K(2,2)/M(2,2)-K(1,2)/M(1,2))/(M(1,1)/M(1,2)-M(2,1)/M(2,2));
    stab(3,3)=(C(2,1)/M(2,2)-C(1,1)/M(1,2))/(M(1,1)/M(1,2)-M(2,1)/M(2,2));
    stab(3,4)=(C(2,2)/M(2,2)-C(1,2)/M(1,2))/(M(1,1)/M(1,2)-M(2,1)/M(2,2));
    stab(4,1)=(K(2,1)/M(2,1)-K(1,1)/M(1,1))/(M(1,2)/M(1,1)-M(2,2)/M(2,1));
    stab(4,2)=(K(2,2)/M(2,1)-K(1,2)/M(1,1))/(M(1,2)/M(1,1)-M(2,2)/M(2,1));
    stab(4,3)=(C(2,1)/M(2,1)-C(1,1)/M(1,1))/(M(1,2)/M(1,1)-M(2,2)/M(2,1));
    stab(4,4)=(C(2,2)/M(2,1)-C(1,2)/M(1,1))/(M(1,2)/M(1,1)-M(2,2)/M(2,1));
    eigval(1:4,i)=eig(stab);
end
%-------------------Organize Eigenvalue Matrix
% set first column in organized matrix to the first column of the original matrix
eigorg(:,1)=eigval(:,1);
% rearrange the eigenvalue columns by calcuating the absolute value of the
% difference between values of the successicve column and the preceeding
% column
for i = 1:length(v)-1
    for j = 1:4
        first = eigorg(j,i);
        for k =1:4
            second = eigval(k,i+1);
            diff(k) = abs(second - first);
        end
        [min indice]=min(diff);
        eigorg(j,i+1)=eigval(indice,i+1);
    end
end
%-------------------Extract Weave, Capsize, and Caster Eigenmodes
found = 0;
for i = 1:size(eigorg,1)
    if isreal(eigorg(i,:)) == 0 & found==0
        weave(1,:)=real(eigorg(i,:));
        found = 1;
    elseif isreal(eigorg(i,:)) == 0 & found==1
        weave(2,:)=real(eigorg(i,:));
    elseif isreal(eigorg(i,:)) == 1 & (eigorg(i,size(eigorg,2))-eigorg(i,1))>0
        capsize=real(eigorg(i,:));
    else
        caster=real(eigorg(i,:));
    end
end
%-------------------Plot Eigenmodes vs Velocity
figure(3)
hold on
axis([0 10 -10 10])
title('Eigenmodes')
xlabel('Velocity [m/s]')
ylabel('Eigenvalues [1/s]')
plot(v,weave(1,:),'.b')
plot(v,weave(2,:),'.b')
plot(v,capsize,'.y')
plot(v,caster,'.r')
plot(v,zeros(length(v),1),'k') %plot horizontal line at zero
hold off
%-------------------Calculate Critical Velocities
% find the velocity at which the weave mode becomes stable
for i=1:length(weave(1,:))
    if weave(1,i)<=0
        index=i;
        break
    end
end
vw=mean([v(index-1),v(index)]) % weave critical velocity
% find the velocity at which the capsize mode becomes unstable
for i=1:length(capsize)
    if capsize(i)>=0
        index=i;
        break
    end
end
vc=mean([v(index-1),v(index)]) % capsize critical velocity
%-------------------Plot the Critical Velocities
hold on
plot(vw,0,'ok',vw,caster(length(caster)):0.1:0,'k',vc,0,'ok',vc,caster(length(caster)):0.1:0,'k')
hold off