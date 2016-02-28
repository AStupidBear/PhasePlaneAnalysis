function [Xs,stable]=pplane(Fx,Fy,tRange, xRange, yRange,quiverScale) 
% X0=[x0 y0]--size(X0,1) stable points
% params--parameters
% t--time span
% x,y--grids coordinates in x and y directions
%% set parameters

    function y=Fun(x)
        y=[Fx(x(1),x(2));...
            Fy(x(1),x(2))];
    end
    function nullclines(x,y)
        figure(1); %figure handle
        ezplot(Fx,[min(x) max(x) min(y) max(y)]);
        hold on;
        ezplot(Fy,[min(x) max(x) min(y) max(y)]);
        title('phase plane analysis');
    end
    function dy=fun_ode(t,y)
        dy=zeros(2,1);
        dy(1)=Fx(y(1),y(2));
        dy(2)=Fy(y(1),y(2));
    end
 function stable=classify_fxd_point(eigvals)  % Classify fixed point based on eigenvalues of Jacobian.
        lam1=eigvals(1);
        lam2=eigvals(2);
        if ((isreal(lam1) && isreal(lam2)) && (lam1~=lam2))
            if(lam1 < 0 && lam2 < 0)
                stable = 'stable node';
            elseif (lam1 > 0 && lam2 > 0)
                stable = 'unstable node';
            elseif ((lam1 < 0 && lam2 > 0) || (lam1 > 0 && lam2 < 0))
                stable = 'saddle point';
            end;
        elseif(imag(lam1)~=0)
            if(real(lam1) < 0)
                stable = 'stable spiral';
            elseif(real(lam1) > 0)
                stable = 'unstable';
            end;
        else
            stable = 'center or spiral';
        end;
    end
%% Determine nullclines for the system
nullclines(xRange,yRange);
%% Calculate intersections using fsolve
options = optimoptions('fsolve','Display','off');
T=20;
Xs=zeros(T,2);
exitflag=zeros(T,1);
Jacobian=cell(T,1);
for ii=1:T
    init(1,1)=min(xRange)+max(xRange)-min(xRange)*rand;
    init(1,2)=min(yRange)+max(yRange)-min(yRange)*rand;%randomly pick points in the given range and solve equations
[Xs(ii,:),~,exitflag(ii,1),~,Jacobian{ii}]=fsolve(@Fun,init,options);
end
Xs(exitflag~=1)=[];
Jacobian(exitflag~=1)=[];
Xs=round(Xs*10000)/10000;
[Xs,ind,~]=unique(Xs,'rows'); %select unique steady points
Jacobian=Jacobian(ind);
%% Determine the vector field for the phase plane.
[X, Y] = meshgrid(xRange, yRange);
Y=flipud(Y);
U=Fx(X,Y);
V=Fy(X,Y);
figure(1);
quiver(X,Y,U,V,quiverScale);
%% Find trajectories and plot them. Use 4th order Runge-Kutta to simulate system.
    function init=init_point(~,~,t)
        
        pt=get(gca,'CurrentPoint');
        init=[pt(1,1) pt(1,2)];
        [~,Y]=ode45(@fun_ode,[min(t) max(t)],init);
        figure(1);
        plot(Y(:,1),Y(:,2));
    end
set(gcf,'WindowButtonDownFcn',{@init_point,tRange});
%% Determine the stability of one of the fixed points.
if(isempty(Xs))
    stable='no fixed points';
else
    if(size(Xs,1) ~= 1)
        display('system has multiple fixed points');
    end
    stable=cell(size(Xs,1),1);
    %Calculate Jacobian of system in question.
    for k=1:size(Xs,1)
        eigvals = eig(Jacobian{k});
        stable{k}=classify_fxd_point(eigvals);
    end
end
end