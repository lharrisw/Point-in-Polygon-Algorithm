clear
close all
clc

% parameterized curve (can be any parameterized curve)
N = 900;
tt = linspace(0,2*pi,N);
dt = tt(2)-tt(1);
x = 1.5 + 3*sin(tt) + 5*cos(2*tt);
y = 7*cos(tt);
s = [x;y];

% number of points in the polygon
aa = round(max([x y])); % grid size (aa x aa)
[xx,yy] = meshgrid(linspace(-aa,aa,2*aa+1)); % formation of the grid
[inside,on] = inpolygon(xx,yy,x,y); % matlab point-in-polygon algorithm

% % area enclosed by the curve
% xx2 = @(y) -7/2+3*sqrt(1-y.^2/49)+10*y.^2/49;
% xx1 = @(y) -7/2-3*sqrt(1-y.^2/49)+10*y.^2/49;
% area = integral(@(y)xx2(y)-xx1(y),-7,7);

% unit tangent and outward unit normal
dxdt = 3*cos(tt)-10*sin(2*tt); % derivative of x-component
dydt = -7*sin(tt); % derivative of y-component
t = [dxdt;dydt]... % unit tangent
     ./vecnorm([dxdt;dydt]);
n = [-dydt;dxdt]./vecnorm([dxdt;dydt]); % unit normal

% vector field
ds = sqrt(dxdt.^2+dydt.^2)*dt; % arc length
in = 0; % initialization of inside points
out = 0; % initilization of outside points
Flux = zeros(2*aa+1); % matrix to store Flux values

% % debugging code
% a = -1;
% b = 0;
% q = 1;
% E1x = q*(x-a)./((x-a).^2+(y-b).^2)/(2*pi);
% E1y = q*(y-b)./((x-a).^2+(y-b).^2)/(2*pi);
% E1 = [E1x;E1y];
% Flux1 = dot(E1,n)*ds';

% iterate over each grid point and compute the flux through the surface
for j = 1:2*aa+1
    for k = 1:2*aa+1
        q = 2*pi; % define charge with value of 2pi
        E = q.*(s-[xx(1,k);yy(j,1)])... % E field eval. at s when q = (xx(1,j),yy(k,1))
            ./(2*pi*vecnorm(s-[xx(1,k);yy(j,1)]).^2);
        Flux(j,k)= dot(E,n)*ds'; % Flux of E 2*pi*winding number
        if abs(Flux(j,k)) > q/2 
            in = in + 1; % flux inside = q
        else
            out = out + 1; % flux outside = 0
        end
    end
end

disp(['Points inside the curve (algorithm) = ',num2str(in)]);
disp(['Points inside the curve (MATLAB) = ',num2str(numel(yy(inside)))]);
disp(['Points outside the curve = ',num2str(out)]);
 
figure
hold on
grid on
plot(x,y); % plot of the curve
plot(xx,yy,'.'); % plot of the vector field
% plot(xx1(yy),yy,'ro'); % plot of left half of curve
% plot(xx2(yy),yy,'ro'); % plot of right half of curve
% quiver(x,y,n(1,:),n(2,:)); % unit normal
% quiver(x,y,t(1,:),t(2,:)); % unit tangent
% quiver(x,y,E1(1,:),E1(2,:),'color','k'); % electric field evaluated on the curve
axis equal
