tic
world = World;
world=world.setLimits(10,5);
world=world.setNodes(101,51);
ni = world.ni;
nj = world.nj;
dx = world.dx;
dv = world.dv;
dt = 1/2.0;

disp(['dx: ',num2str(dx),'dv: ',num2str(dv)]);

world.periodic = true;


f = zeros(world.ni,world.nj);
fs = zeros(world.ni,world.nj);
fss = zeros(world.ni,world.nj);
ne = zeros(ni,1);
b = zeros(ni,1);
E = zeros(ni,1);
phi = zeros(ni,1);
for i = 1:ni
    for j = 1:nj
        x = world.getX(i);
        v = world.getV(j);

        vth2 = 0.02;
        vs1 = 1.6;
        vs2 = -1.4;
        f(i,j) = 0.5/sqrt(vth2*pi)*exp(-(v-vs1)*(v-vs1)/vth2);
        f(i,j) = f(i,j)+0.5/sqrt(vth2*pi)*exp(-(v-vs2)*(v-vs2)/vth2)*(1+0.02*cos(3*pi*x/(world.L)));
    end
end

%set some constant e field
for i = 1:ni
    E(i) = 0;
end

for it = 0:1000
    if (mod(it,100)==0)
        disp(it)
    end
%     if (mod(it,50)==0)
        %             world.saveVTK(it,world,scalar2D,scatter1D)
%         figure;
        imagesc(f')
        drawnow
%     end

    %compute f*
    for i = 1:ni
        for j = 1:nj
            v = world.getV(j);
            x = world.getX(i);

            fs(i,j) = world.interp(f,x-v*0.5*dt,v);
        end
    end

    fs = world.applyBC(fs);

    %compute number density by integrating f with the trapezoidal rule
    for i = 1:ni
        ne(i) = 0;
        for j = 1:nj-1
            ne(i) = ne(i) + 0.5*(fs(i,j+1)+fs(i,j))*dv;
        end
    end

    %compute the right hand side, -rho = (ne-1)
    for i = 1:ni
        b(i) = ne(i)-1;
    end
    b(1) = 0.5*(b(1)+b(ni));
    b(ni) = b(1);

    %solution of the Poisson's equation
    [b,phi,E] = solvePoissonsEquationGS(world,b,phi,E);

    %compute f**
    for i = 1:ni
        for j = 1:nj
            v = world.getV(j);
            x = world.getX(i);
            fss(i,j) = world.interp(fs,x,v+E(i)*dt);
        end
    end

    fss = world.applyBC(fss);

    for i = 1:ni
        for j = 1:nj
            v = world.getV(j);
            x = world.getX(i);
            f(i,j) = world.interp(fss,x-v*0.5*dt,v);
        end
    end

    f = world.applyBC(f);
end
toc

%     saveVTK(it,world,scalars2D,scarlars1D);

% solves Poisson's equation with Dirichlet boundaries using the direct Thomas algorithm and returns the electric field
function [b,phi,E] = solvePoissonsEquationGS(obj,b,phi,E)
dx2 = obj.dx*obj.dx;
dx = obj.dx;
ni = obj.ni;
tol = 1e-3;
for i = 1:ni
    phi(i) = 0;
end
for it = 1:10000
    phi(1) = 0.5*(phi(ni-1)+phi(2)-dx2*b(1));
    for i = 2:ni-1
        g = 0.5*(phi(i-1)+phi(i+1)-dx2*b(i));
        phi(i) = phi(i)+1.4*(g-phi(i));
    end
    phi(ni) = 0.5*(phi(ni-1)+phi(2)-dx2*b(ni));

    % check for convergence
    if (mod(it,50)==0)
        R_sum = 0;
        for i =2:ni-1
            dR = (phi(i-1)-2.*phi(i)+phi(i+1))/dx2-b(i);
            R_sum = R_sum+dR*dR;
        end
        dR = (phi(ni-1)-2*phi(1)+phi(2))/dx2-b(1);
        R_sum = R_sum+dR*dR;
        dR = (phi(ni-1)-2*phi(ni)+phi(2))/dx2-b(ni);
        R_sum = R_sum+dR*dR;
        norm = sqrt(R_sum/ni);
        if(norm<tol)
            break;
        end
    end
end

if(norm>tol)
    disp(['GS failed to converge, norm = ',num2str(norm)])
else
    disp(['OK, norm = ',num2str(norm)])
end

%set periodic boundary
phi(1) = 0.5*(phi(1)+phi(ni));
phi(ni) = phi(1);

%compute electric field
for i = 2:obj.ni-1
    E(i) = -(phi(i+1)-phi(i-1))/(2*dx);
end
E(1) = -(phi(2)-phi(obj.ni-1))/(2*dx);
E(obj.ni) = E(1);
end


