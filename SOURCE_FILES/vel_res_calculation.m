allall=(sum(abs(Vx))+sum(abs(Vy)))/(length(Vx)+length(Vy));

if timestep>1
    pfpf=allall-allall0;
    maxdiffVx=max(abs(Vx-Vx0));
    maxdiffVy=max(abs(Vy-Vy0));    
    iterations(1,niter+maxiter*(timestep-1))=abs(pfpf);
    pfpf=max([maxdiffVx maxdiffVy]);
    iterations(2,niter+maxiter*(timestep-1))=pfpf;
end

allall0=allall;
Vx0=Vx;
Vy0=Vy;

if timestep==1
    P0 = P;
end

dP_node(1,niter+maxiter*(timestep-1)) = max(abs(P-P0));

P0=P;