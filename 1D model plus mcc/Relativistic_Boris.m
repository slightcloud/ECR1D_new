function v_new = Relativistic_Boris(m,q,E,B,v_old,dt)
global c
gamma_old = sqrt(1./(1-sum((v_old'/c).^2)));
u_old = gamma_old'.*v_old;
u_minus = u_old+q*E*dt/(2*m);
gamma_middle = sqrt(1+sum((u_minus'/c).^2));
t = q*B*dt./(2*gamma_middle'*m);
u_prime = u_minus + cross(u_minus,t);
s = 2*t./(1+sum(t'.^2))';
u_plus = u_minus + cross(u_prime,s);
u_new = u_plus +q*dt*E/(2*m);
gamma_new = sqrt(1+sum((u_new/c)'.^2));
v_new = u_new./gamma_new';