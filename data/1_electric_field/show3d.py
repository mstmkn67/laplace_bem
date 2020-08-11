from math import *
phi=$output.face[].phi
face=$input.face[].vertex[]
phi_min,phi_max=min(phi),max(phi)
E=$input.external_field.gradient_phi
a=1.0
for i in range(len(face)):
	r0=get(getLocation("Vertex",face[i][0])+".position")
	r1=get(getLocation("Vertex",face[i][1])+".position")
	r2=get(getLocation("Vertex",face[i][2])+".position")
	#show simulation result
	c=(phi[i]-phi_min)/(phi_max-phi_min)
	polygon([r0,r1,r2],[c,0,1.-c,1])
	#show analitical solution
	r=[(r0[0]+r1[0]+r2[0])/3.0,(r0[1]+r1[1]+r2[1])/3.0,(r0[2]+r1[2]+r2[2])/3.0]
	l2=r[0]*r[0]+r[1]*r[1]+r[2]*r[2]
	c=(1.0+0.5*a*a*a/l2/sqrt(l2))*(r[0]*E[0]+r[1]*E[1]+r[2]*E[2])
	print(phi[i],c)
	c=(phi[i]+1.5*a)/(3*a)
	r0[0]+=2.5*a;r1[0]+=2.5*a;r2[0]+=2.5*a
	polygon([r0,r1,r2],[c,0,1.-c,1])
text([-0.5*a,-1.5*a,0],"simulation",[0,0,0,1,36])
text([2.5*a-0.5*a,-1.5*a,0],"analytical solution",[0,0,0,1,36])
