phi=$output.edge[].phi
edge=$input.edge[].vertex[]
phi_min,phi_max=min(phi),max(phi)
E=$input.external_field.gradient_phi
a=1.0
for i in range(len(edge)):
	r0=get(getLocation("Vertex",edge[i][0])+".position")
	r1=get(getLocation("Vertex",edge[i][1])+".position")
	#show simulation result
	c=(phi[i]-phi_min)/(phi_max-phi_min)
	line([r0[0],r0[1],0.0],[r1[0],r1[1],0.0],[c,0,1.-c,1])
	#show analitical solution
	r=[0.5*(r0[0]+r1[0]),0.5*(r0[1]+r1[1]),0.0]
	r2=r[0]*r[0]+r[1]*r[1]
	c=(1.0+a*a/r2)*(r[0]*E[0]+r[1]*E[1])
	print(phi[i],c)
	c=(phi[i]+2*a)/(4*a)
	line([r0[0]+2.5*a,r0[1],0.0],[r1[0]+2.5*a,r1[1],0.0],[c,0,1.-c,1])
text([-0.5*a,-1.5*a,0],"simulation",[0,0,0,1,36])
text([2.5*a-0.5*a,-1.5*a,0],"analytical solution",[0,0,0,1,36])
