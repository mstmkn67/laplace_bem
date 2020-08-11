phi=$output.face[].phi
max_phi,min_phi=max(phi),min(phi)

face=[]
region=$input.region_condition[]
for i in region:
	if i[0]!="zp":
		face+=i[1]

for i in face:
	index=Location(getLocation("Face",i)).getIndex()[0]
	v=$input.face[index].vertex[]
	r0=get(getLocation("Vertex",v[0])+".position")
	r1=get(getLocation("Vertex",v[1])+".position")
	r2=get(getLocation("Vertex",v[2])+".position")
	c=(phi[index]-min_phi)/(max_phi-min_phi)
	polygon([r0,r1,r2],[c,0,1-c,1])
