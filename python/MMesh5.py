from numpy import cross,array,dot,sqrt,cos,sin,pi,abs,zeros,append
from sys import float_info
from copy import copy

#EPS=float_info.epsilon
EPS=1e-10

##tools for meshing
def on_sphere(vertex_):
	new_vertex=[]
	for i in vertex_:
		new_vertex.append(i/sqrt(i[0]*i[0]+i[1]*i[1]+i[2]*i[2]))
	return new_vertex

def on_circle(vertex_):
	new_vertex=[]
	for i in vertex_:
		s=1
		if i[2]!=0:
			s=i[2]/sqrt(i[0]*i[0]+i[1]*i[1])
		new_vertex.append([s*i[0],s*i[1],i[2]])
	return array(new_vertex)

def get_1st_generation_octahedron():
	vertex_=[[-1,-1,0],[ 1,-1,0],[0,0, 1],
	         [ 1,-1,0],[ 1, 1,0],[0,0, 1],
	         [ 1, 1,0],[-1, 1,0],[0,0, 1],
	         [-1, 1,0],[-1,-1,0],[0,0, 1],
	         [-1,-1,0],[0,0,-1],[ 1,-1,0],
	         [ 1,-1,0],[0,0,-1],[ 1, 1,0],
	         [ 1, 1,0],[0,0,-1],[-1, 1,0],
	         [-1, 1,0],[0,0,-1],[-1,-1,0]]
	return array(vertex_)

def get_1st_generation_icosahedron():
	phi=0.5*(1.0+sqrt(5.0))
	pos=[[0,-1,phi],[0,1,phi],[0,-1,-phi],[0,1,-phi],
	     [phi,0,1],[phi,0,-1],[-phi,0,1],[-phi,0,-1],
	     [1,phi,0],[1,-phi,0],[-1,phi,0],[-1,-phi,0]]
	vertex_=[copy(pos[0]),copy(pos[4]),copy(pos[1]),
	         copy(pos[0]),copy(pos[9]),copy(pos[4]),
	         copy(pos[9]),copy(pos[5]),copy(pos[4]),
	         copy(pos[4]),copy(pos[5]),copy(pos[8]),
	         copy(pos[4]),copy(pos[8]),copy(pos[1]),
	         copy(pos[8]),copy(pos[10]),copy(pos[1]),
	         copy(pos[8]),copy(pos[3]),copy(pos[10]),
	         copy(pos[5]),copy(pos[3]),copy(pos[8]),
	         copy(pos[5]),copy(pos[2]),copy(pos[3]),
	         copy(pos[2]),copy(pos[7]),copy(pos[3]),
	         copy(pos[7]),copy(pos[10]),copy(pos[3]),
	         copy(pos[7]),copy(pos[6]),copy(pos[10]),
	         copy(pos[7]),copy(pos[11]),copy(pos[6]),
	         copy(pos[11]),copy(pos[0]),copy(pos[6]),
	         copy(pos[0]),copy(pos[1]),copy(pos[6]),
	         copy(pos[6]),copy(pos[1]),copy(pos[10]),
	         copy(pos[9]),copy(pos[0]),copy(pos[11]),
	         copy(pos[9]),copy(pos[11]),copy(pos[2]),
	         copy(pos[9]),copy(pos[2]),copy(pos[5]),
	         copy(pos[7]),copy(pos[2]),copy(pos[11])]
	return array(vertex_)

def get_1st_generation_2d():
	vertex_=[[ 1,-1,1],[-1,-1,1],[0,0,0],
	         [ 1, 1,1],[ 1,-1,1],[0,0,0],
	         [-1, 1,1],[ 1, 1,1],[0,0,0],
	         [-1,-1,1],[-1, 1,1],[0,0,0]]
	return array(vertex_)

def get_1st_generation_2d_2():
	vertex_=[[ 0, 1,1],[ 1, 0,1],[0,0,0],
	         [ 1, 0,1],[ 0,-1,1],[0,0,0],
	         [ 0,-1,1],[-1, 0,1],[0,0,0],
	         [-1, 0,1],[ 0, 1,1],[0,0,0]]
	return array(vertex_)

def get_1st_generation_2d_2b():
	vertex_=[[ 1, 0,1],[ 0, 1,1],[0,0,0],
	         [ 0,-1,1],[ 1, 0,1],[0,0,0],
	         [-1, 0,1],[ 0,-1,1],[0,0,0],
	         [ 0, 1,1],[-1, 0,1],[0,0,0]]
	return array(vertex_)

def get_1st_generation_2d_hexagon():
	vertex_=[[cos(pi/3),sin(pi/3), 1],[ 1, 0, 1],[0,0,0],
	         [cos(2*pi/3),sin(2*pi/3), 1],[cos(pi/3),sin(pi/3), 1],[0,0,0],
	         [-1, 0,1],[cos(2*pi/3),sin(2*pi/3), 1],[0,0,0],
	         [cos(4*pi/3),sin(4*pi/3), 1],[-1, 0, 1],[0,0,0],
	         [cos(5*pi/3),sin(5*pi/3), 1],[cos(4*pi/3),sin(4*pi/3), 1],[0,0,0],
	         [1,0,1],[cos(5*pi/3),sin(5*pi/3), 1],[0,0,0]]
	return array(vertex_)

def get_1st_generation_box():
	vertex_=[[-1,-1,-1],[ 0, 0,-1],[ 1,-1,-1],
	         [ 1,-1,-1],[ 0, 0,-1],[ 1, 1,-1],
	         [ 1, 1,-1],[ 0, 0,-1],[-1, 1,-1],
	         [-1, 1,-1],[ 0, 0,-1],[-1,-1,-1],#
	         [ 1,-1, 1],[ 0, 0, 1],[-1,-1, 1],
	         [ 1, 1, 1],[ 0, 0, 1],[ 1,-1, 1],
	         [-1, 1, 1],[ 0, 0, 1],[ 1, 1, 1],
	         [-1,-1, 1],[ 0, 0, 1],[-1, 1, 1],#perpendicular to z axis
	         [-1, 1,-1],[ 0, 1, 0],[ 1, 1,-1],
	         [ 1, 1,-1],[ 0, 1, 0],[ 1, 1, 1],
	         [ 1, 1, 1],[ 0, 1, 0],[-1, 1, 1],
	         [-1, 1, 1],[ 0, 1, 0],[-1, 1,-1],#
	         [ 1,-1,-1],[ 0,-1, 0],[-1,-1,-1],
			 [-1,-1, 1],[ 0,-1, 0],[ 1,-1, 1],
	         [ 1,-1, 1],[ 0,-1, 0],[ 1,-1,-1],
	         [-1,-1,-1],[ 0,-1, 0],[-1,-1, 1],#perpendicular to y axis
	         [-1,-1,-1],[-1, 0, 0],[-1, 1,-1],
	         [-1, 1,-1],[-1, 0, 0],[-1, 1, 1],
	         [-1, 1, 1],[-1, 0, 0],[-1,-1, 1],
	         [-1,-1, 1],[-1, 0, 0],[-1,-1,-1],#
	         [ 1, 1,-1],[ 1, 0, 0],[ 1,-1,-1],
	         [ 1, 1, 1],[ 1, 0, 0],[ 1, 1,-1],
	         [ 1,-1, 1],[ 1, 0, 0],[ 1, 1, 1],
	         [ 1,-1,-1],[ 1, 0, 0],[ 1,-1, 1]]#perpendicular to x axis
	return array(vertex_)

def get_1st_generation_square_pyramid_without_bottom():
	vertex_=[[-1,-1, 0],[ 1,-1, 0],[0,0, 1],
	         [ 1,-1, 0],[ 1, 1, 0],[0,0, 1],
	         [ 1, 1, 0],[-1, 1, 0],[0,0, 1],
	         [-1, 1, 0],[-1,-1, 0],[0,0, 1]]
	return array(vertex_)

def create_next_generation(vertex_):
	new_vertex_=[]
	for i in range(0,len(vertex_),3):
		i0,i1,i2=i,i+1,i+2
		r3,r4,r5=0.5*(vertex_[i0]+vertex_[i1]),0.5*(vertex_[i1]+vertex_[i2]),0.5*(vertex_[i2]+vertex_[i0])
		for j in vertex_[i0],r3,r5, r3,vertex_[i1],r4, r4,vertex_[i2],r5, r5,r3,r4:
			new_vertex_.append(copy(j))
	return new_vertex_

###########
def compare_vertex_(v0,v1):
	if v1[0]-v0[0]>EPS:
		return -1
	elif v0[0]-v1[0]>EPS:
		return 1
	else:
		if v1[1]-v0[1]>EPS:
			return -1
		elif v0[1]-v1[1]>EPS:
			return 1
		else:
			if v1[2]-v0[2]>EPS:
				return -1
			elif v0[2]-v1[2]>EPS:
				return 1
			else:
				return 0

def sort_vertex_(l,r,vertex_):
	if l<r:
		mv=copy(vertex_[int((l+r)/2)])
		i,j=l-1,r+1
		while True:
			while True:
				i+=1
				if compare_vertex_(vertex_[i],mv)!=-1:
					break
			while True:
				j-=1
				if compare_vertex_(vertex_[j],mv)!=1:
					break
			if i>=j:
				break
			swap=copy(vertex_[i])
			vertex_[i]=copy(vertex_[j])
			vertex_[j]=copy(swap)
			#swap=vertex_[i]
			#vertex_[i]=vertex_[j]
			#vertex_[j]=swap
		sort_vertex_(l,i-1,vertex_)
		sort_vertex_(j+1,r,vertex_)

def number_vertex_(vertex_):
	new_vertex_=[]
	for i in range(len(vertex_)):
		new_vertex_.append(append(vertex_[i],i))
	return array(new_vertex_)

def reduce_vertex(vertex_):
	vertex_=number_vertex_(vertex_)
	sort_vertex_(0,len(vertex_)-1,vertex_)
	n=len(vertex_)
	vertex,face=[],[[-1,-1,-1] for i in range(int(n/3))]
	refv=[float_info.max,float_info.max,float_info.max]
	start,end=0,0
	for i in range(n):
		if (i==n-1 or 
		   abs(vertex_[i][0]-refv[0])>EPS or
		   abs(vertex_[i][1]-refv[1])>EPS or
		   abs(vertex_[i][2]-refv[2])>EPS):
			end=i
			if i!=0:
				vertex.append(copy(refv))
				for j in range(start,end+1):
					fi=int(vertex_[j][3]/3)
					vi=int(vertex_[j][3]%3)
					face[fi][vi]=len(vertex)-1
			if end!=n-1:
				refv=vertex_[end][0:3]
				start=end
	return vertex,face

#def get_edge_old(face):############ 
#	edge=[]
#	for i in range(len(face)):
#		t=[[face[i][0],face[i][1]],[face[i][1],face[i][2]],[face[i][2],face[i][0]]]
#		for j in t:
#			flag=0
#			for k in range(len(edge)):
#				if (j[0]==edge[k][0][0] and j[1]==edge[k][0][1]) or (j[0]==edge[k][0][1] and j[1]==edge[k][0][0]):
#					edge[k][1].append(i)
#					flag=1
#					continue
#			if flag==0:
#				edge.append([j,[i]])
#	return edge

def compare_edge_(e0,e1):
	if e1[0]>e0[0]:
		return -1
	elif e0[0]>e1[0]:
		return 1
	else:
		if e1[1]>e0[1]:
			return -1
		elif e0[1]>e1[1]:
			return 1
		else:
			return 0

def sort_edge_(l,r,edge_):
	if l<r:
		me=copy(edge_[(l+r)//2])
		i,j=l-1,r+1
		while True:
			while True:
				i+=1
				if compare_edge_(edge_[i],me)!=-1:
					break
			while True:
				j-=1
				if compare_edge_(edge_[j],me)!=1:
					break
			if i>=j:
				break
			swap=copy(edge_[i])
			edge_[i]=copy(edge_[j])
			edge_[j]=copy(swap)
		sort_edge_(l,i-1,edge_)
		sort_edge_(j+1,r,edge_)

def get_edge(face):
	edge_=[]
	for i in face:
		for j in range(3):
			if i[j%3] > i[(j+1)%3]:
				edge_.append([i[j%3],i[(j+1)%3]])
			else:
				edge_.append([i[(j+1)%3],i[j%3]])
	sort_edge_(0,len(edge_)-1,edge_)
	n=len(edge_)
	refe=edge_[0]
	edge=[]
	for i in range(n):
		if (i==n-1 or 
		   edge_[i][0]!=refe[0] or edge_[i][1]!=refe[1]):
			edge.append(copy(refe))
			refe=edge_[i]
	return edge

def euler_number(vertex,face):
	edge=get_edge(face)
	e,v,f=len(edge),len(vertex),len(face)
	return v-e+f

def draw_face(vertex,face,udf):
	for i in face:
		r0,r1,r2=vertex[i[0]],vertex[i[1]],vertex[i[2]]
		r=[[r0[0],r0[1],r0[2]],[r1[0],r1[1],r1[2]],[r2[0],r2[1],r2[2]]]
		udf.polygon(r,[0,0,1,1])

def draw_face_with_normal(vertex,face,udf):
	for i in face:
		r0,r1,r2=vertex[i[0]],vertex[i[1]],vertex[i[2]]
		r=[[r0[0],r0[1],r0[2]],[r1[0],r1[1],r1[2]],[r2[0],r2[1],r2[2]]]
		udf.polygon(r,[0,0,1,1])
		r=array(r)
		n=array(cross(r[1]-r[0],r[2]-r[0]))
		n=n/sqrt(dot(n,n))
		r0=(r[0]+r[1]+r[2])/3.0
		udf.line(r0.tolist(),(r0+n).tolist(),1)
#############################
def sphere(radius,discretization_level,mother="octahedron"):
	if mother=="octahedron":
		vertex_=get_1st_generation_octahedron()
	elif mother=="icosahedron":
		vertex_=get_1st_generation_icosahedron()
	else:
		print("error")
		exit(1)
	vertex_=on_sphere(vertex_)
	for i in range(discretization_level):
		vertex_=create_next_generation(vertex_)
		vertex_=on_sphere(vertex_)
	vertex,face=reduce_vertex(vertex_)
	for i in range(len(vertex)):
		r=vertex[i]
		vertex[i]=[radius*r[0],radius*r[1],radius*r[2]]
	return array(vertex),face

def hemisphere_without_bottom(radius,discretization_level):
	vertex_=get_1st_generation_square_pyramid_without_bottom()
	vertex_=on_sphere(vertex_)
	for i in range(discretization_level):
		vertex_=create_next_generation(vertex_)
		vertex_=on_sphere(vertex_)
	vertex,face=reduce_vertex(vertex_)
	for i in range(len(vertex)):
		r=vertex[i]
		vertex[i]=[radius*r[0],radius*r[1],radius*r[2]]
	return array(vertex),face

def hemi_ellipsoid_without_bottom(a,b,c,discretization_level):
	vertex,face=hemisphere_without_bottom(1,discretization_level)
	for i in range(len(vertex)):
		r=vertex[i]
		vertex[i]=[a*r[0],b*r[1],c*r[2]]
	return array(vertex),face

def hemisphere_with_bottom(radius,discretization_level):
	vertex_=get_1st_generation_square_pyramid_without_bottom()
	vertex_=on_sphere(vertex_)
	vertex2_=get_1st_generation_2d()
	vertex2_=on_circle(vertex2_)
	for i in range(discretization_level):
		vertex_=create_next_generation(vertex_)
		vertex_=on_sphere(vertex_)
		vertex2_=create_next_generation(vertex2_)
		vertex2_=on_circle(vertex2_)
	for i in range(len(vertex_)):
		r=vertex_[i]
		vertex_[i]=[radius*r[0],radius*r[1],radius*r[2]]
	for i in range(len(vertex2_)):
		r=vertex2_[i]
		vertex2_[i]=[radius*r[0],radius*r[1],0.0]
	vertex_=array(array(vertex_).tolist()+array(vertex2_).tolist())
	vertex,face=reduce_vertex(vertex_)
	return array(vertex),face

def hemi_ellipsoid_with_bottom(a,b,c,discretization_level):
	vertex,face=hemisphere_with_bottom(1,discretization_level)
	for i in range(len(vertex)):
		r=vertex[i]
		vertex[i]=[a*r[0],b*r[1],c*r[2]]
	return array(vertex),face

def ellipsoid(a,b,c,discretization_level,mother="octahedron"):
	vertex,face=sphere(1,discretization_level,mother)
	for i in range(len(vertex)):
		r=vertex[i]
		vertex[i]=[a*r[0],b*r[1],c*r[2]]
	return array(vertex),face

def biconcave(a,c1,c2,c3,discretization_level,mother="octahedron"):
	vertex,face=sphere(1,discretization_level,mother)
	for i in range(len(vertex)):
		r=vertex[i]
		x2z2=r[0]*r[0]+r[2]*r[2]
		if x2z2>1.0 or r[1]==0.0:
			r[1]=0.0
		else:
			p=r[1]/abs(r[1])
			r[1]=p*0.5*(c1+c2*x2z2+c3*x2z2*x2z2)*sqrt(1.0-x2z2)
		vertex[i]=[a*r[0],a*r[1],a*r[2]]
	return array(vertex),face

def rect3d(a,b,c,discretization_level):
	vertex_=get_1st_generation_box()
	for i in range(discretization_level):
		vertex_=create_next_generation(vertex_)
	vertex,face=reduce_vertex(vertex_)
	for i in range(len(vertex)):
		r=vertex[i]
		vertex[i]=[a*r[0],b*r[1],c*r[2]]
	return array(vertex),face

#radius a
#pitch 2 pi b
#length 2l
def twist_ribbon(a,b,l,Nx,Nz,ch=1):
	###rectangular lattice
	sx,sz=1.0/Nx,1.0/Nz
	r=[0]*((Nx+1)*(Nz+1))
	for j in range(Nz+1):
		for i in range(Nx+1):
			index=j*(Nx+1)+i
			x,z=2*(i*sx-0.5),2*(j*sz-0.5)
			r[index]=[x,0.0,z]
	ele=[]
	for j in range(Nz):
		for i in range(Nx):
			index=j*(Nx+1)+i
			ele.append([index,index+1,index+1+Nx+1])
			ele.append([index+Nx+1,index,index+1+Nx+1])
	#right_helicoid
	for i in range(len(r)):
		r[i]=[a*r[i][0]*cos(l*r[i][2]/b),ch*a*r[i][0]*sin(l*r[i][2]/b),l*r[i][2]]
	return array(r),ele

def disk(radius,discretizetion_level,mother="square"):
	if mother=="square":
		vertex_=get_1st_generation_2d()
	elif mother=="hexagon":
		vertex_=get_1st_generation_2d_hexagon()
	for i in range(discretizetion_level):
		vertex_=create_next_generation(vertex_)
	vertex_=on_circle(vertex_)
	vertex,face=reduce_vertex(vertex_)
	for i in range(len(vertex)):
		r=vertex[i]
		vertex[i]=[radius*r[0],radius*r[1],0.0]
	return array(vertex),face

def ellipse(a,b,discretizetion_level,mother="square"):
	vertex,face=disk(1,discretizetion_level,mother)
	for i in range(len(vertex)):
		r=vertex[i]
		vertex[i]=[a*r[0],b*r[1],0.0]
	return array(vertex),face

def rect2d(a,b,discretizetion_level):
	vertex_=get_1st_generation_2d()
	for i in range(discretizetion_level):
		vertex_=create_next_generation(vertex_)
	vertex,face=reduce_vertex(vertex_)
	for i in range(len(vertex)):
		r=vertex[i]
		vertex[i]=[a*r[0],b*r[1],0.0]
	return array(vertex),face

def cylinder(radius,half_length,n_theta,n_height):
	vertex=[]
	for i in range(n_height+1):
		z=-1.0+2.0*i/n_height
		for j in range(n_theta):
			th=2.*pi*j/n_theta
			vertex.append([radius*cos(th),radius*sin(th),half_length*z])
	vertex_=[]
	for i in range(n_height):
		for j in range(n_theta):
			vertex_.append(vertex[i*n_theta+j])
			vertex_.append(vertex[i*n_theta+(j+1)%n_theta])
			vertex_.append(vertex[(i+1)*n_theta+j])
			vertex_.append(vertex[(i+1)*n_theta+j])
			vertex_.append(vertex[i*n_theta+(j+1)%n_theta])
			vertex_.append(vertex[(i+1)*n_theta+(j+1)%n_theta])
	vertex,face=reduce_vertex(vertex_)
	return vertex,face

def cylinder_with_caps(radius,half_length,cap_discretizetion_level,n_height):
	n_theta=int(pow(2,cap_discretizetion_level+2))
	vertex=[]
	for i in range(n_height+1):
		z=-1.0+2.0*i/n_height
		for j in range(n_theta):
			th=2.*pi*j/n_theta
			vertex.append([radius*cos(th),radius*sin(th),half_length*z])
	vertex_=[]
	for i in range(n_height):
		for j in range(n_theta):
			vertex_.append(vertex[i*n_theta+j])
			vertex_.append(vertex[i*n_theta+(j+1)%n_theta])
			vertex_.append(vertex[(i+1)*n_theta+j])
			vertex_.append(vertex[(i+1)*n_theta+j])
			vertex_.append(vertex[i*n_theta+(j+1)%n_theta])
			vertex_.append(vertex[(i+1)*n_theta+(j+1)%n_theta])
	#cap
	vertex_cap0_=get_1st_generation_2d_2()
	vertex_cap1_=get_1st_generation_2d_2b()
	for i in range(cap_discretizetion_level):
		vertex_cap0_=create_next_generation(vertex_cap0_)
		vertex_cap1_=create_next_generation(vertex_cap1_)
		vertex_cap0_=on_circle(vertex_cap0_)
		vertex_cap1_=on_circle(vertex_cap1_)
	for i in range(len(vertex_cap0_)):
		vertex_cap0_[i][0]=radius*vertex_cap0_[i][0]
		vertex_cap0_[i][1]=radius*vertex_cap0_[i][1]
		vertex_cap0_[i][2]=-half_length
		vertex_cap1_[i][0]=radius*vertex_cap1_[i][0]
		vertex_cap1_[i][1]=radius*vertex_cap1_[i][1]
		vertex_cap1_[i][2]=half_length
	vertex,face=reduce_vertex(array(array(vertex_).tolist()+array(vertex_cap0_).tolist()+array(vertex_cap1_).tolist()))
	return vertex,face

def cylinder_with_spherical_caps(radius,half_length,cap_discretizetion_level,n_height):
	n_theta=int(pow(2,cap_discretizetion_level+2))
	vertex=[]
	for i in range(n_height+1):
		z=-1.0+2.0*i/n_height
		for j in range(n_theta):
			th=2.*pi*j/n_theta
			vertex.append([radius*cos(th),radius*sin(th),half_length*z])
	vertex_=[]
	for i in range(n_height):
		for j in range(n_theta):
			vertex_.append(vertex[i*n_theta+j])
			vertex_.append(vertex[i*n_theta+(j+1)%n_theta])
			vertex_.append(vertex[(i+1)*n_theta+j])
			vertex_.append(vertex[(i+1)*n_theta+j])
			vertex_.append(vertex[i*n_theta+(j+1)%n_theta])
			vertex_.append(vertex[(i+1)*n_theta+(j+1)%n_theta])
	#cap
	vertex_cap_=get_1st_generation_square_pyramid_without_bottom()
	vertex_cap_=on_sphere(vertex_cap_)
	for i in range(cap_discretizetion_level):
		vertex_cap_=create_next_generation(vertex_cap_)
		vertex_cap_=on_sphere(vertex_cap_)
	for i in range(len(vertex_cap_)):
		vertex_cap_[i][0]=radius*vertex_cap_[i][0]
		vertex_cap_[i][1]=radius*vertex_cap_[i][1]
		vertex_cap_[i][2]+=half_length
	vertex_cap2_=zeros((len(vertex_cap_),3),float)
	for i in range(len(vertex_cap2_)):
		vertex_cap2_[i][0]=vertex_cap_[i][1]
		vertex_cap2_[i][1]=vertex_cap_[i][0]
		vertex_cap2_[i][2]=-vertex_cap_[i][2]
	vertex,face=reduce_vertex(array(array(vertex_).tolist()+array(vertex_cap_).tolist()+array(vertex_cap2_).tolist()))
	return vertex,face

def cylinder_with_spherical_caps2(radius,half_length,cap_discretizetion_level,n_height):
	n_theta=int(pow(2,cap_discretizetion_level+2))
	vertex=[]
	for i in range(n_height+1):
		z=-1.0+2.0*i/n_height
		for j in range(n_theta):
			th=2.*pi*j/n_theta
			vertex.append([radius*cos(th),radius*sin(th),half_length*z])
	###
	for i in range(n_height):
		z=-1.0+2.0*(i+0.5)/n_height
		for j in range(n_theta):
			th=2.*pi*(j+0.5)/n_theta
			vertex.append([radius*cos(th),radius*sin(th),half_length*z])
	vertex_=[]
	nn=(n_height+1)*n_theta
	for i in range(n_height):
		for j in range(n_theta):
			r0,r1,r2=vertex[i*n_theta+j],vertex[(i+1)*n_theta+j],vertex[(i+1)*n_theta+(j+1)%n_theta]
			r3,r4=vertex[i*n_theta+(j+1)%n_theta],vertex[nn+i*n_theta+j]
			for k in [r0,r4,r1, r1,r4,r2, r2,r4,r3, r3,r4,r0]:
				vertex_.append(k)
	#cap
	vertex_cap_=get_1st_generation_square_pyramid_without_bottom()
	vertex_cap_=on_sphere(vertex_cap_)
	for i in range(cap_discretizetion_level):
		vertex_cap_=create_next_generation(vertex_cap_)
		vertex_cap_=on_sphere(vertex_cap_)
	for i in range(len(vertex_cap_)):
		vertex_cap_[i][0]=radius*vertex_cap_[i][0]
		vertex_cap_[i][1]=radius*vertex_cap_[i][1]
		vertex_cap_[i][2]+=half_length
	vertex_cap2_=zeros((len(vertex_cap_),3),float)
	for i in range(len(vertex_cap2_)):
		vertex_cap2_[i][0]=vertex_cap_[i][1]
		vertex_cap2_[i][1]=vertex_cap_[i][0]
		vertex_cap2_[i][2]=-vertex_cap_[i][2]
	vertex,face=reduce_vertex(array(array(vertex_).tolist()+array(vertex_cap_).tolist()+array(vertex_cap2_).tolist()))
	return vertex,face

def torus(small_radius,large_radius,n_theta,n_circumference):
	vertex=[]
	for i in range(n_circumference+1):
		z=-1.0+2.0*i/n_circumference
		for j in range(n_theta):
			th=2.*pi*j/n_theta
			x,y=cos(th),sin(th)
			xp=(large_radius+small_radius*x)*cos(2*pi*z/2)
			zp=(large_radius+small_radius*x)*sin(2*pi*z/2)
			vertex.append([xp,small_radius*y,zp])
	vertex_=[]
	for i in range(n_circumference):
		for j in range(n_theta):
			vertex_.append(vertex[i*n_theta+j])
			vertex_.append(vertex[i*n_theta+(j+1)%n_theta])
			vertex_.append(vertex[(i+1)*n_theta+j])
			vertex_.append(vertex[(i+1)*n_theta+j])
			vertex_.append(vertex[i*n_theta+(j+1)%n_theta])
			vertex_.append(vertex[(i+1)*n_theta+(j+1)%n_theta])
	vertex,face=reduce_vertex(array(vertex_))
	return vertex,face

def troidal_twist_ribbon(width,large_radius,twist_number,n_width,n_circumference):
	ch=1.0
	if twist_number<0.0:
		ch=-1.0
	#######create twist ribbon
	l,Nx,Nz=1.0,n_width,n_circumference
	###rectangular lattice
	sx,sz=1.0/Nx,1.0/Nz
	r=[0]*((Nx+1)*(Nz+1))
	for j in range(Nz+1):
		for i in range(Nx+1):
			index=j*(Nx+1)+i
			x,z=2*(i*sx-0.5),2*(j*sz-0.5)
			r[index]=[x,0.0,z]
	#right_helicoid
	if twist_number!=0:
		a,b=width,2./abs(twist_number)/pi
		for i in range(len(r)):
			r[i]=[a*r[i][0]*cos(l*r[i][2]/b),ch*a*r[i][0]*sin(l*r[i][2]/b),l*r[i][2]]
	else:
		for i in range(len(r)):
			r[i]=[width*r[i][0],0.0,l*r[i][2]]
	vertex_=[]
	for j in range(Nz):
		for i in range(Nx):
			index=j*(Nx+1)+i
			vertex_.append(r[index])
			vertex_.append(r[index+1])
			vertex_.append(r[index+1+Nx+1])
			vertex_.append(r[index+Nx+1])
			vertex_.append(r[index])
			vertex_.append(r[index+1+Nx+1])
	for i in range(len(vertex_)):
		x,y,z=vertex_[i]
		yp=(large_radius+y)*cos(2*pi*z/2)
		zp=(large_radius+y)*sin(2*pi*z/2)
		vertex_[i]=[x,yp,zp]
	vertex,face=reduce_vertex(array(vertex_))
	return vertex,face

#radius helix_a
#pitch 2 pi helix_b
def helix(helix_a,helix_b,radius,half_length,n_theta,n_height,ch=1):
	vertex,face=cylinder(radius,half_length,n_theta,n_height)
	for i in range(len(vertex)):
		s=vertex[i][2]
		a,b=helix_a,helix_b
		pos=array([a*cos(s/b),ch*a*sin(s/b),s])
		c=sqrt(a*a+b*b)
		#e1=array([-a*sin(s/b),ch*a*cos(s/b),b])/c
		e2=array([-cos(s/b),-ch*sin(s/b),0])
		e3=array([b*ch*sin(s/b),-b*cos(s/b),a*ch])/c
		vertex[i]=pos+vertex[i][0]*e2+vertex[i][1]*e3
	return vertex,face

#radius helix_a
#pitch 2 pi helix_b
def helix_with_caps(helix_a,helix_b,radius,half_length,cap_discretizetion_level,n_height,ch=1):
	vertex,face=cylinder_with_caps(radius,half_length,cap_discretizetion_level,n_height)
	for i in range(len(vertex)):
		s=vertex[i][2]
		a,b=helix_a,helix_b
		pos=array([a*cos(s/b),ch*a*sin(s/b),s])
		c=sqrt(a*a+b*b)
		#e1=array([-a*sin(s/b),ch*a*cos(s/b),b])/c
		e2=array([-cos(s/b),-ch*sin(s/b),0])
		e3=array([b*ch*sin(s/b),-b*cos(s/b),a*ch])/c
		vertex[i]=pos+vertex[i][0]*e2+vertex[i][1]*e3
	return vertex,face

def slightlyDeformedSphere(radius,discretization_level,beta,monopole,dipole,quadrupole,mother="octahedron"):
	vertex,face=sphere(radius,discretization_level,mother)
	normal_vertex=zeros((len(vertex),3),float)
	for i in face:
		r0,r1,r2=array(vertex[i[0]]),array(vertex[i[1]]),array(vertex[i[2]])
		l0,l1,l2=dot(r2-r1,r2-r1),dot(r2-r0,r2-r0),dot(r1-r0,r1-r0)
		w0,w1,w2=l1+l2,l0+l2,l0+l1
		n=array(cross(r1-r0,r2-r0))
		n=n/sqrt(dot(n,n))
		normal_vertex[i[0]]+=n/w0
		normal_vertex[i[1]]+=n/w1
		normal_vertex[i[2]]+=n/w2
	for i in range(len(normal_vertex)):
		normal_vertex[i]/=sqrt(dot(normal_vertex[i],normal_vertex[i]))
	dipole,quadrupole=array(dipole),array(quadrupole)
	for i in range(len(vertex)):
		n=normal_vertex[i]
		f=monopole+3.*dot(dipole,n)+2.5*dot(dot(n,quadrupole),n)
		vertex[i]=vertex[i]*(1.0+beta*f)
	return vertex,face

#import MMesh5
#vertex,face=MMesh5.sphere(1,3,"icosahedron")
#d,q=array([0.0,0.0,2./3]),0.4*array([[-1.,0,0],[0.,-1,0],[0,0,2]])
#vertex,face=MMesh5.slightlyDeformedSphere(1,3,0.2,0,d,q,"icosahedron")
#vertex,face=MMesh5.hemisphere_with_bottom(1,3)
#vertex,face=MMesh5.ellipsoid(4,2,1,3)
#vertex,face=MMesh5.twist_ribbon(1,1,10,3,100,1)
#vertex,face=MMesh5.biconcave(1.38581895,0.207,2.003,-1.123,3)
#vertex,face=MMesh5.hemi_ellipsoid_with_bottom(1,1,3,3)
#vertex,face=MMesh5.hemi_ellipsoid_without_bottom(1,1,3,4)
#vertex,face=MMesh5.rect3d(2,2,2,3)
#vertex,face=MMesh5.disk(2,3)
#vertex,face=MMesh5.rect2d(2,3,3)
#vertex,face=MMesh5.ellipse(2,3,3)
#vertex,face=MMesh5.cylinder(1,10,20,20)
#vertex,face=MMesh5.cylinder_with_caps(1,10,3,20)
#vertex,face=MMesh5.cylinder_with_spherical_caps(1,2.5,3,15)
#vertex,face=MMesh5.torus(1.0,3.0,10,10)
#vertex,face=MMesh5.troidal_twist_ribbon(1.0,3.0,1,5,50)
#vertex,face=MMesh5.helix(3.0,1.0,1.0,15,10,50)
#vertex,face=MMesh5.helix_with_caps(3.0,1.0,1.0,15,2,50)

#MMesh5.draw_face_with_normal(vertex,face,_udf_)
#MMesh5.draw_face(vertex,face,_udf_)
#print(MMesh5.euler_number(vertex,face))
