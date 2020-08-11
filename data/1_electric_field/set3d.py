import MMesh5
from numpy import *
size=1.0
radius=1.0
$input.vertex[]=[]
$input.face[]=[]
$input.region_condition[]=[]
$input.evaluation_point[]=[]

vertex,face=MMesh5.sphere(size,3,"icosahedron")
nv=len(vertex)
for i in range(nv):
	$input.vertex[i].id=i
	$input.vertex[i].position=(radius*vertex[i]).tolist()
nf=len(face)
$input.region_condition[0].name="sphere"
$input.region_condition[0].type="q"
$input.region_condition[0].q=0
for i in range(nf):
	$input.face[i].id=i
	$input.face[i].vertex[]=[face[i][0],face[i][1],face[i][2]]
	$input.region_condition[0].face[i]=i
$input.external_field.phi=0.0
$input.external_field.gradient_phi=[1.0,0.0,0.0]
