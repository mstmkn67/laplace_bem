from math import *
$input.vertex[]=[]
$input.edge[]=[]
$input.region_condition[]=[]
$input.evaluation_point[]=[]
radius=1.0
division=5000
#### mesh ###
index=0
dphi=2.*pi/division
for i in range(division):
	phi=i*dphi
	x,y=radius*cos(phi),radius*sin(phi)
	$input.vertex[index]=[index,[x,y]]
	index+=1
for i in range(index-1):
	$input.edge[i]=[i,[i,i+1]]
	$input.region_condition[0].edge[i]=i
$input.edge[index-1]=[index-1,[index-1,0]]
$input.region_condition[0].edge[index-1]=index-1
$input.region_condition[0].type="q"
$input.region_condition[0].q=0.0
$input.external_field.phi=0.0
$input.external_field.gradient_phi=[1.0,0.0]

