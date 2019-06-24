#!/usr/bin/env python
from __future__ import print_function

import vtk
from sys import argv
import math
import collections
from collections import defaultdict
import csv
import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt


MIN_LEN=0.5     
MAX_LEN=100
SHRINK_FACTOR=10
ASPECT_RATIO_FACTOR=100000
MAX_SPEED=100


def area_of_triangle(a, b, c):
    ab = b - a
    ac = c - a
    bc = c - b

    a_len = np.linalg.norm(ab)
    b_len = np.linalg.norm(ac)
    c_len = np.linalg.norm(bc)

    s = (a_len + b_len + c_len) / 2
    return math.sqrt(s * (s - a_len) * (s - b_len) * (s - c_len))

print('Reading mesh...')
sreader = vtk.vtkPolyDataReader()
sreader.SetFileName(argv[1])
#sreader.SetFileName("pt_053.vtk")
sreader.Update()
print('Finished reading mesh.')
test = sreader.GetOutput()

geofilter = vtk.vtkGeometryFilter()
geofilter.SetInputData(test)
geofilter.Update()
mesh= geofilter.GetOutput()

featureedges=vtk.vtkFeatureEdges()
featureedges.SetInputData(mesh)
featureedges.BoundaryEdgesOn()
featureedges.ManifoldEdgesOn()
featureedges.NonManifoldEdgesOff()
featureedges.FeatureEdgesOff()
featureedges.Update()

connectivity=vtk.vtkConnectivityFilter()
connectivity.SetInputData(featureedges.GetOutput())
connectivity.SetExtractionModeToAllRegions()
connectivity.Update()

triangles=vtk.vtkStripper()
triangles.SetInputData(connectivity.GetOutput())
triangles.Update()

triangles2=vtk.vtkTriangleFilter()
triangles2.SetInputData(triangles.GetOutput())
triangles2.Update()

writer=vtk.vtkPolyDataWriter()
writer.SetInputData(triangles2.GetOutput())
writer.SetFileName('triangless.vtk')
writer.Write()

original = sreader.GetOutput()
pointScalarsRetrieved = mesh.GetPointData().GetArray('LAT')
list_of_scalars=[]
size_pts=mesh.GetNumberOfPoints()
""" for k in range(0, size_pts):
    point = [0.0, 0.0, 0.0]
    pS=[0.0]
    mesh.GetPoint(k, point)
    pointScalarsRetrieved.GetTuple(k, pS)
    list_of_scalars.append(pS[0]) """
    
vtkcells = vtk.vtkCellData()
vtkcells=mesh.GetCellData()


def GetGeomCenterOfMass(mesh):


    cellcenter=vtk.vtkCellCenters()
    cellcenter.SetInputData(mesh)
    cellcenter.Update()
    cell_centeroutput=cellcenter.GetOutput()

    cell_centeroutput_num=cell_centeroutput.GetNumberOfPoints()

    center=[0.0,0.0,0.0]
    centerOfMassFilter =vtk.vtkCenterOfMass()
    centerOfMassFilter.SetInputData(mesh)
    centerOfMassFilter.SetUseScalarsAsWeights(0)
    centerOfMassFilter.Update()
    xx=centerOfMassFilter.GetCenter()
    return xx,cell_centeroutput

def IsNormFlipped(COM_Total,pt, cell_normal):

    
    vec2=np.subtract(COM_Total,pt)

    a = np.dot(vec2,cell_normal)
    

    #normal is inwards if the dot produt is more than 0
    print(a)
    return a>0


def is_inside_mesh(triangles2, point):
    selectEnclosedPoints=vtk.vtkSelectEnclosedPoints()  
    selectEnclosedPoints.SetSurfaceData(triangles2)
    selectEnclosedPoints.CheckSurfaceOn()
    selectEnclosedPoints.SetTolerance(0.001)

    testpt=vtk.vtkPoints()
    pointpoly=vtk.vtkPolyData()
    testpt.InsertNextPoint(point)
    pointpoly.SetPoints(testpt)

    selectEnclosedPoints.SetInputData(pointpoly)
    selectEnclosedPoints.Update()
    _is_inside = selectEnclosedPoints.IsInside(0)
    selectEnclosedPoints.Update()
    return _is_inside

def move_point(point, normal, shrink_factor):
    return np.array(point)-shrink_factor*np.array(normal)

def createShrunkVolume(mesh):
        
    pts_shrink=vtk.vtkPoints()
    smoother = vtk.vtkSmoothPolyDataFilter()
    smoother.SetInputData(mesh)
    smoother.SetNumberOfIterations(300)
    smoother.Update()

    numberOfCells = mesh.GetNumberOfCells()

    normFilterForPoints = vtk.vtkPolyDataNormals()
    normFilterForPoints.ComputePointNormalsOn()
    normFilterForPoints.SetInputData(smoother.GetOutput())
    normFilterForPoints.SplittingOff()  # Prevents us ending up with more points than started with
    normFilterForPoints.ConsistencyOff()
    normFilterForPoints.AutoOrientNormalsOff()
    normFilterForPoints.NonManifoldTraversalOn()
    normFilterForPoints.Update()
    normFilterAfter=normFilterForPoints.GetOutput()
    pointNormalsRetrieved = normFilterForPoints.GetOutput().GetPointData().GetNormals()

    normFilterForCells = vtk.vtkPolyDataNormals()
    normFilterForCells.ComputeCellNormalsOn()
    normFilterForCells.SetInputData(mesh)
    normFilterForCells.SplittingOff()  # Prevents us ending up with more points than started with
    normFilterForCells.ConsistencyOn()
    normFilterForCells.Update()
    cellNormalsRetrieved = normFilterForCells.GetOutput().GetCellData().GetNormals()

    COM_Total,COM_each_cell=GetGeomCenterOfMass(mesh)

    cN=[0.0, 0.0, 0.0]
    p=[0.0, 0.0, 0.0]
    cellt=vtk.vtkGenericCell()
    new_point=[0.0, 0.0, 0.0]
    size_points = normFilterForPoints.GetOutput().GetNumberOfPoints()

    size_cells = normFilterForCells.GetOutput().GetNumberOfCells()

    for i in range(0, size_points):
    
        normFilterForPoints.GetOutput().GetPoint(i, p)
        pointNormalsRetrieved.GetTuple(i, cN)
        current_shrink_factor = SHRINK_FACTOR
        inside_mesh = False

        while not inside_mesh and current_shrink_factor >=  0.1 * SHRINK_FACTOR:
            # if i == 4266:
            #     print('something')
            new_point = move_point(p, cN, current_shrink_factor)
            inside_mesh = bool(is_inside_mesh(triangles2.GetOutput(), list(new_point)))
            # print("**", current_shrink_factor, inside_mesh, new_point)
            current_shrink_factor = float(current_shrink_factor) * 0.9
  
        if not inside_mesh:
            # print('could not place point inside the mesh with shrink factor %s' % current_shrink_factor * 2)
            pts_shrink.InsertNextPoint(p)
        else:
            pts_shrink.InsertNextPoint(new_point)

    cellarray=vtk.vtkCellArray()

    for j in range(0, numberOfCells):
        sreader.GetOutput().GetCell(j, cellt)
        cellarray.InsertNextCell(cellt)

    # import pdb; pdb.set_trace()
    pd=vtk.vtkUnstructuredGrid()
    pd.SetCells(vtk.VTK_TRIANGLE, cellarray)
    pd.SetPoints(pts_shrink)

    writer=vtk.vtkUnstructuredGridWriter()
    writer.SetInputData(pd)
    writer.SetFileName('pd.vtk')
    writer.Write()

    geofilter = vtk.vtkGeometryFilter()
    geofilter.SetInputData(pd)
    geofilter.Update()
    pd_2= geofilter.GetOutput() 

    ShrunkCOM_Total,a=GetGeomCenterOfMass(pd_2)
    COM_diff=np.asarray(COM_Total)-np.asarray(ShrunkCOM_Total)
    print(COM_diff)
    transform2=vtk.vtkTransform()
    transform2.Translate(COM_diff)
    transform2.Update()

    transformPD1 = vtk.vtkTransformPolyDataFilter()
    transformPD1.SetTransform(transform2)
    transformPD1.SetInputData(pd_2)
    transformPD1.Update()
    pd_3=transformPD1.GetOutput()

    
    writer=vtk.vtkPolyDataWriter()
    writer.SetInputData(pd_3)
    writer.SetFileName('pd_3.vtk')
    writer.Write()  
    return pd, pd_2, pd_3, pts_shrink, size_points

pd,pd_2,pd_3,pts_shrink, size_points=createShrunkVolume(mesh)


pointsPolydata =vtk.vtkPolyData()
pointsPolydata.SetPoints(pts_shrink)
selectEnclosedPoints = vtk.vtkSelectEnclosedPoints()
selectEnclosedPoints.SetInputData(pointsPolydata)
selectEnclosedPoints.SetSurfaceData(pd_2)
selectEnclosedPoints.Update()

outputisinside=[]

pt_isinside = vtk.vtkPoints()
list_of_points2=[]
inside_or_outside=[]
for kk in range (0,size_points):

    xx=selectEnclosedPoints.IsInside(kk)
    inside_or_outside.append(xx)


list_of_points=[]
size_pts=original.GetNumberOfPoints()

# for k in range(0, size_pts):
#     point = [0.0, 0.0, 0.0]
#     original.GetPoint(k, point)
#     list_of_points.append(point)
columns = defaultdict(list)

with open(argv[2]) as f:
    reader = csv.reader(f,delimiter="\t")
    next(reader)
    for row in reader:
        for (i,v) in enumerate(row):
         
            columns[i].append(v)

print('Read in car file')

with open(argv[3]) as f:
    triang_matlab = [list(map(float,rec)) for rec in csv.reader(f, delimiter=',')]
       
triang_mat_int = [tuple(int(x) for x in tup) for tup in triang_matlab]        


A=columns[4]
B=columns[5]
C=columns[6]
E=columns[12]
unipolar=columns[10]
bipolar=columns[11]

A=[float(i) for i in A]
B=[float(i) for i in B]
C=[float(i) for i in C]
list_of_scalars=[float(i) for i in E]
list_of_unipolar=[float(i) for i in unipolar]
list_of_bipolar=[float(i) for i in bipolar]
D=A,B,C
D=np.transpose(D)









# mlab.triangular_mesh([D[:,0]], [D[:,1]], [D[:,2]], indices)
# mlab.show()





car_points_vtk=vtk.vtkPoints()
for i in range(0, D.shape[0]):
    car_points_vtk.InsertNextPoint(D[i])

car_pointsPD=vtk.vtkPolyData()
car_pointsPD.SetPoints(car_points_vtk)

transform2=vtk.vtkTransform()
# transform2.SetMatrix([0.30, -0.49, -0.82, 401.09, -0.60, -0.75, 0.23, 23.97, -0.73, 0.43, -0.53, 312.95, 0,0,0,1]) #pt102
# transform2.SetMatrix([0.47, -0.43, -0.76, 389.67, -0.70, -0.71, -0.02, 59.48, -0.53, 0.55, -0.64, 320.31,  0,0,0,1])#pt099 not working
# transform2.SetMatrix([0.26, -0.3, -0.92, 412.67, -0.45, -0.88, 0.17, -98.34, -0.85, 0.38, -0.37, 387.06, 0,0,0,1])#pt089 
# transform2.SetMatrix([0.71, 0.01, -0.68, 257.57, -0.02, -1.0, -0.06, 25.33, -0.69, 0.07, -0.72, 336.48, 0,0,0,1]) #pt073
# transform2.SetMatrix([0.83, -0.12, -0.55, 194.91, -0.22, -0.96, -0.13, 0.68, -0.51, 0.23, -0.83, 451.68, 0,0,0,1]) #pt068
# transform2.SetMatrix([0.31, -0.20, -0.93, 383.71, -0.54, -0.84, 0, -6.88, -0.78, 0.50, -0.37, 188.47, 0,0,0,1]) #pt060
# transform2.SetMatrix([0.65, -0.09, -0.75, 309.38, -0.40, -0.90, -0s.19, -58.44, -0.66, 0.40, -0.62, 326.25, 0,0,0, 1])#pt021
# transform2.SetMatrix([0.73, 0.68, 0.0, 147.74, -0.39, -0.42, -0.82, 72.83, -0.56, 0.6, -0.57, 49.22,0 ,0,0,1])#pt053
# transform2.SetMatrix([0.39, 0.05, -0.93, 364.26, -0.33, -0.93, -0.19, 42.65,-0.86, 0.37, -0.35, 286.46, 0,0,0,1 ]) #pt002
# transform2.SetMatrix([-0.57, -0.03, -0.80, 342.17, -0.13, -0.99, -0.06, -17.7, -0.81, 0.14, -0.6, 323.16, 0,0,0,1]) #pt066
# transform2.SetMatrix([0.34, -0.26, -0.91, 321.02, -0.43, -0.91, 0.05, -21.94, -0.83, 0.38, -0.42, 330.29, 0, 0, 0, 1])  #pt024
transform2.SetMatrix([0.37, -0.37, -0.86, 276.97, -0.47, -0.87, 0.16, -8.60, -0.81, 0.34, -0.49, 260.60, 0,0,0,1]) #pt046
# transform2.SetMatrix([0.49, -0.31, -0.81, 333.10, -0.51, -0.86, 0.02, -64.51, -0.71, 0.41, -0.58, 409.35, 0,0,0,1]) #pt003

transform2.Update()



transformPD1 = vtk.vtkTransformPolyDataFilter()
transformPD1.SetTransform(transform2)

transformPD1.SetInputData(car_pointsPD)
transformPD1.Update()
car_pointsPD_translate=transformPD1.GetOutput()

writer=vtk.vtkPolyDataWriter()
writer.SetInputData(car_pointsPD_translate)
writer.SetFileName('car_points.vtk')
writer.Write()



list_of_points=[]
size_pts_car=car_pointsPD_translate.GetNumberOfPoints()


for k in range(0,size_pts_car):
    point=[0.0,0.0, 0.0]
    car_pointsPD_translate.GetPoint(k,point)
    list_of_points.append(point)

list_of_points=np.array(list_of_points)

hull=ConvexHull(list_of_points)
indices=hull.simplices

indices=np.array(indices)
indices=np.int64(indices)



tree=vtk.vtkOBBTree()
tree.SetDataSet(pd_2)
tree.BuildLocator()
intersectPoints=vtk.vtkPoints()

def is_inside(idx,inside_or_outside):
    return not bool(inside_or_outside[idx])

def pt_lat_value(idx,list_of_scalars):
    return not bool(list_of_scalars[idx])



def to_points(triangle,idx_to_point):
    return list(map(lambda x: idx_to_point[x], triangle))


def to_scaler(triangle,list_of_scalars):
    return list(map(lambda x: list_of_scalars[x], triangle))

def FindValidPoints(list_of_points):
    valid_points = collections.defaultdict(set)

    idx_to_point = {
        i: p
        for i, p in
        enumerate(list_of_points)
    }

    for i1, p1 in idx_to_point.items():
        for i2, p2 in idx_to_point.items():
        
            if is_inside(i1,inside_or_outside) or is_inside(i2,inside_or_outside):
                continue
            
            if pt_lat_value(i1,list_of_scalars) or pt_lat_value(i2,list_of_scalars):
                continue

            if np.linalg.norm(p1 - p2) < MIN_LEN:
                continue

            if np.linalg.norm(p1 - p2) > MAX_LEN:
                continue
        
            tree.IntersectWithLine(p1, p2, intersectPoints, None)
            # if i1 in (24, 42) and i2 in (24, 42):
            #     
            if not intersectPoints.GetNumberOfPoints():
                # print('test')
                # print(p1)
                # print(p2)
                valid_points[i1].add(i2)

        print("finding intersecting points: %s/%s done" % (i1, len(idx_to_point)), end='\r')

    print("\n")


    sum_of_lens = 0
    max_len = 0
    min_len = 0

    for key, value in valid_points.items():
        value_len = len(value)
        sum_of_lens += value_len
        if value_len > max_len:
            max_len = value_len
        if value_len < min_len:
            min_len = value_len

    print(sum_of_lens / len(valid_points), max_len, min_len)

    return valid_points, idx_to_point

valid_points, idx_to_point = FindValidPoints(list_of_points)

def generate_triangles(points):
    visited_ids = set()  # remember the nodes that we have tested already
    for point_a_id in points:
        for point_b_id in valid_points[point_a_id]:
            if point_b_id == point_a_id:
                raise ValueError  # nodes shouldn't point to themselves
            if point_b_id in visited_ids:
                continue  # we should have already found b->a->??->b
            for point_c_id in points[point_b_id]:
                if point_c_id in visited_ids:
                    continue  # we should have already found c->a->b->c
                if point_a_id in points[point_c_id]:
                    yield(point_a_id, point_b_id, point_c_id)
        visited_ids.add(point_a_id)  # don't search a - we already have all those cycles

print("generating triangles...")

# all_triangles_ids_list = [
#     triangle
#     for triangle in
#     generate_triangles(valid_points)
# ]

# all_triangles_ids_list =[tuple(l) for l in indices]

all_triangles_ids_list=triang_mat_int

print(len(all_triangles_ids_list))
all_triangles_ids_list_np=np.array(all_triangles_ids_list)
all_triangles_ids_list_np.sort(axis=1)
all_triangles_ids_list_np_unique=np.unique(all_triangles_ids_list_np,axis=0)
all_triangles_ids_list_np_unique_list=np.array(all_triangles_ids_list_np_unique).tolist()
all_triangles_ids_list=[tuple(l) for l in all_triangles_ids_list_np_unique_list]
print(len(all_triangles_ids_list))





all_triangle_points = np.array([
    to_points(triangle,idx_to_point)
    for triangle in all_triangles_ids_list
])


all_triangle_scalers = np.array([
    to_scaler(triangle,list_of_scalars)
    for triangle in all_triangles_ids_list
])

all_triangle_unipolar = np.array([
    to_scaler(triangle,list_of_unipolar)
    for triangle in all_triangles_ids_list
])

all_triangle_bipolar = np.array([
    to_scaler(triangle,list_of_bipolar)
    for triangle in all_triangles_ids_list
])

all_triangles_ids = np.array(all_triangles_ids_list)

del all_triangles_ids_list  # del to free memory

print("generated %s triangles" % len(all_triangles_ids))



def get_aspect_ratio_filter(triangles, cutoff=ASPECT_RATIO_FACTOR):
    A = triangles[:, 0, :]
    B = triangles[:, 1, :]
    C = triangles[:, 2, :]

    AB = B - A
    AC = C - A
    BC = C - B

    AB_x_AC = np.cross(AB, AC)

    NOM = (np.cross(AB_x_AC, AB).T * (np.linalg.norm(AC, axis=1) ** 2)).T + (np.cross(AC, AB_x_AC).T * (np.linalg.norm(AB, axis=1) ** 2)).T
    DOM = (2.0 * (np.linalg.norm(AB_x_AC, axis=1) ** 2))
    M = A + (NOM.T / DOM).T
    R = np.linalg.norm(A - M, axis=1)
    AREA_CIRCUMCIRCLE = np.pi * (R**2)

    A_len = np.linalg.norm(AB, axis=1)
    B_len = np.linalg.norm(AC, axis=1)
    C_len = np.linalg.norm(BC, axis=1)
    S = (A_len + B_len + C_len) / 2

    AREA_TRI = (S * (S - A_len) * (S - B_len) * (S - C_len)) ** 0.5

    ASPECT_FILTER = (AREA_CIRCUMCIRCLE / AREA_TRI) > cutoff

    #import pdb; pdb.set_trace()

    return ASPECT_FILTER, AREA_TRI

print("filtering out triangles with bad aspect ratio...")


def get_scaler_filter(scalers, triangles):
    scaler_sort = scalers.argsort(axis=1)

    sorted_scalers = scalers[np.arange(scalers.shape[0])[:, None], scaler_sort]
    sorted_triangles = triangles[np.arange(triangles.shape[0])[:, None], scaler_sort]

    LAT_O = sorted_scalers[:, 0]
    LAT_A = sorted_scalers[:, 1]
    LAT_B = sorted_scalers[:, 2]

    t_OA = LAT_A - LAT_O
    t_OB = LAT_B - LAT_O
    t_AB = LAT_B - LAT_A

    print(t_AB,t_OB,t_OA)
    
    #filter out LATs that are too small
    min_time=MAX_LEN/MAX_SPEED

    #filter1 = np.logical_or(abs(t_AB) <min_time, abs(t_OB) <min_time , abs(t_OA) <min_time)
    filter1=[]
    

    t_OA[filter1]=np.nan
    t_OB[filter1]=np.nan
    t_AB[filter1]=np.nan

    O = sorted_triangles[:, 0]
    A = sorted_triangles[:, 1]
    B = sorted_triangles[:, 2]

    OA = A - O
    OB = B - O
    AB = B - A



    OA_len = np.linalg.norm(OA, axis=1)
    OB_len = np.linalg.norm(OB, axis=1)
    AB_len = np.linalg.norm(AB, axis=1)

    theta = np.arccos(
        ((OA_len**2 + OB_len**2) - AB_len**2) /
        (2* OA_len * OB_len)
    )

    alpha = np.arctan(
        (((t_OB * OA_len) / (t_OA * OB_len)) - np.cos(theta)) /
        np.sin(theta)
    )

    V = (OA_len / t_OA) * np.cos(alpha)

    
    #Finding v vector
    print('finding x')
    X= np.arcsin(np.sin(theta) * OB_len/AB_len)
 
    R= np.pi-alpha-X
    gamma=OA_len* (np.sin(alpha)/np.sin(R))
    print('gamma {}'.format(gamma))
    J_vec= A + ((B-A)/AB_len[:,None])*gamma[:,None]
    print('J_vec')
    print(J_vec)
    #import pdb; pdb.set_trace()
    J_norm_factor=V[:,None]/np.linalg.norm(J_vec-O,axis=1)[:,None]
    V_vec=(J_vec-O)*J_norm_factor
    
    
    print('V_vec')

  

    print('V Shape', V.shape)

    return V, filter1, V_vec


V_s, V_FILTER, V_VEC = get_scaler_filter(all_triangle_scalers, all_triangle_points)  #V_s: matrix of velocities which has been filtered,

# print(V_s.shape)
# print(V_FILTER.shape)
ASPECT_FILTER, AREA_TRI= get_aspect_ratio_filter(all_triangle_points)



def check_if_filt_tri_inside_PD(tripts, mesh):

    points1 =vtk.vtkPoints()
    points1.InsertNextPoint (triangle_points[0])
    points1.InsertNextPoint (triangle_points[1])
    points1.InsertNextPoint (triangle_points[2])
    
    

    triangle =vtk.vtkTriangle()
    triangle.GetPointIds().SetId ( 0, 0 )
    triangle.GetPointIds().SetId ( 1, 1 )
    triangle.GetPointIds().SetId ( 2, 2 )
    
    triangles = vtk.vtkCellArray()
    triangles.InsertNextCell(triangle)
    
    trianglePolyData =vtk.vtkPolyData()
    trianglePolyData.SetPoints( points1 )
    trianglePolyData.SetPolys( triangles )

    A=triangle_points[0]
    B=triangle_points[1]
    C=triangle_points[2]

    a_minus_b=A-B
    a_minus_C=A-C
    normtoplane=np.cross(a_minus_b,a_minus_C)

    centerOfMassFilter =vtk.vtkCenterOfMass()
    centerOfMassFilter.SetInputData(trianglePolyData)
    centerOfMassFilter.Update()
    center=centerOfMassFilter.GetCenter()


    plane = vtk.vtkPlane()
    plane.SetOrigin(center)
    plane.SetNormal(normtoplane)
    planeCut = vtk.vtkCutter()
    planeCut.SetInputData(mesh)
    planeCut.SetCutFunction(plane)
    planeCut.Update()


    num_inter=  planeCut.GetOutput().GetNumberOfPoints()
    
    return num_inter


len_of_ids=all_triangles_ids.shape[0]
intersection_num=[]

for k in range(0, len_of_ids):
   
    pt_ids=all_triangles_ids[k]
    triangle_points=list_of_points[pt_ids]
    
    num_inter= check_if_filt_tri_inside_PD(triangle_points, pd_2)

    intersection_num.append(num_inter)
    print("finding intersecting points: %s/%s done" % (k, len_of_ids), end='\r')


intersection_num=np.array(intersection_num)
PLANE_INTERSECTION_MASK=(intersection_num) > 0

V_FILTER=np.isinf(V_s)



# all_filters= V_FILTER | PLANE_INTERSECTION_MASK | ASPECT_FILTER

all_filters= V_FILTER | ASPECT_FILTER

filter_omitting=~all_filters


V_s_filtered= V_s[filter_omitting]


filterd_tri = all_triangle_points[filter_omitting]
all_triangles_ids_filtered = all_triangles_ids[filter_omitting]

V_vec_filtered=V_VEC[filter_omitting]
AREA_TRI_filtered=AREA_TRI[filter_omitting]


list_of_bipolar_filtered=all_triangle_bipolar[filter_omitting]
list_of_unipolar_filtered=all_triangle_unipolar[filter_omitting]


print("finished filtering")

print("printing a random sample:")
print(all_triangles_ids_filtered)


print('writing triangles to file')

all_triangles_ids_filtered.tofile(argv[4])


list_of_points.tofile(argv[5])

V_s_filtered.tofile(argv[6])

AREA_TRI_filtered.tofile(argv[7])
V_vec_filtered.tofile(argv[8])
list_of_unipolar_filtered.tofile(argv[9])
list_of_bipolar_filtered.tofile(argv[10])

list_of_scalars=np.array(list_of_scalars)
list_of_scalars.tofile(argv[11])


np.savetxt('all_tris.csv', all_triangles_ids_filtered)
np.savetxt('pts.csv', list_of_points)
np.savetxt('vs_filt.csv', V_s_filtered)
np.savetxt('area_filt.csv', AREA_TRI_filtered)
np.savetxt('vvec_filt.csv', V_vec_filtered)
np.savetxt('unipolar.csv', list_of_unipolar_filtered)
np.savetxt('bipolar.csv', list_of_bipolar_filtered)
np.savetxt('lats.csv',list_of_scalars)






