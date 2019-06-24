#!/usr/bin/env python

import vtk
from sys import argv
import numpy as np
from vtk.util import numpy_support
import math
import csv
import numpy


def euclidean2(v1,v2):  
    if type(v1)==list:
        v1=np.array(v1)  
    if type(v2)==list:
        v2=np.array(v2)
    diff = v1 - v2
    squareDistance = np.dot(diff.T, diff)
    dist=math.sqrt(squareDistance)
    return dist

#==================================================================
# MAIN 
#==================================================================

# Read in mesh files


print('Reading mesh...')
sreader = vtk.vtkPolyDataReader()
sreader.SetFileName(argv[1])
sreader.Update()
print('Finished reading mesh.')
epicardial_wall= sreader.GetOutput()

print('Reading mesh...')
sreader = vtk.vtkPolyDataReader()
sreader.SetFileName(argv[2])
sreader.Update()
print('Finished reading mesh.')
endocardial_wall= sreader.GetOutput()


# # Read in image file
print('Reading image')
ireader = vtk.vtkMetaImageReader()
ireader.SetFileName(argv[3])
ireader.Update()
print('Finished reading image.')



# Create surface normals filter

normFilter = vtk.vtkPolyDataNormals()
normFilter.ComputePointNormalsOn()
normFilter.SplittingOff() # Prevents us ending up with more points than started with
normFilter.ConsistencyOn()
if vtk.VTK_MAJOR_VERSION <=5:
    normFilter.SetInput(endocardial_wall)
else:
    normFilter.SetInputData(endocardial_wall)
normFilter.Update()

pointNormalsRetrieved = normFilter.GetOutput().GetPointData().GetNormals()





cN=[0.0,0.0,0.0] 
p=[0.0,0.0,0.0]
size_points = normFilter.GetOutput().GetNumberOfPoints() 
pointNormalsArray =vtk.vtkDoubleArray()
pointNormalsArray.SetNumberOfComponents(3)
pointNormalsArray.SetNumberOfTuples(10000000)

pts_Inner=vtk.vtkPoints()
for i in range(0,size_points):
  normFilter.GetOutput().GetPoint(i,p)
  pointNormalsRetrieved.GetTuple(i, cN)
  pts_Inner.InsertNextPoint(p)
  pointNormalsArray.SetTuple(i,cN)


kdTree =vtk.vtkKdTreePointLocator()
kdTree.SetDataSet(epicardial_wall)
kdTree.BuildLocator()

ps=[0.0,0.0,0.0]
pt=[0.0,0.0,0.0]
x=[0.0,0.0,0.0] 
distmat=[]
pts_Intersect =vtk.vtkPoints()
size_pts=endocardial_wall.GetNumberOfPoints()
for i in range(0,size_pts):
        
        pts_Inner.GetPoint(i,ps)
       
        iD= kdTree.FindClosestPoint(ps)
     
        kdTree.GetDataSet().GetPoint(iD,x)
        dist = euclidean2(x,ps)
        distmat.append(dist)

        pts_Intersect.InsertNextPoint(x)


for max_wall_thickness in np.arange (2.5, 3.5, 2.5):
    print(max_wall_thickness)
   
    # Setting maximum wall thickness
    distmat_filt=np.array(distmat)
    distmat_filt[distmat_filt>max_wall_thickness]=max_wall_thickness
    distmat_filt=distmat_filt.tolist()

  
    endocardial_wall.GetPointData().SetNormals(pointNormalsArray)
    pointNormalsRetrievedendo = endocardial_wall.GetPointData().GetNormals()


    f = open("dist_%s.txt" %max_wall_thickness, "w")
    f.write("\n".join(map(lambda x: str(x), distmat_filt)))
    f.close()



    list_of_nodes=[]
    list_of_normals=[]
    size_pts=endocardial_wall.GetNumberOfPoints()
    for k in range(0, size_pts):
        point = [0.0,0.0,0.0]
        pN=[0.0,0.0,0.0] 
        endocardial_wall.GetPoint(k, point)
        pointNormalsRetrievedendo.GetTuple(k,pN)
        list_of_nodes.append(point)
        list_of_normals.append(pN)
        
 


    with open("dist_%s.txt" %max_wall_thickness) as f:
        list_of_advances=f.readlines()
        
        list_of_advances=[float(i.strip("/n")) for i in list_of_advances]

    assert len(list_of_nodes)==len(list_of_advances)
# print("list of advances:",  list_of_advances[:10])


    stepsize=1.7
    list_of_advances=[a//stepsize for a in list_of_advances]

# print("list of advances:",  list_of_advances[:50])
    nodes_advances = zip(list_of_nodes, list_of_advances,list_of_normals)


    min_advance = int(min(list_of_advances))
    max_advance = int(max(list_of_advances))
    print(max_advance)
    min_advance=0
   
    current_advance=min_advance



    idstoremove=[]


    iterations = 1
    advanced_points_final=[]



    while current_advance<max_advance:

  # nodestoexclude= filter(lambda i: i[1] < current_advance, nodes_advances)
        new_nodes_advances = []
        for node, advance, normal in nodes_advances:
            if advance >= current_advance:
             new_nodes_advances.append((node, advance, normal))
  # print(len(new_nodes_advances))
  # print(len(nodes_advances))
  
 # assert len(new_nodes_advances)<=   len(nodes_advances)     
  # take new nodes advances and get new mesh from it
        new_points = [n[0] for n in new_nodes_advances]
        new_normals= [n[2] for n in new_nodes_advances]

        new_points = np.array(new_points)
        new_normals = np.array(new_normals)

  # print(new_points.shape, new_normals.shape)
  # print("new points:", new_points[:10])
  # print("new normals:", new_normals[:10])  
  # print("iteration:", iterations)
  # print("length:", stepsize*iterations)

  # print("advancement",((stepsize*iterations)*new_normals))
  # advanced_points=new_points + ((stepsize*iterations)*new_normals)
        advanced_points=new_points + ((stepsize*iterations)*new_normals)
        advanced_points_final.append((advanced_points))

    #   advancedpoints_polyVTK=vtk.vtkPolyData()
    #   advancedpoints_VTK=vtk.vtkPoints()

    #   for p in list(advanced_points):
    #     advancedpoints_VTK.InsertNextPoint(list(p))

    #   advancedpoints_polyVTK.SetPoints(advancedpoints_VTK)
    #   writer=vtk.vtkPolyDataWriter()
    #   writer.SetInputData(advancedpoints_polyVTK)
    #   writer.SetFileName('advancedvtk_%s.vtk' % str(current_advance))
    #   writer.Write()

        nodes_advances = new_nodes_advances
        current_advance+=1
    # print("current_advance",current_advance)
        iterations += 1
    # print("iteration",iterations)

        advancedpointsf_polyVTK=vtk.vtkPolyData()
        advancedpointsf_VTK=vtk.vtkPoints()
        # print("advancepoints_final",advanced_points_final[:10])
        adv_points=[]
        new_list=[l.tolist() for l in advanced_points_final]
        i=1
        for list in new_list:
            for number in list:
                adv_points.append(number)
                advancedpointsf_VTK.InsertNextPoint((number))

        advancedpointsf_polyVTK.SetPoints(advancedpointsf_VTK)
        probeArray =vtk.vtkDoubleArray()
        probeArray.SetNumberOfComponents(1)
        probeArray.SetNumberOfTuples(10000000)




        probeFilter= vtk.vtkProbeFilter()
        probeFilter.SetInputData(advancedpointsf_polyVTK)
        probeFilter.SetSourceData(ireader.GetOutput())
        probeFilter.Update()

        probeArray=probeFilter.GetOutput().GetPointData().GetArray('MetaImage')
        advancedpointsf_polyVTK.GetPointData().SetScalars(probeArray)


        writer=vtk.vtkPolyDataWriter()
        writer.SetInputData(advancedpointsf_polyVTK)
        writer.SetFileName('advancedvtk_endo_%s.vtk' %max_wall_thickness )
        writer.Write()

        probeArray2 =vtk.vtkDoubleArray()
        probeArray2.SetNumberOfComponents(1)
        probeArray2.SetNumberOfTuples(10000000)


