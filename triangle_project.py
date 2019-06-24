#!/usr/bin/env python

import vtk
from sys import argv
import numpy as np
from vtk.util import numpy_support
import math
import csv
import numpy
# import matplotlib.pyplot as plt
# import seaborn

# import pickle

RADIUS_TIMES_AVERAGE_EDGE_LENGTH=1

def ClipMRIWithCylinder(mri_data,cylinderTransformed):
    implicitPolyDataDistance = vtk.vtkImplicitPolyDataDistance()
    implicitPolyDataDistance.SetInput(cylinderTransformed)
    mri_data_copy=mri_data
    # Create an array to hold distance information
    signedDistances = vtk.vtkFloatArray()
    signedDistances.SetNumberOfComponents(1)
    signedDistances.SetName("SignedDistances")


    # Evaluate the signed distance function at all of the grid points
    for pointId in range(mri_data.GetNumberOfPoints()):
        p = mri_data.GetPoint(pointId)
        signedDistance = implicitPolyDataDistance.EvaluateFunction(p)
        signedDistances.InsertNextValue(signedDistance)
    
    # add the SignedDistances to the grid
    mriDistanceData=np.zeros([int(mri_data.GetNumberOfPoints()),1])

    for i in range(mri_data.GetNumberOfPoints()):
        mriDist=signedDistances.GetTuple(i)[0]
        mriDistanceData[i]=mriDist
       
        
    mriDistVtk=vtk.vtkDoubleArray()
    mriDistVtk.SetArray(mriDistanceData,mri_data.GetNumberOfPoints(),1)
    mri_data_copy.GetPointData().SetScalars(mriDistVtk)

    clipper = vtk.vtkClipDataSet()
    clipper.SetInputData(mri_data_copy)
    clipper.InsideOutOn()
    clipper.SetValue(0.0)
    clipper.Update()

    return clipper


    
    

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
# sreader.SetFileName("target.vtk")
sreader.Update()
print('Finished reading mesh.')
trianglepd= sreader.GetOutput()

print('Reading lge mesh data...')
sreader1 = vtk.vtkPolyDataReader()
sreader1.SetFileName(argv[2])
# sreader1.SetFileName("advanced_all.vtk")
sreader1.Update()
print('Finished reading mesh.')
mri_data= sreader1.GetOutput()

mri_inten=numpy_support.vtk_to_numpy(mri_data.GetPointData().GetArray('MetaImage'))
max_intn=np.amax(mri_inten)
print(max_intn)
delaunay3D = vtk.vtkDelaunay3D()
delaunay3D.SetInputData(mri_data)
delaunay3D.Update()

geofilter = vtk.vtkGeometryFilter()
geofilter.SetInputData(delaunay3D.GetOutput())
geofilter.Update()
MRI_delaunay3D_pd= geofilter.GetOutput()

# # Read in image file
print('Reading image')
ireader = vtk.vtkMetaImageReader()
ireader.SetFileName(argv[3])
# ireader.SetFileName("pt_002.mha")
ireader.Update()
print('Finished reading image.')


normFilter = vtk.vtkPolyDataNormals()
normFilter.ComputeCellNormalsOn()
normFilter.SetInputData(trianglepd)
normFilter.SplittingOff()  # Prevents us ending up with more points than started with
normFilter.ConsistencyOn()
normFilter.Update()
cellNormalsRetrieved = normFilter.GetOutput().GetCellData().GetNormals()

cellcenter=vtk.vtkCellCenters()
cellcenter.SetInputData(normFilter.GetOutput())
cellcenter.Update()
cell_centeroutput=cellcenter.GetOutput()

cell_centeroutput_num=cell_centeroutput.GetNumberOfPoints()

center=[0.0,0.0,0.0]
centerOfMassFilter =vtk.vtkCenterOfMass()
centerOfMassFilter.SetInputData(trianglepd)
centerOfMassFilter.SetUseScalarsAsWeights(0)
centerOfMassFilter.Update()
xx=centerOfMassFilter.GetCenter()


p=[0.0,0.0,0.0]
epsilon=0.00001
cN=[0.0,0.0,0.0] 
x=[0.0,0.0,0.0] 
pcoords=x=[0.0,0.0,0.0] 
subId = vtk.mutable(0)
inf=100000000
t=vtk.mutable(0)

size_normals = normFilter.GetOutput().GetNumberOfCells()
print(size_normals)
value=[]
value1=[]
normalvec=[]
cell_center=[]
cell_cent_shift=[]
for i in range(0,size_normals):
    cellcenter.GetOutput().GetPoint(i,p)
    cellNormalsRetrieved.GetTuple(i, cN)
    
    vec2=np.subtract(xx,p)

    a = np.dot(vec2,cN)
    cN=np.array(cN)
    if a > 0:
        cN=-1*cN
    normalvec.append((list(cN)))
    cell_center.append((list(p)))

    cell_center_shiftedontocylin=p+10*cN
    cell_cent_shift.append((list(cell_center_shiftedontocylin)))

    
cellId = vtk.mutable(0)
pointIdList=vtk.vtkIdList()
dataList=vtk.vtkIdList()


outputCellIds = []
outputPointIds = []


cN=[0.0, 0.0, 0.0]
p=[0.0, 0.0, 0.0]
cellt=vtk.vtkGenericCell()
new_point=[0.0, 0.0, 0.0]
cellarray=vtk.vtkCellArray()
pts=vtk.vtkPoints()
numberOfCells = trianglepd.GetNumberOfCells()
p=[0.0,0.0,0.0]
averag_ed_length=[]

def averag_length(tri_pts):
    A = np.array(tri_pts[0])
    B = np.array(tri_pts[1])
    C = np.array(tri_pts[2])  
    AB = B - A
    AC = C - A
    BC = C - B
    A_len = np.linalg.norm(AB)
    B_len = np.linalg.norm(AC)
    C_len = np.linalg.norm(BC)
    mean_len=(A_len+B_len+C_len)/3

    return mean_len

for j in range(0, numberOfCells):

    sreader.GetOutput().GetCell(j, cellt)
    pts= cellt.GetPoints() 
    tri_pts=[]
    avg_len=[]
    p=[0.0,0.0,0.0]
    i=0
    
    for i in range(0,3):       
        pts.GetPoint(i,p)
        tri_pts.append((list(p)))

    avg_len=averag_length(tri_pts)  

    averag_ed_length.append(avg_len)

test=[]


cylin_axis=[0, 1, 0]
dotp=np.matmul(normalvec,cylin_axis)

crossprod=np.zeros((size_normals,3))
theta_deg=np.zeros((size_normals,1))

# crossprod=np.cross(normalvec,cylin_axis)
# theta=np.arccos(dotp)
# theta_deg=360-numpy.rad2deg(theta) #subtyract 360 from angle to realign arccos

append = vtk.vtkAppendPolyData()
averageIntensities=np.zeros([size_normals,1])
prop_scar=np.zeros([size_normals,1])
SIZE_points=np.zeros([size_normals,1])

 
# dataListToSave=[]



for ii in range (2525,2526): 

    print(ii)
    crossprod[ii]=np.cross(normalvec[ii],cylin_axis)
    theta=np.arccos(dotp[ii])
    theta_deg[ii]=360-numpy.rad2deg(theta) #subtyract 360 from angle to realign arccos

  
    cylinder = vtk.vtkCylinderSource()
    cylinder.SetRadius(0.4*averag_ed_length[ii])
    cylinder.SetHeight(20)
    cylinder.SetResolution(360)
    cylinder.Update()

    transform=vtk.vtkTransform()
    transform.PostMultiply()
    transform.RotateWXYZ(theta_deg[ii],crossprod[ii])
    transform.Translate(cell_cent_shift[ii])
    transform.Update()
    transformPD1 = vtk.vtkTransformPolyDataFilter()
    transformPD1.SetTransform(transform)
    transformPD1.SetInputData(cylinder.GetOutput())
    transformPD1.Update()


    append.AddInputData(transformPD1.GetOutput())
    append.Update()

    
    # clip=ClipMRIWithCylinder(mri_data,transformPD1.GetOutput())


    # size_points = clip.GetOutput().GetNumberOfPoints()
    # SIZE_points[ii]=size_points
    # print(size_points)
    # if size_points==0:
    #     averageIntensities[ii]=np.nan
    #     allIntensities=np.nan
    #     prop_scar[ii]=np.nan
    #     continue

    # ptsClipped=vtk.vtkPoints()
    # ptsClipped=clip.GetOutput().GetPoints()
    # ptsClippedPD=vtk.vtkPolyData()
    # ptsClippedPD.SetPoints(ptsClipped)
    
  


    # probeFilter=vtk.vtkProbeFilter()
    # probeFilter.SetInputData(ptsClippedPD)
    # probeFilter.SetSourceData(ireader.GetOutput())
    # probeFilter.Update()
    
    # allIntensities=numpy_support.vtk_to_numpy(probeFilter.GetOutput().GetPointData().GetArray('MetaImage')) 

    # # import pdb;pdb.set_trace()
    # # num_true= (np.count_nonzero(allIntensities>(max_intn/2)))
    # # print(num_true)
    # prop_scar[ii]= float((np.count_nonzero(allIntensities>(max_intn/2))))/size_points
    # print(prop_scar[ii])
    
    
    # averageIntensities[ii]=np.mean(allIntensities)
    # print(averageIntensities[ii])
    # # dataListToSave.append(allIntensities)

writer=vtk.vtkPolyDataWriter()
writer.SetInputData(append.GetOutput())
writer.SetFileName('cylinder.vtk')
writer.Write()
import pdb; pdb.set_trace()
# averageIntensities.tofile(argv[4])
# prop_scar.tofile(argv[5])
   

# 
allCVs=np.fromfile(argv[4])

filt=[~np.isinf(allCVs)]
allCVs=allCVs[filt]
# import pdb; pdb.set_trace()
CVs_reshape=allCVs[:size_normals]
CVs_reshape=CVs_reshape.reshape(averageIntensities.shape)
avgints_filt=averageIntensities[~np.isnan(averageIntensities)]
prop_scar_filt=prop_scar[~np.isnan(averageIntensities)]
CVs_reshape_filt=CVs_reshape[~np.isnan(averageIntensities)]

size_pts_np=np.array(SIZE_points)
SIZE_points_filt=size_pts_np[~np.isnan(averageIntensities)]



avgints_filt.tofile(argv[5])
prop_scar_filt.tofile(argv[6])
CVs_reshape_filt.tofile(argv[7])
SIZE_points_filt.tofile(argv[8])

# getCorrelationCVIntensity(allCVs,averageIntensities)
# cvs_filt=allCVs.reshape(averageIntensities.shape)[~np.isnan(averageIntensities)]


# ax = seaborn.regplot(x=avgints_filt, y=allCVs)
# plt.show()


# SIZE_points_filt2=SIZE_points_filt[SIZE_points_filt>50]
# prop_scar_filt2=prop_scar_filt[SIZE_points_filt>50]
# CVs_reshape_filt2=CVs_reshape_filt[SIZE_points_filt>50]
# avgints_filt2=avgints_filt[SIZE_points_filt>50]     
    
  