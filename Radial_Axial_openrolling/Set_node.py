# -*- coding: mbcs -*-
# Do not delete the following import lines
session.journalOptions.setValues(replayGeometry=INDEX,recoverGeometry=INDEX)
from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
s='Model-1'
points=24

#上端外圈点SA
point_one=1
point_two=559
diff=558
for i in range(1,points+1):
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)#在视图中显示网格

    a = mdb.models[s].rootAssembly
    n1 = a.instances['Ring-1'].nodes
    nodes1 = n1[(point_one-1+diff*(i-1)):(point_one+diff*(i-1))]
    a.Set(nodes=nodes1, name='SA%s'%i)

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
    regionDef=mdb.models[s].rootAssembly.sets['SA%s'%i]

    mdb.models[s].HistoryOutputRequest(name='XA%s'%i,
    createStepName='Step-1', variables=('COOR1', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)
    regionDef=mdb.models[s].rootAssembly.sets['SA%s'%i]
    
    mdb.models[s].HistoryOutputRequest(name='YA%s'%i,
    createStepName='Step-1', variables=('COOR3', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)

#上圈内点SB
point_one_in = 6
point_two_out = point_one_in +diff
for i in range(1,points+1):
    a = mdb.models[s].rootAssembly
    n1 = a.instances['Ring-1'].nodes

    nodes1 = n1[(point_one_in-1+diff*(i-1)):(point_one_in+diff*(i-1))]#选择节点范围
    a.Set(nodes=nodes1, name='SB%s'%i)

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(adaptiveMeshConstraints=ON)

    regionDef=mdb.models[s].rootAssembly.sets['SB%s'%i]
    #创建x方向的历程输出请求
    mdb.models[s].HistoryOutputRequest(name='XB%s'%i,
    createStepName='Step-1', variables=('COOR1', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)
    #获取刚刚创立的集合
    regionDef=mdb.models[s].rootAssembly.sets['SB%s'%i]
    #创建y方向的历程输出请求
    mdb.models[s].HistoryOutputRequest(name='YB%s'%i,
    createStepName='Step-1', variables=('COOR3', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)

#SC
point_three_in =61
point_four_out =point_three_in + diff 
for i in range(1,points+1):
    a = mdb.models[s].rootAssembly
    n1 = a.instances['Ring-1'].nodes

    nodes1 = n1[(point_three_in-1+diff*(i-1)):(point_three_in+diff*(i-1))]#选择节点范围
    a.Set(nodes=nodes1, name='SC%s'%i)

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(adaptiveMeshConstraints=ON)

    regionDef=mdb.models[s].rootAssembly.sets['SC%s'%i]
    #创建x方向的历程输出请求
    mdb.models[s].HistoryOutputRequest(name='XC%s'%i,
    createStepName='Step-1', variables=('COOR1', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)
    #获取刚刚创立的集合
    regionDef=mdb.models[s].rootAssembly.sets['SC%s'%i]
    #创建y方向的历程输出请求
    mdb.models[s].HistoryOutputRequest(name='YC%s'%i,
    createStepName='Step-1', variables=('COOR3', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)

 #SD
point_five_in = 66
point_six_out = point_five_in + diff
for i in range(1,points+1):
    a = mdb.models[s].rootAssembly
    n1 = a.instances['Ring-1'].nodes

    nodes1 = n1[(point_five_in-1+diff*(i-1)):(point_five_in+diff*(i-1))]#选择节点范围
    a.Set(nodes=nodes1, name='SD%s'%i)

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(adaptiveMeshConstraints=ON)

    regionDef=mdb.models[s].rootAssembly.sets['SD%s'%i]
    #创建x方向的历程输出请求
    mdb.models[s].HistoryOutputRequest(name='XD%s'%i,
    createStepName='Step-1', variables=('COOR1', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)
    #获取刚刚创立的集合
    regionDef=mdb.models[s].rootAssembly.sets['SD%s'%i]
    #创建y方向的历程输出请求
    mdb.models[s].HistoryOutputRequest(name='YD%s'%i,
    createStepName='Step-1', variables=('COOR3', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)
  #SE  
point_7 = 121
point_8 = point_7 + diff
for i in range(1,points+1):
    a = mdb.models[s].rootAssembly
    n1 = a.instances['Ring-1'].nodes

    nodes1 = n1[(point_7-1+diff*(i-1)):(point_7+diff*(i-1))]#选择节点范围
    a.Set(nodes=nodes1, name='SE%s'%i)

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(adaptiveMeshConstraints=ON)

    regionDef=mdb.models[s].rootAssembly.sets['SE%s'%i]
    #创建x方向的历程输出请求
    mdb.models[s].HistoryOutputRequest(name='XE%s'%i,
    createStepName='Step-1', variables=('COOR1', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)
    #获取刚刚创立的集合
    regionDef=mdb.models[s].rootAssembly.sets['SE%s'%i]
    #创建y方向的历程输出请求
    mdb.models[s].HistoryOutputRequest(name='YE%s'%i,
    createStepName='Step-1', variables=('COOR3', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)

#SF
point_9 = 126
point_10 = point_9 + diff
for i in range(1,points+1):
    a = mdb.models[s].rootAssembly
    n1 = a.instances['Ring-1'].nodes

    nodes1 = n1[(point_9-1+diff*(i-1)):(point_9+diff*(i-1))]#选择节点范围
    a.Set(nodes=nodes1, name='SF%s'%i)

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(adaptiveMeshConstraints=ON)

    regionDef=mdb.models[s].rootAssembly.sets['SF%s'%i]
    #创建x方向的历程输出请求
    mdb.models[s].HistoryOutputRequest(name='XF%s'%i,
    createStepName='Step-1', variables=('COOR1', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)
    #获取刚刚创立的集合
    regionDef=mdb.models[s].rootAssembly.sets['SF%s'%i]
    #创建y方向的历程输出请求
    mdb.models[s].HistoryOutputRequest(name='YF%s'%i,
    createStepName='Step-1', variables=('COOR3', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)

#SG
point_11 = 181
point_12 = point_11 + diff
for i in range(1,points+1):
    a = mdb.models[s].rootAssembly
    n1 = a.instances['Ring-1'].nodes

    nodes1 = n1[(point_11-1+diff*(i-1)):(point_11+diff*(i-1))]#选择节点范围
    a.Set(nodes=nodes1, name='SG%s'%i)

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(adaptiveMeshConstraints=ON)

    regionDef=mdb.models[s].rootAssembly.sets['SG%s'%i]
    #创建x方向的历程输出请求
    mdb.models[s].HistoryOutputRequest(name='XG%s'%i,
    createStepName='Step-1', variables=('COOR1', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)
    #获取刚刚创立的集合
    regionDef=mdb.models[s].rootAssembly.sets['SG%s'%i]
    #创建y方向的历程输出请求
    mdb.models[s].HistoryOutputRequest(name='YG%s'%i,
    createStepName='Step-1', variables=('COOR3', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)

  #SH  
point_13 = 186
point_14 = point_13 + diff
for i in range(1,points+1):
    a = mdb.models[s].rootAssembly
    n1 = a.instances['Ring-1'].nodes

    nodes1 = n1[(point_13-1+diff*(i-1)):(point_13+diff*(i-1))]#选择节点范围
    a.Set(nodes=nodes1, name='SH%s'%i)

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(adaptiveMeshConstraints=ON)

    regionDef=mdb.models[s].rootAssembly.sets['SH%s'%i]
    #创建x方向的历程输出请求
    mdb.models[s].HistoryOutputRequest(name='XH%s'%i,
    createStepName='Step-1', variables=('COOR1', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)
    #获取刚刚创立的集合
    regionDef=mdb.models[s].rootAssembly.sets['SH%s'%i]
    #创建y方向的历程输出请求
    mdb.models[s].HistoryOutputRequest(name='YH%s'%i,
    createStepName='Step-1', variables=('COOR3', ), numIntervals=200,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)
