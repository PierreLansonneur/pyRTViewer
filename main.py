# -*- coding: utf-8 -*-
############################################
# Visualize dosimetry files with this viewer
#       launch with : python main.py
# 		   ---
#           P. Lansonneur 2019
############################################

import matplotlib.pyplot as P
import numpy as np
import SimpleITK as sitk
import pydicom, pydicom.uid
from pydicom.dataset import Dataset, FileDataset
import datetime, time
import os
import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="No contour levels were found")
warnings.filterwarnings("ignore", message="converting a masked element to nan")

from matplotlib.path import Path
import matplotlib.patches as patches
from scipy.ndimage import rotate, zoom, map_coordinates
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

import sys
from fbs_runtime.application_context.PyQt5 import ApplicationContext
from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui, Qt
from ROI import ROISet
    
class RTGeneralVolume():
    
    def __init__(self):

        self.open = False					        # is opened?
        self.dim_x = self.dim_y = self.dim_z = 100		        # scan dimensions
        self.volume = np.zeros((self.dim_x,self.dim_y,self.dim_z))	# volume
        self.spacing = [1, 1, 1]				        # spacing
        self.origin = [0, 0, 0]					        # origin
        self.filename = None                                            # filename

    def image(self, ax, slice_=0, rot=False):

        if(len(np.shape(self.volume))==3):
          
            if(ax==1):   im = self.volume[slice_,:,:].T
            if(ax==2):   im = self.volume[:,slice_,:].T
            if(ax==3):   im = self.volume[:,:,slice_].T

        if(len(np.shape(self.volume))==2):  im = self.volume[:,:]

        if rot: im = np.rot90(im)
        
        return im
        
    def extent_(self, ax, rot=False):

        if(len(np.shape(self.volume))==3):
            if(ax==1):   ext_=np.array([self.origin[1],self.origin[1]-self.dim_y*self.spacing[1],self.origin[2]-self.dim_z*self.spacing[2],self.origin[2]])
            if(ax==2):   ext_=np.array([self.origin[0],self.origin[0]+self.dim_x*self.spacing[0],self.origin[2]-self.dim_z*self.spacing[2],self.origin[2]])
            if(ax==3):   ext_=np.array([self.origin[0],self.origin[0]+self.dim_x*self.spacing[0],self.origin[1]-self.dim_y*self.spacing[1],self.origin[1]])

        if(len(np.shape(self.volume))==2):  ext_=np.array([self.origin[1],self.origin[1]-self.dim_y*self.spacing[1],self.origin[0]-self.dim_x*self.spacing[0],self.origin[0]])
        
        if rot: ext_ = np.array([ext_[3],ext_[2],ext_[0],ext_[1]])
        #if rot: ext_ = [ext[2],ext[3],ext[1],ext[0]] #dosi ??
            
        return ext_
            
class RTCTVolume(RTGeneralVolume):
    
    def __init__(self):
        
        super().__init__()
        self.x1_lim = self.x2_lim = self.x3_lim =(0,self.dim_x)	        # projections range
        self.y1_lim = self.y2_lim = self.y3_lim = (0,self.dim_x)	# projections range
        #self.im1 = self.im2 = self.im3 = self.im = self.volume[0,:,:]   # scan projections
        #self.ext7 = self.ext8 = self.ext9 = self.extent = [0,100,0,100] # scan projections dimensions
        self.colormap = P.get_cmap('Greys')                             # colormap
        
class RTDosiVolume(RTGeneralVolume):
    
    def __init__(self):
        
        super().__init__()
        self.levels = np.array([0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]) 	# isodose levels
        self.D_PTV = 1 #np.nan                                          # PTV average dose
        self.colormap = P.get_cmap('jet')                               # colormap
        self.show = False                                               # show isodoses
        self.showOption = 1                                             # display option (1:,2:,3:)
        
    def PrintColorMap(self):
        print(self.colormap)

    def slice_(self, ct, ax, slice_ct):
        """
	return the slice along axis 'ax' in the dosimetry volume 
	corresponding to the slice "slice_ct" in the CT volume
	"""
        if(ax==1):  condition = (ct.origin[0] + slice_ct*ct.spacing[0] >= self.origin[0] + self.dim_x*self.spacing[0])
        if(ax==2):  condition = (ct.origin[1] + slice_ct*ct.spacing[1] >= self.origin[1] + self.dim_y*self.spacing[1])
        if(ax==3):  condition = (ct.origin[2] + slice_ct*ct.spacing[2] >= self.origin[2] + self.dim_z*self.spacing[2])
        
        if (ct.origin[ax-1] + slice_ct*ct.spacing[ax-1] <= self.origin[ax-1]):  return 0
        if condition:   return -1
        else:
            if(ax==1):  return int((slice_ct-int((self.origin[ax-1] - ct.origin[ax-1])/ct.spacing[ax-1]))*ct.spacing[ax-1]/self.spacing[ax-1])
            if(ax==2)or(ax==3): return int((slice_ct+int((self.origin[ax-1] - ct.origin[ax-1])/ct.spacing[ax-1]))*ct.spacing[ax-1]/self.spacing[ax-1])
            
class RTMainWindow(QMainWindow):
   
    from OpenFile import OpenFile, OpenDicomSerie, OpenDosi, OpenROI
    from ROI import ROISet

    def __init__(self):
        
        self.myCTVolume = RTCTVolume()
        self.myDosiVolume = RTDosiVolume()
        self.w1_moved = self.w2_moved = self.w3_moved = True
        self.rot1 = self.rot2 = self.rot3 = True
        self.inv_scale = False
        #self.c_scale = 'HOUNSFIELD'
        self.c_scale = 'AUTO'
        
        super().__init__()
        self.initUI()
        
    def initUI(self):
        
        global sb        
        self.setWindowTitle('PyRTViewer')              
        self.resize(700, 750)
        sb = self.statusBar()
        self.createGridLayout()
        self.createMenuBar()

    def createMenuBar(self):
        
        menubar = self.menuBar()
        ### File menu
        fileMenu = menubar.addMenu('File')
        fileMenu.addAction(QAction('Open', self, shortcut="Ctrl+O", triggered=self.OpenFile))        
        fileMenu.addAction(QAction('Open a DICOM serie', self, triggered=self.OpenDicomSerie))
        fileMenu.addAction(QAction('Import dosimetry', self, shortcut="Ctrl+I", triggered=self.OpenDosi))
        fileMenu.addAction(QAction('Import structures...', self, triggered=self.OpenROI))
        #fileMenu.addAction(QAction('Import a RP file...', self))
        fileMenu.addSeparator()
        fileMenu.addAction(QAction('Save as...', self))
        fileMenu.addAction(QAction('Exit', self, shortcut="Ctrl+Q",triggered=QApplication.instance().quit))

        ### Edit menu
        from functools import partial
        
        editMenu = menubar.addMenu('Edit')
        rot_menu = QMenu('Rotate', self)
        rot_menu.addAction(QAction('ax1', self, triggered=partial(self.rotateAxes,1)))
        rot_menu.addAction(QAction('ax2', self, triggered=partial(self.rotateAxes,2)))
        rot_menu.addAction(QAction('ax3', self, triggered=partial(self.rotateAxes,3)))
        editMenu.addMenu(rot_menu)
        
        ### Scale menu
        scalemenu = menubar.addMenu('Scale')
        ag = QActionGroup(self, exclusive=True)
        scalemenu.addAction(ag.addAction(QAction('Auto', self, checkable=True, triggered=partial(self.SetScale,'AUTO'))))
        scalemenu.addAction(ag.addAction(QAction('Hounsfield', self, checkable=True, triggered=partial(self.SetScale,'HOUNSFIELD'))))
        scalemenu.addSeparator()
        scalemenu.addAction(QAction('Invert scale', self, triggered=self.InvertScale))
        #scalemenu.invoke(0)

        ### Tool menu
        toolmenu = menubar.addMenu('Tools')
        toolmenu.addAction(QAction('Isodoses', self))
        toolmenu.addAction(QAction('Profile', self))
        toolmenu.addAction(QAction('Histogram', self))
        toolmenu.addAction(QAction('Gamma', self))
        toolmenu.addAction(QAction('Iso-surface', self))

        ### ROI menu
        roimenu = menubar.addMenu('Structures')
        roimenu.addAction(QAction('Compute DVH...', self, enabled=True))
        roimenu.addAction(QAction('show ROI', self, triggered=self.showROI))
        roimenu.addAction(QAction('Normalize dosi to PTV dose', self))

        ### Infos action
        menubar.addAction(QAction('Infos', self))
        
    def createGridLayout(self):
        global w1, w2, w3, fc1, fc2, fc3, vtkWidget
        
        self.horizontalGroupBox = QGroupBox()
        layout = QGridLayout()
        layout.setColumnStretch(0, 10);
        layout.setColumnStretch(1, 10);
        fc1 = FigureCanvasQTAgg(Figure())
        fc2 = FigureCanvasQTAgg(Figure())
        fc3 = FigureCanvasQTAgg(Figure())

        for fig in [fc1, fc2, fc3]:
            ax = fig.figure.gca()
            ax.axis('equal')
            ax.set_xticks([])
            ax.set_yticks([])
            for feat in ['top','bottom','left','right']:    ax.spines[feat].set_linewidth(0)
            fig.figure.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    
        vtkWidget = QVTKWidget()
        w1 = QSlider(QtCore.Qt.Horizontal, minimum = 0, maximum = self.myCTVolume.dim_x)
        w2 = QSlider(QtCore.Qt.Horizontal, minimum = 0, maximum = self.myCTVolume.dim_y)
        w3 = QSlider(QtCore.Qt.Horizontal, minimum = 0, maximum = self.myCTVolume.dim_z)

        w1.setValue(int(self.myCTVolume.dim_x/2))
        w2.setValue(int(self.myCTVolume.dim_y/2))
        w3.setValue(int(self.myCTVolume.dim_z/2))
        
        w1.valueChanged.connect(self.w1move)
        w2.valueChanged.connect(self.w2move)
        w3.valueChanged.connect(self.w3move)

##        check1 = QCheckBox("isodoses")
##        check1.setChecked(False)
##        check1.toggled.connect(self.show_isodose)
##
##        check2 = QCheckBox("structures")
##        check2.setChecked(False)
##        check2.toggled.connect(self.show_ROI)

        #pbar = QProgressBar(self)
        #pbar.setValue(50.)
             
        layout.addWidget(fc1,1,0)
        layout.addWidget(fc2,1,1)
        layout.addWidget(fc3,2,0)
        layout.addWidget(w1,0,0)
        layout.addWidget(w2,0,1)
        layout.addWidget(w3,3,0)
        layout.addWidget(vtkWidget,2,1)
        #layout.addWidget(pbar, 3, 1)
        #layout.addWidget(check1, 4, 0)
        #layout.addWidget(check2, 4, 1)        

        self.horizontalGroupBox.setLayout(layout)
        wid = QWidget(self)
        self.setCentralWidget(wid)
        windowLayout = QVBoxLayout()
        windowLayout.addWidget(self.horizontalGroupBox)
        wid.setLayout(windowLayout)

        fc1.figure.canvas.mpl_connect('scroll_event', self.on_mousewheel1)
        fc2.figure.canvas.mpl_connect('scroll_event', self.on_mousewheel2)
        fc3.figure.canvas.mpl_connect('scroll_event', self.on_mousewheel3)

        #if(self.myDosiVolume.showOption==3):    self.cbaxes = fc1.figure.add_axes([0.78,0.67,0.02,0.31])

    def show_isodose(self):
        print("show_isodose: " + str(check1.isChecked()))

    def show_ROI(self):
        check2 = self.sender()
        print("show_ROI: " + str(check2.isChecked()))        
        
    def on_mousewheel1(self,event):
        w1.setValue(w1.value()-event.step)
	
    def on_mousewheel2(self,event):
        w2.setValue(w2.value()-event.step)
        
    def on_mousewheel3(self,event):
        w3.setValue(w3.value()-event.step)
        
    def w1move(self):
        self.w1_moved = True
        self.update()
        self.w1_moved = False
	
    def w2move(self):
        self.w2_moved = True
        self.update()
        self.w2_moved = False
        
    def w3move(self):
        self.w3_moved = True
        self.update()
        self.w3_moved = False

    def SetScales(self):
        """
	Set the limits of sliders
	"""
        w1.setMaximum(self.myCTVolume.dim_x-1)
        w2.setMaximum(self.myCTVolume.dim_y-1)
        w3.setMaximum(self.myCTVolume.dim_z-1)
        w1.setValue(int(self.myCTVolume.dim_x/2))
        w2.setValue(int(self.myCTVolume.dim_y/2))
        w3.setValue(int(self.myCTVolume.dim_z/2))

    def Clear_axes(self, ClearAll = False):
        """
	Clear axes
	"""
	#if ((ROI_open==True)and(ROI_show==True)) or ((dosi_open == True)and(isodose_show == True)):
        if self.w1_moved:  fc1.figure.gca().clear()
        if self.w2_moved:  fc2.figure.gca().clear()                        
        if self.w3_moved:  fc3.figure.gca().clear()
            
        if(ClearAll == True):
            for fig in [fc1, fc2, fc3]: fig.figure.gca().clear()

    def UpdateAll(self):
        """
	Update each axes
	"""
        self.w1_moved = self.w2_moved = self.w3_moved = True
        self.update()
        self.statusBar().showMessage('')
        self.w1_moved = self.w2_moved = self.w3_moved = False
	
    def update(self):
        
        #print(w1.value(),w2.value(),w3.value())

        ax1 = fc1.figure.gca()
        ax2 = fc2.figure.gca()
        ax3 = fc3.figure.gca()

        self.Clear_axes()
        
        ### CT        
        if(len(np.shape(self.myCTVolume.volume))==3):

            if self.w1_moved:	ax1.imshow(self.myCTVolume.image(1, w1.value(), self.rot1), extent=self.myCTVolume.extent_(1, self.rot1), cmap=self.myCTVolume.colormap)
            if self.w2_moved:	ax2.imshow(self.myCTVolume.image(2, w2.value(), self.rot2), extent=self.myCTVolume.extent_(2, self.rot2), cmap=self.myCTVolume.colormap)
            if self.w3_moved:	ax3.imshow(self.myCTVolume.image(3, w3.value(), self.rot3), extent=self.myCTVolume.extent_(3, self.rot3), cmap=self.myCTVolume.colormap)            
	    
        if(len(np.shape(self.myCTVolume.volume))==2):   ax1.imshow(self.myCTVolume.image(1, 0, self.rot1), extent=self.myCTVolume.extent_(1, self.rot1), cmap=self.myCTVolume.colormap)


		### profile tool
##		if profile_show:
##			profile_ax.plot([xp1, xp2],[yp1, yp2], 'bo-')
##			Update_profile()

        ### Dosimetry
        if self.myDosiVolume.open and self.myDosiVolume.show:	

            dosemap = self.myDosiVolume.colormap
            levels = self.myDosiVolume.levels
            option = self.myDosiVolume.showOption
            D_PTV = self.myDosiVolume.D_PTV
            
            if self.w1_moved:
                dos1 = self.myDosiVolume.image(1, self.myDosiVolume.slice_(self.myCTVolume, 1, w1.value())-1, self.rot1)
                ext1_dosi = self.myDosiVolume.extent_(1, self.rot1)
                if(option==1):  cs1 = ax1.contour(np.flipud(dos1), np.nanmax(dos1)*levels, cmap = dosemap, linewidths=1, extent=ext1_dosi)
                if(option==2):  cs1 = ax1.contour(np.flipud(dos1), D_PTV*levels, cmap = dosemap, linewidths=1, extent=ext1_dosi)
                if(option==3):  cs1 = ax1.imshow(np.ma.masked_where(dos1<0.05*np.nanmax(self.myDosiVolume.volume),dos1), alpha=0.5, cmap = dosemap, extent=ext1_dosi)

            if self.w2_moved:
                dos2 = self.myDosiVolume.image(2, self.myDosiVolume.slice_(self.myCTVolume, 2, w2.value())-1, self.rot2)
                ext2_dosi = self.myDosiVolume.extent_(2, self.rot2)
                if(option==1):  cs2 = ax2.contour(np.flipud(dos2), np.nanmax(dos2)*levels, cmap = dosemap, linewidths=1, extent=ext2_dosi)
                if(option==2):  cs2 = ax2.contour(np.flipud(dos2), D_PTV*levels, cmap = dosemap, linewidths=1, extent=ext2_dosi)
                if(option==3):  cs2 = ax2.imshow(np.ma.masked_where(dos2<0.05*np.nanmax(self.myDosiVolume.volume),dos2), alpha=0.5, cmap = dosemap, extent=ext2_dosi)

            if self.w3_moved:
                dos3 = self.myDosiVolume.image(3, self.myDosiVolume.slice_(self.myCTVolume, 3, w3.value())-1, self.rot3)
                ext3_dosi = self.myDosiVolume.extent_(3, self.rot3)
                if(option==1):  cs3 = ax3.contour(np.flipud(dos3), np.nanmax(dos3)*levels, cmap = dosemap, linewidths=1, extent=ext3_dosi)
                if(option==2):  cs3 = ax3.contour(np.flipud(dos3), D_PTV*levels, cmap = dosemap, linewidths=1, extent=ext3_dosi)
                if(option==3):  cs3 = ax3.imshow(np.ma.masked_where(dos3<0.05*np.nanmax(self.myDosiVolume.volume),dos3), alpha=0.5, cmap = dosemap, extent=ext3_dosi)
                
            if(option==1 or 2):
                try:
                    #for g in range(len(cs1.collections)):	cs1.collections[g].set_label('{0} Gy'.format(int(100*levels[g])))
                    for g in range(len(cs1.collections)):	cs1.collections[g].set_label('{0} %'.format(int(100*levels[g])))
                    if(self.w1_moved and option==1):	leg = ax1.legend(frameon=False,title="$D_{max}$")
                    if(self.w1_moved and option==2):	leg = ax1.legend(frameon=False,title="$D_{PTV}$")
                    if self.inv_scale and self.w1_moved:	
                            for text in leg.get_texts():	text.set_color("white")
                            leg.get_title().set_color("white")
                except Exception:   pass

            if(option==3):
                try:    
                        cbar = P.colorbar(cs1, cax=cbaxes)
                        if self.inv_scale and w1_moved:      cbaxes.tick_params(colors = "white")#,grid_color = "white")
                except Exception:	pass

	### Structures
        try:
            if self.w1_moved:
                
                for ROI_index in range(self.myROISet.N_ROI):
                    
                    x = self.myROISet.infos[ROI_index,1]
                    y = self.myROISet.infos[ROI_index,2]
                    z = self.myROISet.infos[ROI_index,3]
                    N_vert_cumul = self.myROISet.infos[ROI_index,4]
                    ROI_color = self.myROISet.infos[ROI_index,5]

                    if self.rot1: 	y,x = x,y

                    z_s = ((z - self.myCTVolume.origin[0])/self.myCTVolume.spacing[0]).astype(int) # scale vertices coordinates to voxel unit

                    if not (w1.value() in z_s):   continue
                    
                    index = np.where(z_s == w1.value()) # indices of the slice concerned
                    
                    if(len(N_vert_cumul)>1): ### surface contours
                        for count in range(len(index[0])):
                            if (index[0][count] in N_vert_cumul):	
                                mini = N_vert_cumul[np.where(N_vert_cumul==index[0][count])[0][0]-1]
                                maxi = N_vert_cumul[np.where(N_vert_cumul==index[0][count])[0][0]]  # get 1st and last indices of the contour,
                                if (maxi<=mini):    continue
                                vertices = np.append(y[mini:maxi],x[mini:maxi]).reshape(2,-1).T # get the vertices coordinates of the contour,
                                path = Path(vertices, closed=True)				# make them a path,
                                patch = patches.PathPatch(path, facecolor='none', edgecolor=ROI_color, lw = 1) #zorder=3 to display in front
                                ax1.add_patch(patch)                                                        

                    elif(len(N_vert_cumul)==1):     ax1.plot([y[0]],[x[0]],color=ROI_color,marker = '*',zorder=3) ### points

        except Exception:	pass
        
	#Set_axes_lim()
       
        if self.c_scale=='HOUNSFIELD':  self.SetDisplayRange(-1000,2000)
        if self.c_scale=='USER':	self.SetDisplayRange(0,30)

        if self.w1_moved:   fc1.draw()
        if self.w2_moved:   fc2.draw()
        if self.w3_moved:   fc3.draw()

    def SetScale(self, scale):
        self.c_scale=scale
        self.UpdateAll()
        
    def SetDisplayRange(self, minI, maxI):
        """ 
	Set axes colorscale 
	>>>Usage:
	SetDisplayRange(-1000, 2000) #Hounsfield scale
	"""
        for fig in [fc1, fc2, fc3]:
            for im in fig.figure.gca().get_images():#    im.set_clim(minI, maxI)
                if self.myDosiVolume.open and self.myDosiVolume.show:
                    if(len(fig.figure.gca().get_images())>0):   fig.figure.gca().get_images()[0].set_clim(minI, maxI)
            else:
                for im in fig.figure.gca().get_images():    im.set_clim(minI, maxI)

    def InvertScale(self):
        """
        Reverse grayscale (CT volume)
	"""
        if(self.myCTVolume.colormap==P.get_cmap('Greys')):
            self.myCTVolume.colormap = P.get_cmap('Greys_r')
            for fig in [fc1, fc2, fc3]: fig.figure.gca().set_facecolor('0')
            vtkWidget.ren.SetBackground(0,0,0) # black background
            
        else:
            self.myCTVolume.colormap = P.get_cmap('Greys')
            for fig in [fc1, fc2, fc3]: fig.figure.gca().set_facecolor('1')
            vtkWidget.ren.SetBackground(1,1,1) # white background
            
        self.UpdateAll()
        vtkWidget.update()

    def showROI(self):
        
        for i in range(0,self.myROISet.N_ROI):  vtkWidget.addROI(self.myROISet, i)

        # set camera focal point
        x = self.myROISet.infos[0,1]
        y = self.myROISet.infos[0,2]
        z = self.myROISet.infos[0,3]
        vtkWidget.camera.SetFocalPoint( (np.max(x)+np.min(x))/2., (np.max(y)+np.min(y))/2., (np.max(z)+np.min(z))/2. )
        vtkWidget.ren.ResetCamera()
        
    def rotateAxes(self, ax):
        if(ax==1):  self.rot1 = not self.rot1
        if(ax==2):  self.rot2 = not self.rot2
        if(ax==3):  self.rot3 = not self.rot3
        self.UpdateAll()
        
class QVTKWidget(QVTKRenderWindowInteractor):
    
    def __init__(self):

        QVTKRenderWindowInteractor.__init__(self)
        self.ren = vtk.vtkRenderer()
        self.GetRenderWindow().AddRenderer(self.ren)
        self.ren.SetBackground(1, 1, 1) # white background
   
        self.camera = self.ren.GetActiveCamera()
        self.update()

    def createSphereSource(self):
        
        # Create source
        source = vtk.vtkSphereSource()
        source.SetCenter(0, 0, 0)
        source.SetRadius(5.0)

        # Create a mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())

        # Create an actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        self.ren.AddActor(actor)

    def addROI(self, ROISet, ROI_index):
        
        x = ROISet.infos[ROI_index,1]
        y = ROISet.infos[ROI_index,2]
        z = ROISet.infos[ROI_index,3]
        N_vert_cumul = ROISet.infos[ROI_index,4]
        ROI_color = ROISet.infos[ROI_index,5]
        
        if(len(N_vert_cumul)<=1):   return ### only surface contours

        for i in range(0,len(N_vert_cumul)-1):
            
            mini = N_vert_cumul[i]
            maxi = N_vert_cumul[i+1]  # get 1st and last indices of each contours
            
            if (maxi<=mini): continue
        
            pd = vtk.vtkPolyData()
            points = vtk.vtkPoints()
            lines = vtk.vtkCellArray()

            for j in range(0, maxi-mini):    points.InsertPoint(j, x[mini:maxi][j], y[mini:maxi][j], z[mini:maxi][j])

            vertex_indices = list(range(0, maxi-mini))
            vertex_indices.append(0) # close the contour
            
            lines.InsertNextCell(maxi-mini + 1, vertex_indices)    
        
            pd.SetPoints(points)
            pd.SetLines(lines)

            Mapper = vtk.vtkPolyDataMapper()
            Mapper.SetInputDataObject(pd)
        
            actor = vtk.vtkActor()
            actor.SetMapper(Mapper)
            actor.GetProperty().SetColor(ROI_color) 
            
            self.ren.AddActor(actor)  # add vtkActor for each contour                                                  

        #elif(len(N_vert_cumul)==1):     ax1.plot([y[0]],[x[0]],color=ROI_color,marker = '*',zorder=3) ### points
       
    def update(self):
        
        self.iren = self.GetRenderWindow().GetInteractor()

##        axes = vtk.vtkAxesActor()
##        tmp = vtk.vtkOrientationMarkerWidget()
####        rgba = [0] * 4
####        colors.GetColor("Carrot", rgba)
####        widget.SetOutlineColor(rgba[0], rgba[1], rgba[2])
##        tmp.SetOrientationMarker(axes)
##        tmp.SetInteractor(self.iren)
##        tmp.SetViewport(0.0, 0.0, 0.4, 0.4)
##        tmp.SetEnabled(1)
##        tmp.InteractiveOn()
        
        style = vtk.vtkInteractorStyleTrackballCamera()
        self.iren.SetInteractorStyle(style)
        self.iren.Initialize()
        self.iren.Start()


        
if __name__ == '__main__':
    appctxt = ApplicationContext()       # 1. Instantiate ApplicationContext
    window = RTMainWindow()
    window.show()
    #window.showMaximized()
    exit_code = appctxt.app.exec_()      # 2. Invoke appctxt.app.exec_()
    sys.exit(exit_code)
