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
#import hello
    
class RTGeneralVolume():
    
    def __init__(self):

        self.open = False					        # is opened?
        self.dim_x = self.dim_y = self.dim_z = 100		        # scan dimensions
        self.volume = np.zeros((self.dim_x,self.dim_y,self.dim_z))	# volume
        self.spacing = [1, 1, 1]				        # spacing
        self.origin = [0, 0, 0]					        # origin
        self.filename = None                                            # filename

       
class RTCTVolume(RTGeneralVolume):
    
    def __init__(self):
        
        super().__init__()
        self.x1_lim = self.x2_lim = self.x3_lim =(0,self.dim_x)	        # projections range
        self.y1_lim = self.y2_lim = self.y3_lim = (0,self.dim_x)	# projections range
        self.im1 = self.im2 = self.im3 = self.im = self.volume[0,:,:]   # scan projections
        self.ext7 = self.ext8 = self.ext9 = self.extent = [0,100,0,100] # scan projections dimensions
        self.colormap = P.get_cmap('Greys')                             # colormap
        
class RTDosiVolume(RTGeneralVolume):
    
    def __init__(self):
        
        super().__init__()
        self.isodose_show = False					# show isodoses
        self.levels = np.array([0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]) 	# isodose levels
        #self.Dvar = IntVar(window)                                     # dose display method
        #self.Dvar.set(1)                                               # normalize to the max
        self.coeff_D_PTV = np.nan                                       # PTV average dose
        self.colormap = P.get_cmap('jet')                               # colormap

    def PrintColorMap(self):
        print(self.colormap)
        
class RTMainWindow(QMainWindow):
   
    from OpenFile import OpenFile
    
    def __init__(self):
        
        self.myCTVolume = RTCTVolume()
        self.myDosiVolume = RTDosiVolume()
        self.w1_moved = self.w2_moved = self.w3_moved = True
        self.rot1 = self.rot2 = self.rot3 = True
        super().__init__()
        self.initUI()
        
    def initUI(self):
        
        global sb        
        self.setWindowTitle('PyRTViewer')              
        self.resize(700, 750)
        sb = self.statusBar()
        self.createMenuBar()
        self.createGridLayout()

    def createMenuBar(self):
        
        menubar = self.menuBar()
        ### File menu
        fileMenu = menubar.addMenu('File')
        fileMenu.addAction(QAction('Open', self, shortcut="Ctrl+O", triggered=self.OpenFile))        
        fileMenu.addAction(QAction('Open a DICOM serie', self, triggered=RTDosiVolume.PrintColorMap))
        fileMenu.addAction(QAction('Import dosimetry', self))
        fileMenu.addAction(QAction('Import a RP file...', self))
        fileMenu.addSeparator()
        fileMenu.addAction(QAction('Save as...', self))
        fileMenu.addAction(QAction('Exit', self, shortcut="Ctrl+Q",triggered=QApplication.instance().quit))

        ### Edit menu
        editMenu = menubar.addMenu('Edit')
        revX_menu = QMenu('Flip X axis', self)
        revX_menu.addAction(QAction('ax1', self))
        revX_menu.addAction(QAction('ax2', self))
        revX_menu.addAction(QAction('ax3', self))
        revY_menu = QMenu('Flip Y axis', self)
        revY_menu.addAction(QAction('ax1', self))
        revY_menu.addAction(QAction('ax2', self))
        revY_menu.addAction(QAction('ax3', self))
        rot_menu = QMenu('Rotate', self)
        rot_menu.addAction(QAction('ax1', self))
        rot_menu.addAction(QAction('ax2', self))
        rot_menu.addAction(QAction('ax3', self))
        
        editMenu.addMenu(revX_menu)
        editMenu.addMenu(revY_menu)
        editMenu.addMenu(rot_menu)
        
        ### Scale menu
        scalemenu = menubar.addMenu('Scale')
        ag = QActionGroup(self, exclusive=True)
        scalemenu.addAction(ag.addAction(QAction('Auto', self, checkable=True)))
        scalemenu.addAction(ag.addAction(QAction('Hounsfield', self, checkable=True)))
        scalemenu.addSeparator()
        scalemenu.addAction(QAction('Invert scale', self))
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
        roimenu.addAction(QAction('Import structures...', self))
        roimenu.addAction(QAction('Compute DVH...', self, enabled=True))
        roimenu.addAction(QAction('Crop dosimetry', self))
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
        #layout.addWidget(w1,0,0)
        #layout.addWidget(w2,0,1)
        #layout.addWidget(w3,3,0)
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

    def show_isodose(self):
        print("show_isodose: " + str(check1.isChecked()))

    def show_ROI(self):
        check2 = self.sender()
        print("show_ROI: " + str(check2.isChecked()))
        
        
    def on_mousewheel1(self,event):
        global w1
        w1.setValue(w1.value()-event.step)
	
    def on_mousewheel2(self,event):
        global w2
        w2.setValue(w2.value()-event.step)
        
    def on_mousewheel3(self,event):
        global w3
        w3.setValue(w3.value()-event.step)
        
    def w1move(self):        
        global w1_moved
        self.w1_moved = True
        self.update()
        self.w1_moved = False
	
    def w2move(self):        
        global w2_moved
        self.w2_moved = True
        self.update()
        self.w2_moved = False
        
    def w3move(self):        
        global w3_moved
        self.w3_moved = True
        self.update()
        self.w3_moved = False

    def update(self):
        
        print(w1.value(),w2.value(),w3.value())

        ax1 = fc1.figure.gca()
        ax2 = fc2.figure.gca()
        ax3 = fc3.figure.gca()

        ### CT        
        if(len(np.shape(self.myCTVolume.volume))==3): ### 3D files #############
            
            slice_ax1 = w1.value()
            slice_ax2 = w2.value()
            slice_ax3 = w3.value()

            im1 = self.myCTVolume.volume[slice_ax1,:,:].T
            im2 = self.myCTVolume.volume[:,slice_ax2,:].T
            im3 = self.myCTVolume.volume[:,:,slice_ax3].T

            ext7=[self.myCTVolume.origin[1],self.myCTVolume.origin[1]-self.myCTVolume.dim_y*self.myCTVolume.spacing[1],self.myCTVolume.origin[2]-self.myCTVolume.dim_z*self.myCTVolume.spacing[2],self.myCTVolume.origin[2]]
            ext8=[self.myCTVolume.origin[0],self.myCTVolume.origin[0]+self.myCTVolume.dim_x*self.myCTVolume.spacing[0],self.myCTVolume.origin[2]-self.myCTVolume.dim_z*self.myCTVolume.spacing[2],self.myCTVolume.origin[2]]
            ext9=[self.myCTVolume.origin[0],self.myCTVolume.origin[0]+self.myCTVolume.dim_x*self.myCTVolume.spacing[0],self.myCTVolume.origin[1]-self.myCTVolume.dim_y*self.myCTVolume.spacing[1],self.myCTVolume.origin[1]]

            if self.rot1:	im1, ext7 = np.rot90(im1), [ext7[3],ext7[2],ext7[0],ext7[1]]
            if self.rot2:	im2, ext8 = np.rot90(im2), [ext8[3],ext8[2],ext8[0],ext8[1]]
            if self.rot3:	im3, ext9 = np.rot90(im3), [ext9[3],ext9[2],ext9[0],ext9[1]]

            if self.w1_moved:	ax1.imshow(im1, extent=ext7, cmap=self.myCTVolume.colormap)
            if self.w2_moved:	ax2.imshow(im2, extent=ext8, cmap=self.myCTVolume.colormap)
            if self.w3_moved:	ax3.imshow(im3, extent=ext9, cmap=self.myCTVolume.colormap)
		
        if(len(np.shape(self.myCTVolume.volume))==2): ### 2D files #############
            
            im1 = self.myCTVolume.volume[:,:]
            ext7=[self.myCTVolume.origin[1],self.myCTVolume.origin[1]-self.myCTVolume.dim_y*self.myCTVolume.spacing[1],self.myCTVolume.origin[0]-self.myCTVolume.dim_x*self.myCTVolume.spacing[0],self.myCTVolume.origin[0]]
            if self.rot1:	im1, ext7 = np.rot90(im1), [ext7[3],ext7[2],ext7[0],ext7[1]]
            ax1.imshow(im1, extent=ext7, cmap=self.myCTVolume.colormap)


		### profile tool
##		if profile_show:
##			profile_ax.plot([xp1, xp2],[yp1, yp2], 'bo-')
##			Update_profile()

	#Set_axes_lim()
        for fig in [fc1, fc2, fc3]: fig.draw()
	
	#if c_scale=='HOUNSFIELD':	SetDisplayRange(-1000,2000)
	#if c_scale=='USER':		SetDisplayRange(0,30)

	
class QVTKWidget(QVTKRenderWindowInteractor):
    
    def __init__(self):

        QVTKRenderWindowInteractor.__init__(self)
        self.ren = vtk.vtkRenderer()
        self.GetRenderWindow().AddRenderer(self.ren)
        self.ren.SetBackground(1, 1, 1) # white background
        self.ren.ResetCamera()
        self.iren = self.GetRenderWindow().GetInteractor()
        
        #self.createSphereSource()
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
        
    def update(self):
        
        self.iren.Initialize()
        self.iren.Start()
        
if __name__ == '__main__':
    appctxt = ApplicationContext()       # 1. Instantiate ApplicationContext
    window = RTMainWindow()
    window.show()
    #window.showMaximized()
    exit_code = appctxt.app.exec_()      # 2. Invoke appctxt.app.exec_()
    sys.exit(exit_code)
