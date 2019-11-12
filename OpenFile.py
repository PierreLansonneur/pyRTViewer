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
from PyQt5.QtWidgets import QFileDialog
from PyQt5 import QtCore, QtGui, Qt
    
def OpenFile(self, filepath=None):

    self.statusBar().showMessage('Opening file ...')

    if(filepath==False):
        filepath, _ = QFileDialog.getOpenFileNames(self, "Open file","./","All file (*.dcm *.mhd)")
        filepath = filepath[0]
        filename = QtCore.QFileInfo(filepath[0]).fileName()
        filedir = QtCore.QFileInfo(filepath[0]).path() # +'/'
        #print(dir_ini,tmp)

    ### .dcm file ###
    if(filepath.endswith('.dcm')==True):    OpenDicomFile(self, filepath)
        
def OpenDicomFile(self, filepath):
        
    ds = pydicom.read_file(filepath)
    ds.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian 
    self.myCTVolume.volume = ds.pixel_array

    try:
        self.myCTVolume.spacing[0:1] = ds.PixelSpacing
        self.myCTVolume.origin[0:1] = ds.ImagePositionPatient
        
    except Exception:
        self.myCTVolume.spacing = ds.PixelSpacing
        self.myCTVolume.origin = ds.ImagePositionPatient

    if (ds.Modality == 'RTDOSE'):
        if ("DoseGridScaling" in ds):	self.myCTVolume.volume = float(ds.DoseGridScaling)*self.myCTVolume.volume
        
    else:
        if ("RescaleSlope" in ds):	self.myCTVolume.volume = float(ds.RescaleSlope)*self.myCTVolume.volume
        if ("RescaleIntercept" in ds):	self.myCTVolume.volume = self.myCTVolume.volume + float(ds.RescaleIntercept)

    if(len(np.shape(self.myCTVolume.volume))==3):
        self.myCTVolume.spacing = [ float(ds.GridFrameOffsetVector[1] - ds.GridFrameOffsetVector[0]), float(self.myCTVolume.spacing[1]),float(self.myCTVolume.spacing[0])]
        self.myCTVolume.origin = [float(self.myCTVolume.origin[2]),float(self.myCTVolume.origin[1]),float(self.myCTVolume.origin[0])]

    ct_swapZ =(ds.ImageOrientationPatient[0:3] == [1, 0, 0])
    ct_swapY =(ds.ImageOrientationPatient[3:6] == [0, 1, 0])

    if(len(np.shape(self.myCTVolume.volume))==3):
        self.myCTVolume.dim_x = np.shape(self.myCTVolume.volume)[0]
        self.myCTVolume.dim_y = np.shape(self.myCTVolume.volume)[1]
        self.myCTVolume.dim_z = np.shape(self.myCTVolume.volume)[2]

        ### Dealing with image orientation
        #print '  ct_swapY, ct_swapZ :', ct_swapY, ct_swapZ

        if(ct_swapY == True):
            # flip volume, Y direction
            self.myCTVolume.volume = np.flip(self.myCTVolume.volume,1) 
            self.myCTVolume.origin[1] = self.myCTVolume.origin[1] + self.myCTVolume.dim_y*self.myCTVolume.spacing[1]
            
        if(ct_swapZ == True):
            # flip volume, Z direction
            self.myCTVolume.volume = np.flip(self.myCTVolume.volume,2) 
            self.myCTVolume.origin[2] = self.myCTVolume.origin[2] + self.myCTVolume.dim_z*self.myCTVolume.spacing[2]
            
        if(ct_swapZ == True)and(ct_swapY == True):      self.myCTVolume.spacing[1], self.myCTVolume.spacing[2] = self.myCTVolume.spacing[2], self.myCTVolume.spacing[1]
        
    if(len(np.shape(self.myCTVolume.volume))==2):
        self.myCTVolume.dim_x = np.shape(self.myCTVolume.volume)[0]
        self.myCTVolume.dim_y = np.shape(self.myCTVolume.volume)[1]
        self.myCTVolume.dim_z = 0

    self.myCTVolume.open = True
    
##	Set_axes_lim_init()
##	Set_scales()
    self.update()
        
    self.statusBar().showMessage('file successfully opened!')
