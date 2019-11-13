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
    
    #self.Set_axes_lim_init()
    self.SetScales()
    self.update()
        
    self.statusBar().showMessage('file successfully opened!')
    
def OpenDicomSerie(self, dirname=None):
    """
    Open a dicom serie
    """
    ct_swapY, ct_swapZ = False, False
    self.statusBar().showMessage('Opening DICOM serie ... ')

    # Opening file
    if(dirname==False):
        filepath, _ = QFileDialog.getOpenFileNames(self, "Open file","./","DICOM (*.dcm)")
        filepath = filepath[0]
        filename = QtCore.QFileInfo(filepath[0]).fileName()
        filedir = QtCore.QFileInfo(filepath[0]).path() # +'/'
        filelist = os.listdir(os.path.dirname(filepath))
        #print(filelist)
        
    else:
        filelist = os.listdir(dirname)
        filepath = dirname + filelist[0]

    #filename_CT = file_path
    #dir_ini = str(file_path.rsplit('/', 1)[0])+'/'

    # Getting dimensions
    ds = pydicom.read_file(filepath)
    sp = ds.PixelSpacing
    ds.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
    ct_swapZ =(ds.ImageOrientationPatient[0:3] == [1, 0, 0])
    ct_swapY =(ds.ImageOrientationPatient[3:6] == [0, 1, 0])

    self.myCTVolume.dim_x = 0
    for f in filelist:
        if f.endswith(".dcm"): self.myCTVolume.dim_x += 1        

    self.myCTVolume.dim_y, self.myCTVolume.dim_z = np.shape(ds.pixel_array)[1], np.shape(ds.pixel_array)[0]

    self.myCTVolume.volume = np.zeros((self.myCTVolume.dim_x, self.myCTVolume.dim_y,self.myCTVolume.dim_z))
    slicelocation = np.zeros(self.myCTVolume.dim_x)

    # creating volume
    for f,i in zip(filelist,range(self.myCTVolume.dim_x)):
        if f.endswith(".dcm"):
            ds = pydicom.read_file(os.path.dirname(filepath)+'/'+f)
            ds.file_meta.transfersyntaxuid = pydicom.uid.ImplicitVRLittleEndian 
            self.myCTVolume.volume[i,:,:] = ds.pixel_array
            if('slicelocation' in ds):  slicelocation[i] = ds.SliceLocation
            else:   slicelocation[i] = ds.ImagePositionPatient[2]

    order = np.argsort(slicelocation)
    slicelocation = slicelocation[order] # slicelocation is now sorted

    self.myCTVolume.spacing = [float(slicelocation[1] - slicelocation[0]),float(sp[1]), float(sp[0])]
    self.myCTVolume.origin = [float(slicelocation[0]),float(ds.ImagePositionPatient[1]),float(ds.ImagePositionPatient[0])]
    self.myCTVolume.volume = self.myCTVolume.volume[order,:,:] # volume is now sorted

    if ("RescaleSlope" in ds):	self.myCTVolume.volume *= float(ds.RescaleSlope)
    if ("RescaleIntercept" in ds):	self.myCTVolume.volume += float(ds.RescaleIntercept)

    # Dealing with image orientation
    print('ct_swapY, ct_swapZ :', ct_swapY, ct_swapZ)
    
    if ct_swapY:
        self.myCTVolume.volume = np.flip(self.myCTVolume.volume,1) # flip volume, Y direction
        self.myCTVolume.origin[1] += self.myCTVolume.dim_y*self.myCTVolume.spacing[1]

    if ct_swapZ:
        self.myCTVolume.volume = np.flip(self.myCTVolume.volume,2) # flip volume, Z direction
        self.myCTVolume.origin[2] += self.myCTVolume.dim_z*self.myCTVolume.spacing[2]
        
    if ct_swapZ and ct_swapY:   self.myCTVolume.spacing[1], self.myCTVolume.spacing[2] = self.myCTVolume.spacing[2], self.myCTVolume.spacing[1]

    #self.Set_axes_lim_init()
    self.SetScales()
    #CT_open = True
    self.UpdateAll()

    self.statusBar().showMessage('file successfully opened!')

def OpenDosi(self,filepath=None): 
    """
    Open a dosimetry file with a .dcm or .mhd extension
    """

    dosi_swapY,dosi_swapZ = False, False

    types = [('All files', '*.dcm *.mhd'), ('DCM files', '*.dcm'), ('MHD files', '*.mhd')]

    if(filepath==False):
        filepath, _ = QFileDialog.getOpenFileNames(self, "Open file","./","DICOM (*.dcm)")
        filepath = filepath[0]
        filename = QtCore.QFileInfo(filepath[0]).fileName()
        filedir = QtCore.QFileInfo(filepath[0]).path() # +'/'

    print('Opening RD file ...')

    ### .dcm file ###
    if(filepath.endswith('.dcm')):
        ds = pydicom.read_file(filepath)
        ds.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian 
        self.myDosiVolume.volume = float(ds.DoseGridScaling)*ds.pixel_array
        sp = ds.PixelSpacing
        self.myDosiVolume.spacing = [ float(ds.GridFrameOffsetVector[1] - ds.GridFrameOffsetVector[0]), float(sp[1]),float(sp[0])]
        self.myDosiVolume.origin = ds.ImagePositionPatient
        self.myDosiVolume.origin = [float(self.myDosiVolume.origin[2]),float(self.myDosiVolume.origin[1]),float(self.myDosiVolume.origin[0])]
        dosi_swapZ =(ds.ImageOrientationPatient[0:3] == [1, 0, 0])
        dosi_swapY =(ds.ImageOrientationPatient[3:6] == [0, 1, 0])
    
        self.myDosiVolume.dim_x = np.shape(self.myDosiVolume.volume)[0]
        self.myDosiVolume.dim_y = np.shape(self.myDosiVolume.volume)[1]        
        if(len(np.shape(self.myDosiVolume.volume))==3):   self.myDosiVolume.dim_z = np.shape(self.myDosiVolume.volume)[2]
        if(len(np.shape(self.myDosiVolume.volume))==2):   self.myDosiVolume.dim_z = 1

    #print 'dosi type', dosi.dtype
    
    # Dealing with image orientation
    if(dosi_swapY == True):
        self.myDosiVolume.volume = np.flip(self.myDosiVolume.volume,1) # flip volume
        self.myDosiVolume.origin[1] += self.myDosiVolume.dim_y*self.myDosiVolume.spacing[1]		
    if(dosi_swapZ == True):
        self.myDosiVolume.volume = np.flip(self.myDosiVolume.volume,2) # flip volume
        self.myDosiVolume.origin[2] += self.myDosiVolume.dim_z*self.myDosiVolume.spacing[2]
    if(dosi_swapY == True)and(dosi_swapZ == True):
        self.myDosiVolume.spacing[1], self.myDosiVolume.spacing[2] = self.myDosiVolume.spacing[2], self.myDosiVolume.spacing[1]

    print('dosi_swapY, dosi_swapZ :', dosi_swapY, dosi_swapZ)

    #dosi_open = True
    #isodose_show = True
    #check1.select()
    self.UpdateAll()

    self.statusBar().showMessage('file successfully opened!')
