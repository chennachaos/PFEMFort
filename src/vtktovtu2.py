import sys
from paraview.simple import *

    inputFile = "beam3D-2012020-nproc1.vtk"
    outputFile = inputFile[:-1] + 'u'
    print x,': Converting ', inputFile, '  ->  ', outputFile
    reader = LegacyVTKReader( FileNames= inputFile )
    writer = XMLUnstructuredGridWriter()
    writer.FileName = "beam3D-2012020-nproc1.vtu"
    writer.UpdatePipeline()
    Delete(reader)
    Delete(writer
