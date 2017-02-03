##------------------------------------------------------------------------------------------------------------------
##  Script Name: Stream Cross-Sections [CrossSections.py]
##
##  Creates cross-sections that are perpendicular to a stream channel inside the stream banks, but are
##  'perpendicular' (close-to perpendicular) to the valley walls in which the stream is in.
##
##  Author: Michael Harris
##  Date: 05/19/2016
##
##  Some code snippets and ideas were obtained from Mark Ellefson 'Classify Stream Type' script and
##  and 'Perpendicular Transects' script by Mateus Ferreira.
##------------------------------------------------------------------------------------------------------------------

#Import modules
import arcpy, math, sys, traceback

#Set Environments
arcpy.env.overwriteOutput = True
arcpy.env.XYResolution = "0.00001 Meters"
arcpy.env.XYTolerance = "0.0001 Meters"

##------------------------------------------------------------------------------------------------------------------
## draingeSQMI function calculates the drainage area using flowacc.tif, converts to square miles.
##------------------------------------------------------------------------------------------------------------------
def drainageSQMI(cellNumber):
    drainage = cellNumber*(cellSize**2)/(5280**2)
    return drainage

##------------------------------------------------------------------------------------------------------------------
## thresholdXSW calculates the typical bankfull width of a stream using a regional curve equation.
## This does not take into account the different stream types, obviously F streams will be wider
## than B stream types, but this is a ball park width for the tool to use based on the drainage area.
##------------------------------------------------------------------------------------------------------------------
def thresholdXSW(DA):
    DAstr = str(DA)
    if ratio == "":
        defaultEQ = "10*pow(x, 0.45)" #Default equation if no regional curve equation entered
        expression = defaultEQ.replace("x", DAstr)
        defaultEquationXSW = eval(expression)
        return defaultEquationXSW
    else:
        expression = ratio.replace("x", DAstr)
        equationXSW = eval(expression)
        return equationXSW

##------------------------------------------------------------------------------------------------------------------
## splitline function splits the streamline into segments for use in creating transects
## This function is from the 'Perpendicular Transects' script by Mateus Ferreira
## function is derived from this source: http://nodedangles.wordpress.com/2011/05/01/quick-dirty-arcpy-batch-splitting-polylines-to-a-specific-length/
##------------------------------------------------------------------------------------------------------------------
def splitline (inFC,FCName,alongDist):

    OutDir = arcpy.env.workspace
    outFCName = FCName
    outFC = OutDir+"/"+outFCName
    
    def distPoint(p1, p2):
        calc1 = p1.X - p2.X
        calc2 = p1.Y - p2.Y

        return math.sqrt((calc1**2)+(calc2**2))

    def midpoint(prevpoint,nextpoint,targetDist,totalDist):
        newX = prevpoint.X + ((nextpoint.X - prevpoint.X) * (targetDist/totalDist))
        newY = prevpoint.Y + ((nextpoint.Y - prevpoint.Y) * (targetDist/totalDist))
        return arcpy.Point(newX, newY)

    def splitShape(feat,splitDist):
        # Count the number of points in the current multipart feature
        #
        partcount = feat.partCount
        partnum = 0
        # Enter while loop for each part in the feature (if a singlepart feature
        # this will occur only once)
        #
        lineArray = arcpy.Array()

        while partnum < partcount:
              # Print the part number
              #
              #print "Part " + str(partnum) + ":"
              part = feat.getPart(partnum)
              #print part.count

              totalDist = 0

              pnt = part.next()
              pntcount = 0

              prevpoint = None
              shapelist = []

              # Enter while loop for each vertex
              #
              while pnt:

                    if not (prevpoint is None):
                        thisDist = distPoint(prevpoint,pnt)
                        maxAdditionalDist = splitDist - totalDist

                        #print thisDist, totalDist, maxAdditionalDist

                        if (totalDist+thisDist)> splitDist:
                              while(totalDist+thisDist) > splitDist:
                                    maxAdditionalDist = splitDist - totalDist
                                    #print thisDist, totalDist, maxAdditionalDist
                                    newpoint = midpoint(prevpoint,pnt,maxAdditionalDist,thisDist)
                                    lineArray.add(newpoint)
                                    shapelist.append(lineArray)

                                    lineArray = arcpy.Array()
                                    lineArray.add(newpoint)
                                    prevpoint = newpoint
                                    thisDist = distPoint(prevpoint,pnt)
                                    totalDist = 0

                              lineArray.add(pnt)
                              totalDist+=thisDist
                        else:
                              totalDist+=thisDist
                              lineArray.add(pnt)
                              #shapelist.append(lineArray)
                    else:
                        lineArray.add(pnt)
                        totalDist = 0

                    prevpoint = pnt                
                    pntcount += 1

                    pnt = part.next()

                    # If pnt is null, either the part is finished or there is an
                    #   interior ring
                    #
                    if not pnt:
                        pnt = part.next()
                        if pnt:
                              print "Interior Ring:"
              partnum += 1

        if (lineArray.count > 1):
              shapelist.append(lineArray)

        return shapelist

    if arcpy.Exists(outFC):
        arcpy.Delete_management(outFC)

    arcpy.Copy_management(inFC,outFC)

    #origDesc = arcpy.Describe(inFC)
    #sR = origDesc.spatialReference

    #revDesc = arcpy.Describe(outFC)
    #revDesc.ShapeFieldName

    deleterows = arcpy.UpdateCursor(outFC)
    for iDRow in deleterows:       
         deleterows.deleteRow(iDRow)

    try:
        del iDRow
        del deleterows
    except:
        pass

    inputRows = arcpy.SearchCursor(inFC)
    outputRows = arcpy.InsertCursor(outFC)
    fields = arcpy.ListFields(inFC)

    numRecords = int(arcpy.GetCount_management(inFC).getOutput(0))
    OnePercentThreshold = numRecords // 100

    #printit(numRecords)

    iCounter = 0
    iCounter2 = 0

    for iInRow in inputRows:
        inGeom = iInRow.shape
        iCounter+=1
        iCounter2+=1    
        if (iCounter2 > (OnePercentThreshold+0)):
              #printit("Processing Record "+str(iCounter) + " of "+ str(numRecords))
              iCounter2=0

        if (inGeom.length > alongDist):
              shapeList = splitShape(iInRow.shape,alongDist)

              for itmp in shapeList:
                    newRow = outputRows.newRow()
                    for ifield in fields:
                        if (ifield.editable):
                              newRow.setValue(ifield.name,iInRow.getValue(ifield.name))
                    newRow.shape = itmp
                    outputRows.insertRow(newRow)
        else:
              outputRows.insertRow(iInRow)

    del inputRows
    del outputRows

##------------------------------------------------------------------------------------------------------------------
## GetAzimuthPolyline function; this function is part of the 'Perpendicular Transects' script by Mateus Ferreira
##------------------------------------------------------------------------------------------------------------------
def GetAzimuthPolyline(shape):
    if shape.lastPoint.Y - shape.firstPoint.Y == 0:
        pass
        #When I would enter anything less than 57 (meters while testing) for the 'Distance between transects'
        #I was getting a divide by 0 error; couldn't identify the exact issue.  I believe it has
        #something to do with the transects being non perpendicular to the streamline.
        #I observed a few cross-sections that had the middle 'stream' portion that ran parallel to the streamline.
    else:
        radian = math.atan((shape.lastPoint.X - shape.firstPoint.X)/(shape.lastPoint.Y - shape.firstPoint.Y))
        degrees = radian * 180 / math.pi
        return degrees

##------------------------------------------------------------------------------------------------------------------
## Azimuth function; this function is part of the 'Perpendicular Transects' script by Mateus Ferreira
##------------------------------------------------------------------------------------------------------------------
def Azimuth(direction):
    if direction < 0:
        azimuth = direction + 360
        return azimuth
    else:
        return direction

##------------------------------------------------------------------------------------------------------------------
## findNulls function; this function is part of the 'Perpendicular Transects' script by Mateus Ferreira
##------------------------------------------------------------------------------------------------------------------
def findNulls(fieldValue):
    if fieldValue is None:
        return 0
    elif fieldValue is not None:
        return fieldValue

##------------------------------------------------------------------------------------------------------------------
## Azline1 function; this function is part of the 'Perpendicular Transects' script by Mateus Ferreira
##------------------------------------------------------------------------------------------------------------------
def Azline1(azimuth):
    az1 = azimuth + 90
    if az1 > 360:
        az1-=360
        return az1
    else:
        return az1

##------------------------------------------------------------------------------------------------------------------
## Azline2 function; this function is part of the 'Perpendicular Transects' script by Mateus Ferreira
##------------------------------------------------------------------------------------------------------------------
def Azline2(azimuth):
    az2 = azimuth - 90
    if az2 < 0:
        az2+=360
        return az2
    else:
        return az2

##------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------
###Input data types
##wrkSpace = Workspace
##strmCenterline = Feature Layer
##reachField = Field (obtained from strmCenterline)
##DistanceSplit = Double
##ratio = String
##valleyLeft = Feature Layer
##valleyRight = Feature Layer
##inFlowAcc = Raster Layer
##CrossSections = Feature Class [output]

#Set user-defined inputs
wrkSpace = arcpy.GetParameterAsText(0) #Set working directory (folder) #r"D:\Users\miharris\Desktop\NemadjiTest\ScriptTestFolder"
strmCenterline = arcpy.GetParameterAsText(1) #Input streamline (that has route features already) #r"D:\Users\miharris\Desktop\NemadjiTest\ScriptTestFolder\testStream_meters.shp" 
reachField = arcpy.GetParameterAsText(2) #Select the field that defines the route feature #"Route"
DistanceSplit = float(arcpy.GetParameterAsText(3)) #How often along the streamline you want transects (units are FEET) #100
ratio = arcpy.GetParameterAsText(4).lower() #Regional curve equation for bankfull width #"2.9761*pow(x, 0.9233)"
valleyLeft = arcpy.GetParameterAsText(5) #Input valley wall left #r"D:\Users\miharris\Desktop\NemadjiTest\ScriptTestFolder\testValley_Left.shp"
valleyRight = arcpy.GetParameterAsText(6) #Input valley wall right #r"D:\Users\miharris\Desktop\NemadjiTest\ScriptTestFolder\testValley_Right.shp"
inFlowAcc = arcpy.GetParameterAsText(7) #Input the flow accumulation file (.tif) that was created by Marks tools #r"D:\Users\miharris\Desktop\NemadjiTest\RasterClip\flow_acc.tif"
CrossSections = arcpy.GetParameterAsText(8) #Output cross-sections #r"D:\Users\miharris\Desktop\NemadjiTest\ScriptTestFolder\OutputTransect.shp"

#Set local variables
endPoints = "in_memory\\endPoints"
endptTransctJoin = "in_memory\\endptTransctJoin"
tableLayer = "in_memory\\tableLayer"
smthValleyLeft = "in_memory\\smthValleyLeft"
smthValleyRight = "in_memory\\smthValleyRight"
lineRiverLeft = "in_memory\\lineRiverLeft"
lineRiverRight = "in_memory\\lineRiverRight"
valleyLineMerge = "in_memory\\valleyLineMerge"
valleyTransectJoin = "in_memory\\valleyTransectJoin"

#Set workspace to the user input
arcpy.env.workspace = wrkSpace

#Find "inFlowAcc" cell size
ref2 = arcpy.Describe(inFlowAcc).spatialReference
units2 = ref2.linearUnitName
cellDesc = arcpy.Describe(inFlowAcc)
if units2 == "Meter":
    cellSize = cellDesc.meanCellWidth / 0.3048
else:
    arcpy.AddMessage("linear unit of " + cellDesc.basename + " not 'Meter'...assuming 'Feet'")
    cellSize = cellDesc.meanCellWidth #assumes if linear unit is not meters then it is feet

#Find "strmCenterline" linear units for 'locating features along routes'
spatial_reference = arcpy.Describe(strmCenterline).spatialReference
units3 = spatial_reference.linearUnitName
cellDesc = arcpy.Describe(strmCenterline)
#Set "DistanceSplit" linear units for 'locating features along routes' and change distance units to FEET
if units3 == "Meter":
    locateFeatUnit = " Meters"
    DistanceSplitFT = DistanceSplit * 0.3048
elif units3 == "Feet":
    locateFeatUnit = " Feet"
    DistanceSplitFT = DistanceSplit
elif units3 == "Foot_US":
    locateFeatUnit = " Feet"
    DistanceSplitFT = DistanceSplit
else:
    arcpy.AddMessage("Linear unit for 'Distance between transects' is not feet OR meters 'Locate Features Along Route' will use 'meters' as the distance for finding endpoints... ")
    locateFeatUnit = " Meters"
    DistanceSplitFT = DistanceSplit * 0.3048

try:
    #Setting a geodatabase workspace for the script
    WorkFolder = arcpy.env.workspace
    General_GDB = WorkFolder + "\General.gdb"
    arcpy.CreateFileGDB_management(WorkFolder, "General", "CURRENT")
    arcpy.env.workspace = General_GDB
    
   #Checkout Spatial Analyst extension
    arcpy.AddMessage("Checking license... ")
    if arcpy.CheckExtension("Spatial") == "Available":
        arcpy.CheckOutExtension("Spatial")
        arcpy.AddMessage("Spatial Analyst license checked out... ")
    else:
        arcpy.AddMessage("Spatial Analyst license needed... ")
        raise LicenseError

    #Create progressor for script completion percentage labels in ArcMap for the number of streamline records
    recordCount = int(arcpy.GetCount_management(strmCenterline).getOutput(0))
    arcpy.SetProgressor("step", "Creating Stream ID's...", 0, recordCount, 1)

    #Add field to name and reference cross-section features
    arcpy.AddField_management(strmCenterline, "strm_ID", "LONG")
    objID = arcpy.Describe(strmCenterline).OIDFieldName
    
    #This updateCursor loop is different than the others in this script because it turns out that having the objID
    #var results in an error 'object does not support indexing' the objID var is there to make sure I can obtain
    #either the OID or FID if the streamline is either a shapefile or feature class in a geodatabase
    updateSort = arcpy.UpdateCursor(strmCenterline, "", "", "", objID)
    ID = 1
    for row in updateSort:
        arcpy.SetProgressorLabel("Labeling Stream ID " + str(ID) + " of " + str(recordCount) + "...")
        row.setValue("strm_ID", ID)
        updateSort.updateRow(row)
        ID = ID + 1
        arcpy.SetProgressorPosition()
    del updateSort
    del row
    
    #Add fields for calculation
    arcpy.env.outputMFlag = "Disabled"
    arcpy.env.outputZFlag = "Disabled"
    arcpy.AddMessage("Obtaining maximum drainage area for stream ID... ")
    #Calculate maximum flowAcc value for which each strm_ID passes over
    zonalStats = arcpy.sa.ZonalStatisticsAsTable(strmCenterline, "strm_ID", inFlowAcc, "zonalStatsAsTable", "DATA", "MAXIMUM")
    
    #Join based on attributes
    arcpy.JoinField_management(strmCenterline, "strm_ID", zonalStats, "strm_ID", "MAX")
    
    #Add fields for calculation
    FieldsNames=["DA_sqmi", "W_RegCrve", "Width_125"]
    for fn in FieldsNames:
        arcpy.AddField_management (strmCenterline, fn, "DOUBLE")
    
    #Calculate drainage area, bankfull width, and 25% >BKFw fields
    arcpy.AddMessage("Calculating drainage area and bankfull width... ")
    updateRows = arcpy.da.UpdateCursor(strmCenterline, ["DA_sqmi", "W_RegCrve", "Width_125", "MAX"])
    #fields; index[0] = "DA_sqmi", index[1] = "W_RegCrve", index[2] = "Width_125", index[3] = "MAX"
    iter1 = 1
    for row in updateRows:
        arcpy.SetProgressorLabel("Finding DA and Wbkf on feature " + str(iter1) + " of " + str(recordCount) + "...")
        row[0] = drainageSQMI(row[3])
        row[1] = thresholdXSW(row[0])
        row[2] = row[1] * 1.25 #25% greater than the bankfull width, this will be the width of the transect
        iter1 = iter1 + 1
        arcpy.SetProgressorPosition()
        updateRows.updateRow(row)
    del updateRows
    del row
    
    ##------------------------------------------------------------------------------------------------------------------
    ##------------------------------------------------------------------------------------------------------------------
    #Transect Tool
    
    #Unsplit Line
    arcpy.AddMessage("Creating transects across stream channel... ")
    LineDissolve="LineDissolve"
    arcpy.Dissolve_management (strmCenterline, LineDissolve, "Width_125", "", "SINGLE_PART") #Harris added "Width_125 instead of hard-coded number
    LineSplit="LineSplit" 
    
    #Run SplitLine function
    splitline(LineDissolve, LineSplit, DistanceSplitFT)
    
    #Add fields to LineSplit
    FieldsNames=["LineID", "Direction", "Azimuth", "X_mid", "Y_mid", "AziLine_1", "AziLine_2", "Distance"]
    for fn in FieldsNames:
        arcpy.AddField_management (LineSplit, fn, "DOUBLE")
    
    #Create progressor for script completion percentage labels in ArcMap [Record Count for LineSplit]
    recordCount1 = int(arcpy.GetCount_management(LineSplit).getOutput(0))

    #Calculate fields; Harris updated original code [Mateus'] to use da.cursor modules because they are more efficient
    updateRows = arcpy.da.UpdateCursor(LineSplit, ["SHAPE@", "LineID", "OBJECTID", "Direction", "Azimuth",
                                                   "X_mid", "Y_mid", "AziLine_1", "AziLine_2", "Width_125", "Distance"])
    iter2 = 1
    for row in updateRows:
        arcpy.SetProgressorLabel("Calculating fields for record " + str(iter2) + " of " + str(recordCount1) + "...")
        row[1] = row[2]
        row[3] = GetAzimuthPolyline(row[0])
        row[3] = findNulls(row[3])
        row[4] = Azimuth(row[3])
        row[5] = row[0].positionAlongLine(0.50, True).firstPoint.X
        row[6] = row[0].positionAlongLine(0.50, True).firstPoint.Y
        row[7] = Azline1(row[4])
        row[8] = Azline2(row[4])
        #Set distance to the length of the transect divided by 2, the transect scripts creates the length on either side of the streamline
        #So if the input width it 30 ft, it will create 30 ft on each side (=60 ft), we only want the total to be 30 ft
        row[10] = row[9] / 2
        iter2 = iter2 + 1
        arcpy.SetProgressorPosition()
        updateRows.updateRow(row)
    del updateRows
    del row
    
    #Generate Azline1 and Azline2
    Azline1="Azline1"
    Azline2="Azline2"
    #Hard-coded as FEET below, I don't think I need to change that.
    arcpy.BearingDistanceToLine_management (LineSplit, Azline1, "X_mid", "Y_mid", "Distance", "FEET", "AziLine_1", "DEGREES", "GEODESIC", "LineID", spatial_reference)
    arcpy.BearingDistanceToLine_management (LineSplit, Azline2, "X_mid", "Y_mid", "Distance", "FEET", "AziLine_2", "DEGREES", "GEODESIC", "LineID", spatial_reference)
    
    #Create Azline and append Azline1 and Azline2
    Azline="Azline"
    arcpy.CreateFeatureclass_management(General_GDB, "Azline", "POLYLINE", "", "", "", spatial_reference)
    arcpy.AddField_management (Azline, "LineID", "DOUBLE")
    arcpy.Append_management ([Azline1, Azline2], Azline, "NO_TEST")
    
    #Dissolve Azline
    Azline_Dissolve="Azline_Dissolve"
    arcpy.Dissolve_management (Azline, Azline_Dissolve,"LineID", "", "SINGLE_PART")
    
    #Add Fields to Azline_Dissolve
    FieldsNames2=["x_start", "y_start", "x_end", "y_end"]
    for fn2 in FieldsNames2:
        arcpy.AddField_management (Azline_Dissolve, fn2, "DOUBLE")
    
    #Calculate Azline_Dissolve fields
    updateRows = arcpy.da.UpdateCursor(Azline_Dissolve, ["SHAPE@", "x_start", "y_start", "x_end", "y_end"])
    iter3 = 1
    for row in updateRows:
        arcpy.SetProgressorLabel("Calculating geometry of record " + str(iter3) + " of " + str(recordCount1) + "...")
        row[1] = row[0].positionAlongLine(0, True).firstPoint.X
        row[2] = row[0].positionAlongLine(0, True).firstPoint.Y
        row[3] = row[0].positionAlongLine(1, True).firstPoint.X
        row[4] = row[0].positionAlongLine(1, True).firstPoint.Y
        iter3 = iter3 + 1
        arcpy.SetProgressorPosition()
        updateRows.updateRow(row)
    del updateRows
    del row
    
    #Generate output file
    OutputTransect = arcpy.XYToLine_management (Azline_Dissolve, "OutputTransect", "x_start", "y_start", "x_end", "y_end", "", "", spatial_reference)
    arcpy.AddMessage("Done creating transects across stream channel... ")
    
    ##------------------------------------------------------------------------------------------------------------------
    ##------------------------------------------------------------------------------------------------------------------
    
    #Copy OutputTransect so it doesn't have 'linesegment' field (used for joining later)
    OutputTransectNoLineSeg = arcpy.CopyFeatures_management(OutputTransect, "OutputTransectNoLineSeg")
    
    #Create progressor for script completion percentage labels in ArcMap [Record Count for OutputTransect]
    recordCount2 = int(arcpy.GetCount_management(OutputTransect).getOutput(0))

    #Create endpoints at polyline ends of transects
    arcpy.AddMessage("Creating endpoints on transects... ")
    arcpy.AddField_management(OutputTransect, "Line_Seg", "LONG")
    updateRows = arcpy.da.UpdateCursor(OutputTransect, ["OID", "Line_Seg"]) #FID for shapefile; OID for feature class in geodatabase
    iter4 = 1
    for row in updateRows:
        arcpy.SetProgressorLabel("Calculating RowID " + str(iter4) + " of " + str(recordCount2) + "...")
        row[1] = row[0] #Set Line_Seg equal to the OID
        iter4 = iter4 + 1
        arcpy.SetProgressorPosition()
        updateRows.updateRow(row)
    del updateRows
    del row
    
    #Polyline EndPoints to Points
    arcpy.FeatureVerticesToPoints_management(OutputTransect, endPoints, "BOTH_ENDS")
    
    #Add XY Coordinates to the endPoints
    fldNmes = ["X_UTM", "Y_UTM"]
    for fl in fldNmes:
        arcpy.AddField_management(endPoints, fl, "DOUBLE")
    
    #Create progressor for script completion percentage labels in ArcMap [Record Count for endPoints]
    recordCount3 = int(arcpy.GetCount_management(endPoints).getOutput(0))

    #Calculate XY Coordinates of the endPoints
    updateRows = arcpy.da.UpdateCursor(endPoints, ["SHAPE@X", "SHAPE@Y", "X_UTM", "Y_UTM"])
    iter5 = 1
    for row in updateRows:
        arcpy.SetProgressorLabel("Calculating XY Coordinates for feature " + str(iter5) + " of " + str(recordCount3) + "...")
        row[2] = row[0] #X = X_UTM
        row[3] = row[1] #Y = Y_UTM
        iter5 = iter5 + 1
        arcpy.SetProgressorPosition()
        updateRows.updateRow(row)
        #del row
    del updateRows
    del row
    
    #Use spatial join on target(endPoints) and join(OutputTransect)
    arcpy.SpatialJoin_analysis(endPoints, OutputTransect, endptTransctJoin)
    
    #Locate endpoint features along streamline route
    #Correct left/right orientation depends on direction that the streamlines were produced;
    #either from Mark's tools, or from digitizing direction
    arcpy.AddMessage("Locating endpoints along streamline route... ")
    props = "RID POINT MEAS"
    #Set search radius for locating features, set to the value of transect spacing minus 1
    distnce = str(DistanceSplit) + locateFeatUnit
    #Save the table as "locate_points" back into the wrkSpace user-inputted save location
    arcpy.LocateFeaturesAlongRoutes_lr(endptTransctJoin, strmCenterline, "Route", distnce, "locate_points", props, "FIRST", "DISTANCE", "NO_ZERO", "FIELDS", "NO_M_DIRECTION")
    
    #Create XY event layer from linear referencing table
    arcpy.MakeXYEventLayer_management("locate_points", "X_UTM", "Y_UTM", tableLayer, spatial_reference)
    
    #Find endpoints that are either on the left/right side of the streamline route
    arcpy.AddMessage("Finding endpoints that are river left vs. river right... ")
    rivLeftExpr = '"DISTANCE" <= 0'
    rivRightExpr = '"DISTANCE" >= 0'
    ptRiverLeft = arcpy.FeatureClassToFeatureClass_conversion(tableLayer, General_GDB, "ptRiverLeft", rivLeftExpr)
    ptRiverRight = arcpy.FeatureClassToFeatureClass_conversion(tableLayer, General_GDB, "ptRiverRight", rivRightExpr) 

    #tableLayer was causing lock when trying to delete General_GDB (Arcgis error 000603)
    arcpy.Delete_management(tableLayer)
    
    #Smooth the valley walls
    arcpy.SmoothLine_cartography(valleyLeft, smthValleyLeft, "BEZIER_INTERPOLATION") 
    arcpy.SmoothLine_cartography(valleyRight, smthValleyRight, "BEZIER_INTERPOLATION")
    
    #Find the nearest valley wall from the transect endpoints
    arcpy.Near_analysis(ptRiverLeft, smthValleyLeft, "", "LOCATION", "NO_ANGLE")
    arcpy.Near_analysis(ptRiverRight, smthValleyRight, "", "LOCATION", "NO_ANGLE")
    
    #Create lines from endpoints to valley wall
    arcpy.AddMessage("Creating lines from endpoints to valley wall... ")
    arcpy.XYToLine_management(ptRiverLeft, lineRiverLeft, "X_UTM", "Y_UTM", "NEAR_X", "NEAR_Y", "GEODESIC", "", spatial_reference)
    arcpy.XYToLine_management(ptRiverRight, lineRiverRight, "NEAR_X", "NEAR_Y", "X_UTM", "Y_UTM", "GEODESIC", "", spatial_reference)
    
    #Merge valley lines into one feature class
    arcpy.Merge_management([lineRiverRight, lineRiverLeft, OutputTransectNoLineSeg], valleyLineMerge)
    
    #Join valley lines to stream transect lines
    arcpy.SpatialJoin_analysis(valleyLineMerge, endptTransctJoin, valleyTransectJoin)
    
    #Dissolve by Line_Seg to create completed cross-sections
    arcpy.Dissolve_management(valleyTransectJoin, CrossSections, "Line_Seg")
    xsCount = int(arcpy.GetCount_management(CrossSections).getOutput(0))
    arcpy.AddMessage(str(xsCount) + " cross-sections created!... ")

    #Delete General.gdb
    arcpy.Delete_management(General_GDB)

except arcpy.ExecuteError:
    msgs = arcpy.GetMessages(2)
    arcpy.AddError(msgs)

except:
    #Copied from Mark's streamtype_v2.py script
    #Get the traceback object
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]
    #Concatenate error information into message string
    pymsg = 'PYTHON ERRORS:\nTraceback info:\n{0}\nError Info:\n{1}'\
          .format(tbinfo, str(sys.exc_info()[1]))
    msgs = 'ArcPy ERRORS:\n {0}\n'.format(arcpy.GetMessages(2))
    #Return python error messages
    arcpy.AddError(pymsg)
    arcpy.AddError(msgs)

#Delete items in memory
arcpy.Delete_management("in_memory\\endPoints")
arcpy.Delete_management("in_memory\\endptTransctJoin")
arcpy.Delete_management("in_memory\\tableLayer")
arcpy.Delete_management("in_memory\\smthValleyLeft")
arcpy.Delete_management("in_memory\\smthValleyRight")
arcpy.Delete_management("in_memory\\lineRiverLeft")
arcpy.Delete_management("in_memory\\lineRiverRight")
arcpy.Delete_management("in_memory\\valleyLineMerge")
arcpy.Delete_management("in_memory\\valleyTransectJoin")

print "Done"