# File: redistribute_vf_lines.py
# Name: Matthew V Bilskie, PhD

#----------------------------------------------------------
# M O D U L E S
#----------------------------------------------------------
import numpy as np
import PyAdcirc
import shapefile
from shapely.geometry import Point, LineString
#----------------------------------------------------------

#----------------------------------------------------------
# F U N C T I O N   C U T
#----------------------------------------------------------
#
# Cuts a line in two at a distance from its starting point.
# result = function(line, distance)
#----------------------------------------------------------
def cut(line, distance):
    # Cuts a line in two at a distance from its starting point
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        # line.project() returns the distance along the geometric object to a point nearest the other object
        pd = line.project(Point(p))
        #print pd
        if pd == distance:
            return [
                    LineString(coords[:i+1]),
                    LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance)
            return [
                    LineString(coords[:i] + [(cp.x, cp.y)]),
                    LineString([(cp.x, cp.y)] + coords[i:])]

    # Example on using the cut function
    #line = LineString(pts)
    #print([list(x.coords) for x in cut(line, 1000.0)])
    #out = [list(x.coords) for x in cut(line, 1000.0)]

#----------------------------------------------------------
# F U N C T I O N   R E D I S T  _ V E R T I C I E S
#----------------------------------------------------------
#
# Redistributes points along a line.
# result = function(geometry, distance)
#----------------------------------------------------------
def redist_verticies(geom, distance):
    if geom.geom_type == 'LineString':
        num_vert = int(round(geom.length / distance))
        if num_vert == 0:
            num_vert = 1
        return LineString(
                [geom.interpolate(float(n) / num_vert, normalized=True)
                    for n in range(num_vert + 1)])
    elif geom.geom_type == 'MultiLineString':
        parts = [redist_vertices(part,distance)
                for part in geom]
        return type(geom)([p for p in parts if not p.is_empty])
    else:
        raise ValueError('unhandeled geometry %s', (geom.geom_type,))

#----------------------------------------------------------
# F U N C T I O N   I N T E R P V A L       
#----------------------------------------------------------
#
# Function to obtain z-value from a mesh for a given coordinate pair.
# result = function(mesh,x,y)
#----------------------------------------------------------
def interpval(mesh, x, y):
    # Find element point p is in
    elem = mesh.findElement(x,y)
    # Calculate the value at point p based on the barycentric weights
    try:
        weight = mesh.element(elem).interpolationWeights(x,y)
        z1 = weight[0] * mesh.size[mesh.element(elem).node(0).id()]
        z2 = weight[1] * mesh.size[mesh.element(elem).node(1).id()]
        z3 = weight[2] * mesh.size[mesh.element(elem).node(2).id()]
        return z1 + z2 + z3
    except:
        nnode = mesh.findNearestNode(x,y)
        return mesh.node(nnode).z()

#----------------------------------------------------------
# F U N C T I O N   R E D I S T R I B U T E
#----------------------------------------------------------
#
# Function to redistribute points based on the local mesh
# size function.
# result = function(line, mesh)
#----------------------------------------------------------
def redistribute(line, mesh):
    # Create a list of coordinates along the line
    coords = list(line.coords)

    numNewPts = 0   # Number of verticies on the redistributed line
    lastp = 0       # Index of the last vertex on the original line used to add a point to the rd line

    rd = None

    for p in range(len(coords)):

        # Get the size function value at the first p index to start
        if numNewPts < 1:
            zlength = interpval(mesh, coords[0][0], coords[0][1])

        # Check the distance from the last vertex to point p
        distchk = Point(coords[lastp]).distance(Point(coords[p]))

        # Check the distance exceeds the desired size function, zlength
        if distchk >= zlength:
            if numNewPts < 1:
                # Add the first line segment to the redistributed line
                rd = LineString([coords[lastp], coords[p]])
                numNewPts = numNewPts + 1
                lastp = p
            else:
                # Update tcoords with the lastest rd coordinates
                tcoords = list(rd.coords)
                rd = LineString(tcoords[:] + [coords[p]])
                numNewPts = numNewPts + 1
                lastp = p
        
        # Add the last vertex
        if (p == len(coords)-1) and (rd is not None):
            tcoords = list(rd.coords)
            rd = LineString(tcoords[:] + [coords[p]])
            numNewPts = numNewPts + 1
            lastp = p
        
        # Get the size function value at the last p found 
        zlength = interpval(mesh, coords[lastp][0], coords[lastp][1])

    return rd

#----------------------------------------------------------
# F U N C T I O N   C H E C K R E D I S T
#----------------------------------------------------------
#
# Function to handle the last line segment after redistribution.
# result = function(oline,rdline)
#----------------------------------------------------------
def checkredist(origline, rdline, mesh):

    # Create a tuple of line coordinates
    rdcoords = list(rdline.coords)

    # Number of verticies/points along the line
    numrdpts = len(rdcoords)

    # Calculate the length of the last line segment for the redistributed line
    remainder = Point(rdcoords[numrdpts-2]).distance(Point(rdcoords[numrdpts-1]))

    # Grab the size function that should be used for the last line segment
    sizefunc = interpval(mesh,rdcoords[numrdpts-1][0],rdcoords[numrdpts-1][1])

    # lthres is the lower threshold of the remainder, which at its smallest
    # can be the vertex spacing of the original line multiplied by the number
    # of redistributed segements. This is assuming the original line has
    # constant vertex spacing.
    lthres = (origline.length/(len(list(origline.coords))-1)) * numrdpts-1
    
    maxIter = 5  
    iters = 0
    maxnotreached = True
    while (maxnotreached) and (remainder > lthres) and (iters < maxIter):
    #while (maxnotreached) and (remainder > lthres):
        rdline, remainder, maxnotreached = splitremainder(origline, rdline, remainder, mesh)
        iters = iters + 1

    ratio = min(remainder,sizefunc) / max(remainder,sizefunc)
    if (ratio < 0.2):
        # Add remainder to the last line segment
        rdcoords = list(rdline.coords)
        numrdpts = len(rdcoords)
        rdline = LineString(rdcoords[:numrdpts-2] + [rdcoords[numrdpts-1]])
        return rdline
    #elif (ratio > 0.80):
    else:
        # Cut our losses and trim off the remainder / last line segment
        if (rdline is not None):
            rdcoords = list(rdline.coords)
            numrdpts = len(rdcoords)
            rdline = LineString(rdcoords[:numrdpts-1])
            return rdline
        else:
            return None

#----------------------------------------------------------
# F U N C T I O N   S P L I T R E M A I N D E R
#----------------------------------------------------------
#
# Function to split up the remainder of a line segment
# among the the other line segments in a normalized fashion.
# result, result2, result3 = function(origline, rdline, remainder, mesh)
#----------------------------------------------------------
def splitremainder(origline, rdline, remainder, mesh):

    maxnotreached = True
    rd = None

    # Create a list of coordinates along the lines
    origcoords = list(origline.coords)
    rdcoords = list(rdline.coords)

    numrdpts = len(rdcoords)

    # Compute the weight of each line segment
    weight = []
    #for p in range(len(rdcoords)-1):
    for p in range(len(rdcoords)-2): # -2 b/c we don't need a weight for the last line segement for wish we wish to fix
        weight.append( Point(rdcoords[p]).distance(Point(rdcoords[p+1])) / (rdline.length-remainder) )

    numNewPts = 0   # Number of verticies on the redistributed line
    lsegment = 0    # Keep track of redistributed line segments
    lastp = 0       # Index of the last vertex on the original line used to add a point to the rd line
    tsum = 0
    for p in range(len(origcoords)):

        if lsegment < len(weight):
            # Find the line segment length, which should be near identical to
            # the desired size function on the first iteration
            slength = Point(rdcoords[lsegment]).distance(Point(rdcoords[lsegment+1]))
            sadd = weight[lsegment]*remainder
            
            # Determine the size function at point p
            zlength = interpval(mesh, origcoords[p][0], origcoords[p][1])
            
            # Check to see if we are increasing the line segment length by too much
            # which is defined as 10% of the original size function
            if (sadd / zlength) <= 0.10:
                zlength = zlength + sadd
            else:
                # Still increase, but do not increase beyond 10%
                zlength = zlength + 0.10*zlength
                maxnotreached = False

        # Check the distance from the last vertex to point p
        distchk = Point(origcoords[lastp]).distance(Point(origcoords[p]))

        # Check if the distance exceeds the desired size function, zlength
        if distchk >= zlength:
            #print(distchk,zlength)
            #print(interpval(mesh,origcoords[lastp][0],origcoords[lastp][1]),zlength,weight[lsegment]*remainder)
            if numNewPts < 1:
                # Add the first line segment to the redistributed line
                rd = LineString([origcoords[lastp], origcoords[p]])
                numNewPts = numNewPts + 1
                lastp = p
                lsegment = lsegment + 1
            else:
                # Update tcoords with the lastest rd coordinates
                tcoords = list(rd.coords)
                rd = LineString(tcoords[:] + [origcoords[p]])
                numNewPts = numNewPts + 1
                lastp = p
                lsegment = lsegment + 1

    # Add the last vertex
    if (rd is not None):
        tcoords = list(rd.coords)
        rd = LineString(tcoords[:] + [origcoords[p]])
        numNewPts = numNewPts + 1
        lastp = p
        lsegment = lsegment + 1

        # Final update of the coordinates of the redistributed line
        tcoords = list(rd.coords)
        # Update the remainder for the new redistributed line
        remainder = Point(tcoords[numNewPts-1]).distance(Point(tcoords[numNewPts]))

    return rd, remainder, maxnotreached

#----------------------------------------------------------
# F U N C T I O N   L I N E I N M E S H
#----------------------------------------------------------
#
# Function to check if the entire line segment is in the mesh.
# result = function(line, mesh)
#----------------------------------------------------------
def lineinmesh(line, mesh):

    coords = list(line.coords)

    for p in range(len(coords)):
        try:
            z = interpval(mesh, coords[p][0], coords[p][1])
        except:
            return False

    return True

#----------------------------------------------------------
# M A I N
#----------------------------------------------------------
#
#----------------------------------------------------------
def main():

    # Initialize Variables
    remthreshold = 0.10
    rdl = None
    
    # Load mesh
    #inmesh = raw_input('Name of ADCIRC mesh: ')
    inmesh = 'NGOM_SizeFunction_utm16.grd'
    print('Reading in mesh file and computing size function...')
    mesh = PyAdcirc.Mesh(inmesh)
    ierr = mesh.read()
    if ierr == 0:
        exit(ierr)
    mesh.size = mesh.computeMeshSize()

    # Load shapefile
    #insf = raw_input('Name of shapefile (w/out ext.): ')
    insf = 'PascBig5m_Combined_Ridge_forRedist'
    print('Reading in initial shapefile data...')
    sf = shapefile.Reader(insf)
    shapes = sf.shapes()
    numShapes = len(shapes) # Number of shapes

    #outsf = raw_input('Name of redistributed output shapefile (w/out ext.): ')
    outsf = 'PascBig5m_Combined_Ridge_forRedist_Redistributed'

    # Create shapefile for writing
    print('Creating output shapefile for writing...')
    w = shapefile.Writer(outsf, shapeType = 3)
    w.field('name', 'C')

    for i in range(numShapes):
        
        print('Working on line ',i,' out of ',numShapes)

        s = sf.shape(i) # Read shape i
        line = LineString(s.points) #s.points returns a list of tuples containing (x,y) coordinates of each point in the shape

        # Make sure all points in line are within the meshing domain
        # This should be fixed to see if only a small portion of the lines are out to
        # just remove the line segments that are outside the domain.
        #isinmesh = lineinmesh(line, mesh)
        #if isinmesh == False:
        #    continue
        
        # Redistribute all lines to 1 m spacing
        rdl = redist_verticies(line, 1.0)

        # Re-distribute verticies according to the size function
        rdl = redistribute(rdl, mesh)

        # Now I need a check to examine the last line segment distance
        if (rdl is not None):
            rdl = checkredist(line, rdl, mesh)
       
            if (rdl is not None):
                # Write the redistributed line to the shapefile
                rdl = [rdl.coords[:]]
                w.line(rdl)
                w.record('linestring1')

    # Close shapefile 
    w.close()

    return
#
#----------------------------------------------------------
if __name__=="__main__":
    main()

