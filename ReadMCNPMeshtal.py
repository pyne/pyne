#############################

import sys
import linecache
from itaps import iBase, iMesh

import itertools
from collections import namedtuple, Iterable
from r2s.scdmesh import ScdMesh, ScdMeshError

#############################
    
def find_in_iterable( x, iterable):
    for i, item in enumerate(iterable):
        if item == x:
            return i
    return -1

#############################

class Meshtally:
    
    Number = -1
    Type = 'NA'
    
    sm = -1 ###iMesh.Mesh()   ##########
    
    x_bounds = []
    y_bounds = []
    z_bounds = []
    e_bounds = []
    e_bins = -1
    spatialPoints = -1

    flux = []
    RelE = []

#   a = ['x','y','z','result','rel_error','energy','volume','res*vol']
    a = [-1,-1,-1,-1,-1,-1,-1,-1]

    linecount = 0

    def __init__( self, filename = None, LineCount = 0 ):
        if filename == None:
            pass
        elif LineCount < 0:
            pass
        else:
            self.linecount = LineCount
            self.ReadMeshtallyHead( filename )
            self.ReadBoundaries( filename )
            self.ReadFirstLine( filename )
            self.ReadValues( filename )
            self.createMesh()
            self.tagFluxes()
            
    def __repr__( self ):
        return "Meshtally Number:%s Meshtally Type:%s" % (self.Number, self.Type)

    def getX( self ):
        return self.x_bounds
    def getY( self ):
        return self.y_bounds
    def getZ( self ):
        return self.z_bounds
    def getE( self ):
        return self.e_bins
    def getflux( self ):
        return self.flux
    def getRelE( self ):
        return self.RelE
    def getA( self ):
        return self.a
    def getLineCount( self ):
        return self.linecount

#############################
    
    def ReadMeshtallyHead( self, FileName, LineCount = -1 ):

        #
        #  Action:
        #      Extract Meshtally number and type
        #  Input:
        #      FileName: path of meshtal file
        #      LineCount: Line from which parsing should start
        #  Output:
        #      Returns True after extraction is complete
        #

        if LineCount == -1:
            LineCount = self.linecount

        flag = 2
    
        while flag:
            Line = linecache.getline( FileName, LineCount )
            LineCount = LineCount + 1

            # empty line
            if ( Line.split() == []):
                continue
        
            # ignore comments
            x = Line.strip().find('#')
            if ( x == 0 ):
                continue
            elif (x > 0):
                Line = str(Line.split('#')[:1]).lower().split()
            else:
                Line = Line.lower().split()
 
            # get meshtally number
            if ( self.Number == -1 ) & ('mesh' in Line) & ('tally' in Line) & ('number' in Line):
                for i in Line:
                    if i.replace('.','',1).isdigit():
                        self.Number = i
                        flag = flag-1
                        break
    
            # get meshtally type
            elif ( self.Type == 'NA' ) & ('is' in Line) & ('mesh' in Line) & (('tally.' in Line) | ('tally' in Line)):
                if ('and' in Line ) & ('neutron' in Line ) & ('photon' in Line ):
                    self.Type = 'np'
                    flag = flag-1
                elif ('neutron' in Line):
                    self.Type = 'n'
                    flag = flag-1
                elif ('photon' in Line):
                    self.Type = 'p'
                    flag = flag-1
        
            elif ( LineCount - self.linecount ) > 50:
                return False

        self.linecount =  LineCount
        return True
                
#############################

    def ReadBoundaries( self, FileName, LineCount = -1 ):

        """
        #
        # Action:
        #     Extract x,y,z and energy bin bounds
        # Input:
        #     FileName: path of meshtal file
        #     LineCount: Line from which parsing should start
        # Output:
        #     x_direction: array of x_direction boundary points
        #     y_direction: array of y_direction boundary points
        #     z_direction: array of z_direction boundary points
        #     e_bounds: array of energy bounds
        #     LineCount: Last line used for parsing
        #
        """

        if LineCount == -1:
            LineCount = self.linecount

        flag = 4  
            
        while flag:
            Line = linecache.getline( FileName , LineCount )
            LineCount = LineCount + 1
        
            # empty line
            if ( Line.split() == []):
                continue
        
            # ignore comments
            x = Line.strip().find('#')
            if ( x == 0 ):
                continue
            elif (x > 0):
                Line = str(Line.split('#')[:1]).lower().replace(':','').split()
            else:
                Line = Line.lower().replace(':','').split()

            # get x bounds
            if ( self.x_bounds == []) & ('x' in Line) & ('direction' in Line):
                for i in Line:
                    if i.replace('-','').replace('.','').isdigit():
                        self.x_bounds.append(i)
                flag = flag-1
                    
            # get y bounds
            elif ( self.y_bounds == []) & ('y' in Line) & ('direction' in Line):
                for i in Line:
                    if i.replace('-','').replace('.','').isdigit():
                        self.y_bounds.append(i)
                flag = flag-1
                    
            # get z bounds
            elif ( self.z_bounds == []) & ('z' in Line) & ('direction' in Line):
                for i in Line:
                    if i.replace('-','').replace('.','').isdigit():
                        self.z_bounds.append(i)
                flag = flag-1

            # get energy bounds
            elif ( self.e_bounds == []) & ('energy' in Line) & ('bin' in Line):
                for i in Line:
                    if i.replace('+','').replace('-','').replace('.','').replace('e','').isdigit():
                        self.e_bounds.append(i.replace('e','E')) 
                flag = flag-1               

            elif ( LineCount - self.linecount ) > 50:
                return False
            
        self.linecount = LineCount
        return True
    
#############################

    def ReadFirstLine( self, FileName, LineCount = -1 ):

        """
        # Action:
        #     Extract order in which result and error values are listed
        # Input:
        #     FileName: path of meshtal file
        #     LineCount: Line from which parsing should start
        # Output:
        #     a: array containing order of listing of results
        #     LineCount: Last line used for parsing
        """
   
        if LineCount == -1:
            LineCount = self.linecount

        flag = 1
    
        while flag:
            Line = temp = linecache.getline( FileName, LineCount )
            LineCount = LineCount+1
            
            # empty line
            if ( Line.split() == []):
                continue
        
            # ignore comments
            x = Line.strip().find('#')
            if ( x == 0 ):
                continue
            elif (x > 0):
                Line = str(Line.split('#')[:1]).lower().replace(':','').split()
            else:
                Line = Line.lower().replace(':','').split()

            if ('x' in Line) & ('y' in Line) & ('z' in Line) & ('result' in Line):
                temp = temp.lower().replace('rel ','rel').replace('rslt * ','rslt').strip().split()
            
                self.a[0] = find_in_iterable('x',temp)
                self.a[1] = find_in_iterable('y',temp)
                self.a[2] = find_in_iterable('z',temp)
                self.a[3] = find_in_iterable('result',temp)
                self.a[4] = find_in_iterable('relerror',temp)
                self.a[5] = find_in_iterable('energy',temp)
                self.a[6] = find_in_iterable('volume',temp)
                self.a[7] = find_in_iterable('rsltvol',temp)
                
                flag = flag-1

            if ( LineCount - self.linecount ) > 50:
                return False

        self.linecount = LineCount
        return True
    
#############################

    def ReadValues( self, FileName, LineCount = -1 ):
	
        """
        # Action:
        #     Extract meshtally values
        # Input:
        #     FileName: path of meshtal file
        #     LineCount: Line from which parsing should start
        #     a: array containing order of listing of results
        # Output:
        #     x_val, y_val, z_val, n_val: arrays containing x, y, z and energy values
        #     r_val, e_val: arrays containing result and relerror values
        #     LineCount: Last line used for parsing
        """
        
        if LineCount == -1:
            LineCount = self.linecount

        flag = 1

        while flag==1:
    
            Line = linecache.getline(FileName, LineCount)
            LineCount = LineCount + 1
	    
            if ( Line.split() == [] ):
                if ( linecache.getline(FileName, LineCount+1).split() == [] ):
                    flag = flag - 1
		    break

            if (Line.lower().find('number') > 0):
		break

            else:
                Line = Line.split()
 	        
                """             	
                try:
                    x_val.append(Line[self.a[0]])         # x_val
                except IndexError:
                    print 'x_val not found'
                            
                try:
                    y_val.append(Line[a[1]])         # y_val
                except IndexError:
                    print 'y_val not found'
                            
                try:
                    z_val.append(Line[a[2]])         # z_val
                except IndexError:
                    print 'z_val not found'
                """           
                """
                try:
                    e_val.append(Line[a[5]])         # energy
                except IndexError:
                    print 'e_val not found'
                """

                try:
                    self.flux.append(Line[self.a[3]])         # result
                except IndexError:
                    print 'flux not found'

                try:
                    self.RelE.append(Line[self.a[4]])         # relerror
                except IndexError:
                    print 'RelE not found'

        self.linecount = LineCount
        return True

#############################

    def createMesh( self ):

        """
        #
        # Action:
        #     Create a structured mesh using the ScdMesh class
        #
        """

    #   Calculating pertinent information from meshtal header and input
        self.spatialPoints = ( len ( self.x_bounds )-1 )*(len( self.y_bounds )-1)*(len ( self.z_bounds ) -1)
        
        if len( self.e_bounds ) > 2 :
            self.e_bins = len( self.e_bounds ) #don't substract 1; cancels with totals bin
        elif len( self.e_bounds ) == 2 : #for 1 energy bin, meshtal doesn't have TOTALS group
            self.e_bins = 1 
    
        self.sm = ScdMesh( self.x_bounds, self.y_bounds, self.z_bounds )
        print self.x_bounds
        print self.y_bounds
        print self.z_bounds
        print self.sm     
        
        return True

#############################

    def tagFluxes( self, Norm = 1.0 ):
    
        voxels = list(self.sm.iterateHex('xyz'))

        for E_Group in range(1, self.e_bins + 1): 
            # Create tags if they do not already exist
            if self.e_bins == 1 or E_Group != self.e_bins: # tag name for each E_Bin
                flux_str = '{0}_group_{1:03d}'.format( self.Type, E_Group)
                error_str = '{0}_group_{1:03d}_error'.format( self.Type, E_Group)
            elif E_Group == self.e_bins : # tag name for totals group
                flux_str = self.Type + '_group_total'
                error_str = self.Type + '_group_total_error'
    
            try:
                tag_flux = self.sm.imesh.createTag( flux_str, 1, float)
            except iBase.TagAlreadyExistsError:
                tag_flux = self.sm.imesh.getTagHandle( flux_str)
            try:
                tag_error = self.sm.imesh.createTag( error_str , 1, float)
            except iBase.TagAlreadyExistsError:
                tag_error = self.sm.imesh.getTagHandle( error_str)

            #Create lists of data from meshtal file for energy group 'E_Group'
            flux_data = []
            error_data = []
            for point in range( 0, self.spatialPoints) :
                flux_data.append( float( self.flux[ point + (E_Group-1)* self.spatialPoints ]) * Norm)
                error_data.append( float( self.RelE[ point + (E_Group-1)* self.spatialPoints ]))

            #Tag data for energy group 'E_Group' onto all voxels
            tag_flux[voxels] = flux_data
            tag_error[voxels] = error_data


#############################

class MCNPmeshtal:
    
    SM = list()
    Date = '-1'
    Time = '-1'
    Version = -1
    NHistories = -1

    linecount = 0

    def __init__( self, FileName = None, LineCount = 0):
        if FileName == None:
            pass
        else:
            self.ReadMeshtalHead( FileName, LineCount )
            self.ReadTallies

#############################

    def ReadMeshtalHead( FileName, LineCount = 0 ):

        """
        #
        # Action:
        #     Extract header data
        # Input:
        #     FileName: path of meshtal file
        #     LineCount: Line from which parsing should start
        # Output:
        #     mcnp version and no. of histories
        #
        """

        flag = 2 
    
        while flag:
            Line = linecache.getline( FileName, LineCount )
            LineCount = LineCount+1
        
            # empty line
            if ( Line.split() == []):
                continue
        
            # ignore comments
            x = Line.strip().find('#')
            if ( x == 0 ):
                continue
            elif (x > 0):
                Line = str(Line.split('#')[:1]).lower().split()    
            else:
                Line = Line.lower().split()
        
            # get mcnp version
            if ( Version == -1 ) & ( 'mcnp' in Line ) & ( 'version' in Line ):
                for i in range(Line.index('version')+1,len(Line)):            
                    if Line[i].replace('.','').isdigit():
                        Version = Line[i]
                        Flag = Flag-1
                        break
            
            # get number of histories
            elif ( NHistories == -1 ) & ( 'number' in Line ) & ( 'histories' in Line ):
                for i in Line:
                    if i.replace('.','').isdigit():
                        NHistories = i
                        Flag = Flag-1
                        break
         
            elif (LineCount - self.linecount) > 50:
                return False

        self.linecount = LineCount
        return True
                
#############################

    def ReadTallies( self, FileName, LineCount ):
        
        readNextTally = 1
        linecount = LineCount

        while ( readNextTally ):
            meshtally = Meshtally( FileName, linecount )
            self.SM.append(meshtally)
            linecount = meshtally.getLineCount()
            
            Line = linecache.getline( FileName, linecount )

            if Line.split() == []:
                Line = linecache.getline( FileName, linecount + 1)

            if Line.lower.find('mesh'):
                linecount = linecount + 1
            else:
                readNextTally = readNextTally - 1
 
        return True

################################

m = Meshtally('test1.txt')
print "created mesh m"

################################



