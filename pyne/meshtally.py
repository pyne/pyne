#!/usr/bin/env python

"""

Module for loading an SCDMesh from an MCNP meshtal file

"""

import sys
import linecache
import itertools

from itaps import iBase, iMesh
from collections import namedtuple, Iterable
from r2s.scdmesh import ScdMesh, ScdMeshError

class meshtally:   

    def __init__(self, filename = None, line_count = 1):
        self.line_count = line_count
        self.Number = -1
        self.Type = 'NA'
        self.x_bounds = []
        self.y_bounds = []
        self.z_bounds = []
        self.e_bounds = []
        self.col_idx = {}
        self.Flux = []
        self.RelE = []
        self.e_bins = -1
        self.spatialPoints = -1
        self.sm = -1

        if filename == None:
            return    
        else:
            self.read_meshtally_head(filename, self.line_count)
            self.read_boundaries(filename, self.line_count)
            self.read_column_order(filename, self.line_count)
            self.read_values(filename, self.line_count)
            self.create_mesh()
            self.tag_fluxes()
   
    def read_meshtally_head(self, filename, line_count = -1):
        if line_count == -1:
            line_count = self.line_count

        flag = 2
        while flag:
            Line = linecache.getline(filename, line_count)
            line_count = line_count+1

            # empty line
            if (Line.split() == []):
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
            if (self.Number == -1) & ('mesh' in Line) & \
                                     ('tally' in Line) & ('number' in Line):
                for i in Line:
                    if i.replace('.','',1).isdigit():
                        self.Number = i
                        flag = flag-1
                        break
            # get meshtally type
            elif (self.Type == 'NA') & ('this' in Line) & \
                                       ('is' in Line) & ('mesh' in Line):
                if ('and' in Line) & ('neutron' in Line) & ('photon' in Line):
                    self.Type = 'np'
                    flag = flag-1
                elif ('neutron' in Line):
                    self.Type = 'n'
                    flag = flag-1
                elif ('photon' in Line):
                    self.Type = 'p'
                    flag = flag-1
            elif (line_count-self.line_count) > 50:
                return False

        self.line_count =  line_count
        return True

    def read_boundaries(self, filename, line_count = -1):
        if line_count == -1:
            line_count = self.line_count
     
        flag = 4              
        while flag:
            Line = linecache.getline( filename , line_count )
            line_count = line_count + 1      
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
            if (self.x_bounds == []) & ('x' in Line) & ('direction' in Line):
                for i in Line:
                    if i.replace('-','').replace('.','').isdigit():
                        self.x_bounds.append(i)
                flag = flag-1
            # get y bounds
            elif (self.y_bounds == []) & ('y' in Line) & ('direction' in Line):
                for i in Line:
                    if i.replace('-','').replace('.','').isdigit():
                        self.y_bounds.append(i)
                flag = flag-1
            # get z bounds
            elif (self.z_bounds == []) & ('z' in Line) & ('direction' in Line):
                for i in Line:
                    if i.replace('-','').replace('.','').isdigit():
                        self.z_bounds.append(i)
                flag = flag-1
            # get energy bounds
            elif (self.e_bounds == []) & ('energy' in Line) & ('bin' in Line):
                for i in Line:
                    if i.replace('+','').replace('-','').replace('.','').replace('e','').isdigit():
                        self.e_bounds.append(i.replace('e','E')) 
                flag = flag-1

            elif (line_count-self.line_count) > 50:
                return False
 
        self.line_count = line_count
        return True
    

    def read_column_order( self, filename, line_count = -1 ):
        if line_count == -1:
            line_count = self.line_count
        flag = 1
        while flag:
            Line = temp = linecache.getline( filename, line_count )
            line_count = line_count+1            
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
                colnames = temp.lower().replace('rel ','rel').replace('rslt * ','rslt').strip().split()
                self.col_idx = dict(zip(colnames,range(0,len(colnames))))
                flag = flag-1

            if (line_count - self.line_count) > 50:
                return False

        self.line_count = line_count
        return True
    
    def read_values(self, filename, line_count = -1):  
        if line_count == -1:
            line_count = self.line_count
        while True:
            Line = linecache.getline(filename, line_count)
            line_count = line_count + 1	    
            if (Line.split() == []):
                break
            else:
                Line = Line.split() 	        
                try:
                    self.Flux.append(Line[self.col_idx['result']])
                except IndexError:
                    print 'flux not found'
                try:
                    self.RelE.append(Line[self.col_idx['relerror']])
                except IndexError:
                    print 'RelE not found'
        self.line_count = line_count-1
        return True

    def get_values_line(self, filename, line_count = -1):
        if line_count == -1:
            line_count = self.line_count-1
        i = 50
        Line = linecache.getline(filename, line_count).lower().split()
        while i:
            if ('x' in Line) & ('y' in Line) & ('z' in Line) & ('result' in Line):
                return line_count+1
            i = i-1

    def create_mesh(self):
        #Calculating pertinent information from meshtal header and input
        self.spatialPoints = (len(self.x_bounds)-1) * \
                             (len(self.y_bounds)-1) * (len(self.z_bounds)-1)        
        if len(self.e_bounds) > 2 :
            self.e_bins = len(self.e_bounds) 
            #don't substract 1; cancels with totals bin
        elif len(self.e_bounds) == 2: 
            #for 1 energy bin, meshtal doesn't have TOTALS group
            self.e_bins = 1 
        self.sm = ScdMesh(self.x_bounds, self.y_bounds, self.z_bounds)
        return True

    def tag_fluxes(self, Norm = 1.0):
        voxels = list(self.sm.iterateHex('xyz'))
        for E_Group in range(1, self.e_bins + 1): 
            # Create tags if they do not already exist
            if self.e_bins == 1 or E_Group != self.e_bins: 
                # tag name for each E_Bin
                flux_str = '{0}_group_{1:03d}'.format( self.Type, E_Group)
                error_str = '{0}_group_{1:03d}_error'.format( self.Type, E_Group)
            elif E_Group == self.e_bins : 
                # tag name for totals group
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
                flux_data.append( float( self.Flux[ 
                    point + (E_Group-1)* self.spatialPoints ]) * Norm)
                error_data.append( float( self.RelE[ 
                    point + (E_Group-1)* self.spatialPoints ]))
            #Tag data for energy group 'E_Group' onto all voxels
            tag_flux[voxels] = flux_data
            tag_error[voxels] = error_data


