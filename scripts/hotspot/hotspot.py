from ..distance import Distance
from ..tessellation import GridTessellationTiler

import geopandas as gpd
import numpy as np
from shapely.geometry import Point

class Hotspot():
    
    def __init__(self, percentage=0.03, crs=None):
        """Create hotspots with percentage of study area
        
        Parameters
        ----------
        percentage : float
            Percentage of total to be used as hotspot.
        """
        self.percentage = percentage
        self.hotspots = gpd.GeoDataFrame()
        self.crs = crs
        

    def fit(self, tessellation, scores):
        """Given a tessellation and its , scores creates a geodataframe containing the hotspots
        
        Parameters
        ----------
        tessellation : GridTessellationTiller
            Representation of study area.
        score : numpy array
            Array with scores for each geometry
        """
        self.tessellation = tessellation

        if(isinstance(tessellation, GridTessellationTiler)):
            self.grid_hotspots(scores)
            
            self.area = self.hotspots.area.sum()
            self.total_area = tessellation.geodataframe.area.sum()
        else:
            print('Should be a tessellation instance')
        
        self.hotspots.crs = self.crs
        return self.hotspots

    def grid_hotspots(self, scores):
        area = self.percentage * self.tessellation.geometry.area
        number_of_cells = int(np.ceil(area / (self.tessellation.cell_size*self.tessellation.cell_size)))

        sort_idx = scores.argsort()[::-1][:number_of_cells]

        self.hotspots = self.tessellation.geodataframe.iloc[sort_idx]


    def contains(self, events, method=None, dist=0.0001):

        distance = Distance(self.tessellation, method=method, dist=dist)
        
        nearest_geom_idx = distance.nearest_index(events)
        geom_in_hotspots = np.isin(nearest_geom_idx, self.hotspots.index.values)
        
        if isinstance(events, Point):
            if geom_in_hotspots[0]:
                return [events]
            return []

        return events[geom_in_hotspots]
