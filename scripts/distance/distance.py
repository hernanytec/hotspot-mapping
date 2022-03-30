from  ..tessellation import *

import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import Point
from scipy.spatial import ckdtree
from sklearn.neighbors import BallTree
from sklearn.metrics import pairwise_distances

class Distance():
    
    def __init__(self, tessellation, method=None):
        """Distance class for fastest neighborhood search

        Parameters
        ----------
        tessellation : GridTessellationTiler
            Object that represent the study area.
        
        method : string {None, 'kdtree', 'balltree'}
            Which method to use for finding nearest edge to each point.
            If None, we manually find each edge one at a time. If 'kdtree' we use
            scipy.spatial.cKDTree for very fast euclidean search. Recommended for
            projected graphs. If 'balltree', we use sklearn.neighbors.BallTree for
            fast haversine search. Recommended for unprojected graphs.

        """
        self.tessellation = tessellation
        self.method = method

        if isinstance(tessellation, GridTessellationTiler):
            centroids = tessellation.geodataframe.centroid

            self.tessellation_points = np.array([centroids.x, centroids.y]).T
            self.tree = ckdtree.cKDTree(self.tessellation_points)

            
    def nearest_index(self, events, return_distance=False):
        """Find the index of closest tessellation geometry to the events

        Parameters
        ----------
        events : GeodataFrame or numpy array
            geodataframe or geoseries of shapely Points
            list or numpy array of x and y point coords
                
        return_distance : bool, default is False
            if true, also return the auclidean distance of events to tessellation geometries 
        
        Returns
        -------
        closest_index : ndarray
            for each point, the index of the geometry that are closest to that point.

        closest_index_distance : ndarray
            if return_distances = True, return also the distances in meters, of closest geometry.
        """
        event_x, event_y = _check_events_input(events)
        points = np.array([event_x, event_y]).T
        
        if isinstance(self.tessellation, GridTessellationTiler):
            self.closest_index = find_grid_cells_by_events(self.tessellation, events)
            
            if return_distance:
                self.closest_index_distance, _ = self.tree.query(points, k=1)
        else:
            if self.method is None:
                edges = self.tessellation.geodataframe.geometry.values
                
                idx_dist = []
                for (px, py) in  points:
                    edge_distances = [(idx, Point(px, py).distance(edge)) for idx, edge  in enumerate(edges)]
                    
                    idx_dist.append(min(edge_distances, key=lambda x: x[1]))
                    
                self.closest_index, self.closest_index_distance = list(zip(*idx_dist))

            elif self.method == 'balltree':
                points_rad = np.deg2rad(points)
                
                dist, idx = self.tree.query(points_rad, k=1, return_distance=True)
                
                self.closest_index = self.extended.loc[idx[:, 0], 'index'].values
                self.closest_index_distance = dist
            elif self.method == 'kdtree':
                dist, idx = self.tree.query(points, k=1)  # Returns ids and distance of closest point
                
                self.closest_index = self.extended.loc[idx, 'index'].values
                self.closest_index_distance = dist
            else:
                raise ValueError("method valid are None, 'kdtree' or 'balltree'")
        
        if return_distance:
            return self.closest_index, self.closest_index_distance
        
        return self.closest_index


    def nearest_indexes_within_radius(self, events, radius, return_distances=False):
        """Find all tessellation geometry indexes all to the events

        Parameters
        ----------
        events : GeodataFrame or numpy array
            geodataframe or geoseries of shapely Points
            list or numpy array of x and y point coords
        
        return_distance : bool, default is False
            if true, also return the auclidean distance of events to all tessellation geometries within a radius. 

        Returns
        -------
        indexes_list : ndarray
            for each point, the index of all geometries within a radius.

        distances : ndarray
            if return_distances = True, return also the distances in meters, of all geometries within radius.
        """
        event_x, event_y = _check_events_input(events)

        points = np.array([event_x, event_y]).T

        if self.method is None:
            #TODO implement for none
            geoms = self.tessellation.geodataframe.geometry.values
            
            idx = []
            dist = []

            for (px, py) in  points:
                geom_distances = [(idx, Point(px, py).distance(geom)) for idx, geom  in enumerate(geoms)]
                
                events_within_radius = list(filter(lambda x: x[1] <= radius, geom_distances)) 
            
                i, d = list(zip(*events_within_radius))

                idx.append(np.array(i))
                dist.append(np.array(d))

            self.closest_index, self.closest_index_distance = idx, dist
            
            if return_distances:
                return self.closest_index, self.closest_index_distance
            
            return self.closest_index
        else:
            # for each point get a list which contains indexes of geometries within a radius
            indexes_list = self.tree.query_ball_point(points, r=radius)
            
            # if only one event is passed
            if isinstance(indexes_list, list):
                indexes_list = [indexes_list]
                points = [points]
                
            if return_distances:
                distances = [pairwise_distances([point], self.tessellation_points[idx])
                            if len(idx) > 0
                            else []
                            for point, idx in zip(points, indexes_list)  ]
                
                return indexes_list, distances 
            
            return indexes_list


def _check_events_input(events):
    """Check and validate input of events
    
    Parameters
    ----------
    events : GeodataFrame or numpy array
       geodataframe or geoseries of shapely Points
       list or numpy array of x and y point coords
    
    Returns
    -------
    tuple containing an array for the x coordinates and another for the y coordinates of the events
    """

    if isinstance(events, Point):
        p = np.array([[events.x, events.y]])
        return p[:, 0], p[:, 1]
    elif isinstance(events, (gpd.GeoDataFrame, gpd.GeoSeries)):
        return (events.geometry.x.values, events.geometry.y.values)
    elif isinstance(events, (list, np.ndarray)):
        events = np.array(events).reshape(-1, 2)
        
        return events[:, 0], events[:, 1]
    else:
        raise ValueError('invalid input')


def grid_coord_to_linear_index(grid_coord, grid_shape):
    """Convert grid coord to a linear index

    Parameters
    ----------
    grid_coord : tuple
        represent de x coord and y coord of grid
    
    grid_shape : tuple
        the number of horizontal an vertical cells of grid 

    Returns
    -------
    return index of this grid cell, if cell is out of grid, return -1
    """
    try:
        return np.ravel_multi_index(grid_coord, grid_shape)
    except ValueError:
        return -1

    
def find_grid_cells_by_events(grid_tiler, events):
    """Find the closest grid cells of events

    Parameters
    ----------
    grid_tiles : GridTessellationTiler
        object that represents the study area

    events : GeodataFrame or numpy array
       geodataframe or geoseries of shapely Points
       list or numpy array of x and y point coords
    """
    event_x, event_y = _check_events_input(events)
    
    grid_coord_x = np.floor((event_x - grid_tiler.min_x) / grid_tiler.cell_size)
    grid_coord_y = np.floor((event_y - grid_tiler.min_y) / grid_tiler.cell_size)
    
    grid_coords = np.array([grid_coord_y, grid_coord_x]).T.astype(int)
    
    grid_shape = (grid_tiler.num_cells_y, grid_tiler.num_cells_x)
    
    # convert all grid coords to linear indexes
    idx = np.apply_along_axis(grid_coord_to_linear_index, axis=1, arr=grid_coords, grid_shape=grid_shape)

    return idx



def get_nearest_points(points, events, k=1, return_dist=False):
    """Query for k-edges closest to a list of events.

    Parameters
    ----------
    points : GeoSeries of Point
        Point objects to approximate the events.

    events : GeoDataFrame
        Geopandas gdf containing events that are within the study area.

    k : int
        the number of closest edges to query.

    return_dist : bool
        if True, return also the distance of the k-closest edges
        Default is False.

    Returns
    -------
    edge_idx : ndarray
        for each point,  the k-indexes of the edges that are closest to that point.

    dist : ndarray
        if return_dist = True, return also the distance in meters, of k-closest points.
    """
    points_array = np.array([points.x, points.y]).T
    events_array = np.array([events.geometry.x, events.geometry.y]).T

    # build a k-d tree for euclidean nearest node search
    kd_tree = ckdtree.cKDTree(points_array)

    # query the tree for nearest node to each point
    events_distance, events_index = kd_tree.query(events_array, k)  # Returns ids of k-closest point

    if return_dist:
        return events_index, events_distance

    return events_index
