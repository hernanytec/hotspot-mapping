from ._base_estimator import BaseEstimator
from ..tessellation import GridTessellationTiler
from ._kernel_functions import get_kernel_function, VALID_KERNELS
from ..distance import Distance

import numpy as np
import geopandas as gpd

from sklearn.metrics import pairwise_distances

from shapely.geometry import Point


class KernelDensity(BaseEstimator):

    def __init__(self, bandwidth=0.1, kernel='quartic', method=None):

        """KernelDensityEstimation

        Parameters
        ----------
        bandwidth : float, default 0.1
            The bandwidth of the kernel.

        kernel : str, default is 'quartic'
            The kernel function to use. 
            Valid kernels are ['quartic','linear','exponential']

        method : string {None, 'kdtree', 'balltree'}
            Which method to use for finding nearest edge to each point.
            If None, we manually find each edge one at a time. If 'kdtree' we use
            scipy.spatial.cKDTree for very fast euclidean search. Recommended for
            projected graphs. If 'balltree', we use sklearn.neighbors.BallTree for
            fast haversine search. Recommended for unprojected graphs.

        """
        if bandwidth <= 0:
            raise ValueError("bandwidth must be positive")
        if kernel not in VALID_KERNELS:
            raise ValueError("invalid kernel: '{0}'".format(kernel))

        self.bandwidth = bandwidth
        self.kernel_function = get_kernel_function(kernel)

        self._score = None
        self.method = method
        
    def fit(self, tessellation):
        """
        Parameters
        ----------
        tessellation : GridTessellationTiler 
            Object that represent the study area.
        """
        self.tessellation = tessellation
        
        self.distance = Distance(self.tessellation, method=self.method)
        
        return self

    
    def score(self,  events):
        """Compute kde score for each geometry in tessellation

        Parameters
        ----------
        events : GeodataFrame or numpy array
            geodataframe or geoseries of shapely Points
            list or numpy array of x and y point coords
         
        Returns
        -------
        _score_value : ndarray
            score array with the same size of tessellation, score of all events aggregates
        """
        self.n = len(events)
        
        events = check_events_input(events)

        # for each event get a list which contains indexes of edges inside bandwidth
        indexes, distances = self.distance.nearest_indexes_within_radius(events, radius=self.bandwidth, return_distances=True)
        
        kde_scores = np.zeros(len(self.tessellation.geodataframe))

        for idx, dist in zip(indexes, distances):
            if len(idx) > 0:
                # apply kernel function to get kde scores
                scores = self.kernel_function(dist, self.bandwidth)
                
                kde_scores[np.array(idx).reshape(1, -1)] += scores

        self._score_value = kde_scores * (1/(self.n*self.bandwidth))
        
        return self._score_value


def check_events_input(events):
    """Check and validate input of events
    
    Parameters
    ----------
    events : GeodataFrame or numpy array
       geodataframe or geoseries of shapely Points
       list or numpy array of x and y point coords

    Returns
    -------
    array containing x and coordinates of events
    """
    if isinstance(events, (gpd.GeoDataFrame, gpd.GeoSeries)):
        events = np.array([events.geometry.x.values, events.geometry.y.values]).T
    elif isinstance(events, (list, np.ndarray)):
        try:
            events = np.array(events).reshape(-1, 2)
        except:
            raise ValueError('invalid shape of array')
            
    return events