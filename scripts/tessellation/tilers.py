import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon

class TessellationTilers:

    def __init__(self):
        pass
  
class GridTessellationTiler(TessellationTilers):
    """Creates an object with evenly spaced grid values and the shapely polygons for each grid cell from a cell size and a shapely polygon.
    
    Parameters
    ----------
    cell_size : int
        The size of a cell in the grid in meters.
    
    geometry : Shapely geometry
        Shapely geometry to be used to make the grid on top of it.
    
    Returns
    -------
    tessellation : tesselation object with centroids values and shapely geometry attributes.
        
    Examples
    --------
    >>> from pyhotspot.tesselation import GridTesselationTiler
    >>> tesselation = GridTesselationTiler(100, geometry)
    >>> tesselation.centroids
    array([[[-124680.72022751, 9578735.86374085],...,[-123680.72022751, 9578735.86374085]]])
    >>> tesselation.polygons
    geometry
    0	POLYGON ((-124180.720 9579235.864, -124180.720...)) ...
    """

    
    def __init__(self, cell_size, geometry):
        super().__init__()
        
        self._geodataframe = None
        self._geodataframe_generated = False
        
        self.cell_size = cell_size
        self.geometry = geometry
        
        # Get the limits of the boundbox
        min_x, min_y, max_x, max_y = geometry.bounds
        # Expand the boundaries to fit the proportion
        self.min_x, self.max_x = self._fit_bound(min_x, max_x, self.cell_size)
        self.min_y, self.max_y = self._fit_bound(min_y, max_y, self.cell_size)
        
        # Take the number of cells of each dimension
        self.num_cells_x = round((self.max_x - self.min_x) / self.cell_size)
        self.num_cells_y = round((self.max_y - self.min_y) / self.cell_size)
        # Create the empty grid matrix
        self.centroids = np.zeros((self.num_cells_y, self.num_cells_x, 2))
        # Take the indices of the grid matrix
        indices = np.indices((self.num_cells_y, self.num_cells_x))
        row_indices = indices[0]
        column_indices = indices[1]
        # Initialize all values of the grid (centered)
        self.centroids[:,:,0] = self.min_x + (cell_size / 2) + (column_indices * self.cell_size)
        self.centroids[:,:,1] = self.min_y + (cell_size / 2) + (row_indices * self.cell_size)
    
    # Expand the boundaries to fit the proportion
    def _fit_bound(self, minimum, maximum, proportion):
        diff_rest = (maximum - minimum) % proportion
        if diff_rest != 0.0:
            diff_to_rest = proportion - diff_rest
            minimum -= diff_to_rest / 2.0
            maximum += diff_to_rest / 2.0

        return minimum, maximum

    def _generate_geodataframe(self):

    	# Calculate the size of a half cell to expand
        half_cell = self.cell_size/2
        
        centroids_reshape = self.centroids.reshape(self.centroids.shape[0]*self.centroids.shape[1], 2)
        
        # Construct the polygons of the hotspots cells
        self.cell_poly = [Polygon(((x+half_cell, y+half_cell),
                     (x+half_cell, y-half_cell),
                     (x-half_cell, y-half_cell),
                     (x-half_cell, y+half_cell))) for x, y in centroids_reshape]


        self._geodataframe_generated = True
        # create the GeoDataFrame only with the cell polygons
        self._geodataframe = gpd.GeoDataFrame(geometry=self.cell_poly)
    
    @property
    def geodataframe(self):
    	if not self._geodataframe_generated:
    		self._generate_geodataframe()
    	return self._geodataframe
