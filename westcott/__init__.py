from .westcott import *
from .log_handlers import *
from .tables import *
from .user import *
from .westcott_gfactors import *

__version__='1.0.0'
__author__='David A. Matters and Aaron M. Hurst'

_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    """Method to return absolute path of the data files inside the root of 
    the Python package."""
    return os.path.join(_ROOT, path)

