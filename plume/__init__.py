from ._version import get_versions
from .plume import Plume

__version__ = get_versions()["version"]
del get_versions
