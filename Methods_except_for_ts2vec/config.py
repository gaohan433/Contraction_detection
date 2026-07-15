
class Config:
    PRE_ = 'E:/code/iup/github'
    ROOT_DIR = PRE_

    ICELAND_DIR ='D:/research/data/iceland'
    ICELAND_RESOURCE = PRE_ + '/resource'
    
  
    def __init__(self):
        return
    

# my_project/config.py
from yacs.config import CfgNode as CN


_C = CN()


_C.PARAM = CN()


_C.PARAM.CHANNEL = [9,10]

_C.PARAM.FS = 6
_C.PARAM.FILTER = [0.1,3]
_C.PARAM.FILTER_TYPE = 0
_C.PARAM.NORMALIZE = False


_C.DETECT = CN()
_C.DETECT.METHOD_NAME = 'sample_entropy'

_C.DETECT.WINDOW_SIZE = 60
_C.DETECT.SHIFT = 20
_C.DETECT.BASE_WINDOW_SIZE = 180
_C.DETECT.BASE_SHIFT = 20
_C.DETECT.MEAN_FILTER_SIZE = 10





def get_cfg_defaults():
  """Get a yacs CfgNode object with default values for my_project."""
  # Return a clone so that the defaults will not be altered
  # This is for the "local variable" use pattern
  return _C.clone()


        