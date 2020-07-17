import numpy as np
from defs import *
from values import *

def NaKCl2 (cell,memb_id,xNaKCl2,delmu):

	dJNaKCl2 = cell.area[memb_id[0],memb_id[1]]*xNaKCl2*(delmu[0,memb_id[0],memb_id[1]]+delmu[1,memb_id[0],memb_id[1]]+2*delmu[2,memb_id[0],memb_id[1]])

	return [0,1,2],[dJNaKCl2,dJNaKCl2,2*dJNaKCl2]

