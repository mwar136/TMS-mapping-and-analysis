#Author: Matthew Ward mwar136@aucklanduni.ac.nz

import numpy as np
import quaternionarray as qa

def quaternion_mean(q):
    average = np.array([0,0,0,0])
    for i in q:
        if check_same(q,i) == False:
            i = change_sign(i)
        average = average + i
    return qa.norm(average/q.shape[0])

def check_same(array, single_quaternion):
    dot = np.dot(array[0], single_quaternion)
    if dot < 0.0:
        return False
    else:
        return True

def change_sign(single_quaternion):
    return single_quaternion*-1