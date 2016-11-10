#Author: Matthew Ward mwar136@aucklanduni.ac.nz
import quaternionarray as qa

class TrackedPoint(object):
    def __init__(self, Qi, Ni, Ti):
        self.Qi = Qi
        self.Xi = Ni - Ti
        
    def get_pos(self, Qk, Tk):
        Qk = qa.norm(Qk)
        Qr = (qa.mult(Qk, qa.inv(self.Qi))).flatten()
        pos = (qa.rotate(Qr, self.Xi)).flatten() + Tk
        return pos, Qr