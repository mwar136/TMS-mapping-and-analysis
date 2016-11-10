#Author: Matthew Ward mwar136@aucklanduni.ac.nz

import sys
import os
import logging
import time
import numpy as np

from IPython.config.loader import PyFileConfigLoader
from argparse import ArgumentParser

from PyDragonfly import Dragonfly_Module, MT_EXIT, CMessage, copy_to_msg, copy_from_msg
import Dragonfly_config as rc
from dragonfly_utils import respond_to_ping

from pyface.timer.api import Timer
from traits.api import HasTraits, Instance
from traitsui.api import *
from mayavi.core.ui.api import SceneEditor, MlabSceneModel
from tvtk.api import tvtk
from mayavi import mlab, tools, modules
import mayavi
import threading
import Queue

queue = Queue.Queue()

class PlotHead(threading.Thread):

    def __init__(self, parent, config_file, server):#, parent):
        #HasTraits.__init__(self)
        threading.Thread.__init__(self)
        self.daemon = True
        self.count = 0
        self.parent = parent
        self.plot_vertex_vec = np.array([3,-2,2])
        self.load_config(config_file)
        self.setup_dragonfly(server)
        self.start()
    
    def load_config(self, config_file):
        cfg = PyFileConfigLoader(config_file)
        cfg.load_config()
        self.config = cfg.config
        self.filename = self.config.head_model
        
    def process_message(self, msg):
        # read a Dragonfly message
        msg_type = msg.GetHeader().msg_type
        dest_mod_id = msg.GetHeader().dest_mod_id
        if  msg_type == MT_EXIT:
            if (dest_mod_id == 0) or (dest_mod_id == self.mod.GetModuleID()):
                print 'Received MT_EXIT, disconnecting...'
                self.mod.SendSignal(rc.MT_EXIT_ACK)
                self.mod.DisconnectFromMMM()
                return
        elif msg_type == rc.MT_PING:
            respond_to_ping(self.mod, msg, 'PlotHead')
        elif msg_type == rc.MT_PLOT_POSITION:
            in_mdf = rc.MDF_PLOT_POSITION()
            copy_from_msg(in_mdf, msg)
            tail = np.array(in_mdf.xyz[:])*0.127 + (self.plot_vertex_vec)#Hotspot position
            head = np.array(in_mdf.ori[:3])/4 #Vector head of coil, used to find ori
                 
            if np.any(np.isnan(tail)) == True:
                pass
            elif np.any(np.isnan(head)) == True:
                 pass
            elif np.any(np.isinf(tail)) == True:
                pass
            elif np.any(np.isinf(head)) == True:
                pass
            else:
                queue.put(np.vstack((head, tail)))
                self.count=+1
                print 'sent message'
        elif msg_type == rc.MT_MNOME_STATE:
            in_mdf = rc.MDF_MNOME_STATE()
            copy_from_msg(in_mdf, msg)
            if in_mdf.state == 0:
                print 'got clear'
                self.parent.reset = True
               
                
        
    def setup_dragonfly(self, server):
        subscriptions = [MT_EXIT, \
                         rc.MT_PING, \
                         rc.MT_PLOT_POSITION, \
                         rc.MT_MNOME_STATE]
        self.mod = Dragonfly_Module(0, 0)
        self.mod.ConnectToMMM(server)
        for sub in subscriptions:
            self.mod.Subscribe(sub)
        self.mod.SendModuleReady()
        print "Connected to Dragonfly at ", server
    
   # def timer_event(self, parent):
    def run(self):
        while True:
            msg = CMessage()
            rcv = self.mod.ReadMessage(msg, 0)
            if rcv == 1:
                self.process_message(msg)
             
                
class MainWindow(HasTraits):
    scene = Instance(MlabSceneModel, ())

    view = View(Item('scene', editor=SceneEditor(), resizable=True,
                    show_label=False),
                    resizable=True, height=500.00, width=750.0)

    def __init__(self, config, server):
        HasTraits.__init__(self)
        self.mayavi_view = PlotHead(self, config, server)
        self.tail_data = np.zeros((2,3))
        self.head_data = np.zeros((2,3))
        self.reset = False
        self.init_plot()
        self.timer = Timer(2000, self.update_plot)
        self.timer.Start()
        
    def init_plot(self):
        # create a window with 14 plots (7 rows x 2 columns)
        ## create a window with 8 plots (4 rows x 2 columns)

        reader = tvtk.OBJReader()
        reader.file_name = self.mayavi_view.filename
        mapper = tvtk.PolyDataMapper()
        mapper.input = reader.output
        p = tvtk.Property(opacity=0.99)
        actor = tvtk.Actor(property=p)
        mapper.color_mode = 0x000000
        actor.mapper = mapper
        actor.orientation = (180,0,90)
        self.scene.add_actor(actor)
        mlab.points3d(11.5, -1.3, -14)
        mlab.points3d(-6.6, -1.3, -13.5)
        mlab.points3d(3,-2,2)
        mlab.points3d(2, 9.5, -9.5)
        self.arrows = mlab.quiver3d(self.tail_data[:,0], self.tail_data[:,1], self.tail_data[:,2],
                        self.head_data[:,0], self.head_data[:,1], self.head_data[:,2], scale_mode='vector', scale_factor=1.0)
        self.dots = mlab.points3d(self.tail_data[:,0], self.tail_data[:,1], self.tail_data[:,2],
        color=(1,1,0),opacity=0.5,scale_mode='vector', scale_factor=1.0)
        self.ar = self.arrows.mlab_source
        self.dot_source = self.dots.mlab_source
        '''
        x_arrows = np.ones((300, 6))
        y_arrows = np.ones((300, 6))
        z_arrows = np.ones((300, 6))
        arrows = np.ones((300, 6))
        for i, pos in enumerate(arrows):
            x_arrows[i, 0] = pos[0] + i
            y_arrows[i, 1] = pos[1] + i
            z_arrows[i, 2] = pos[2] + i
        arrows = np.vstack((x_arrows, y_arrows, z_arrows, np.zeros((1,6))))
        self.tail_data = np.vstack((self.tail_data, arrows[:,:3]))
        mlab.points3d(x_arrows[:,0], x_arrows[:,1], x_arrows[:,2], scale_factor = 0.5)
        mlab.quiver3d(x_arrows[:,0], x_arrows[:,1], x_arrows[:,2],
                       x_arrows[:,3], x_arrows[:,4], x_arrows[:,5], scale_factor = 0.5)
          #  z=arrows[:,2]
        
        self.head_data = np.vstack((self.head_data, arrows[:, -3:]))
        self.ar.reset(x=arrows[:,0] +  self.mayavi_view.plot_vertex_vec[0], y=arrows[:,1] + self.mayavi_view.plot_vertex_vec[1], 
           z=arrows[:,2] + self.mayavi_view.plot_vertex_vec[2], u=arrows[:,3]+10, v=arrows[:,4]+10,w=arrows[:,5]+10)
        '''
         
    def update_plot(self):
        print 'started updating'
        if self.reset:
            while not queue.empty():
                queue.get()
            self.head_data = np.zeros((2,3))
            self.tail_data = np.zeros((2,3))
            self.reset = False
        else:
            while not queue.empty():
                new_values = queue.get()
                print new_values.shape
                self.tail_data = np.vstack((self.tail_data, new_values[1]))
                print 'updated position'
                self.head_data = np.vstack((self.head_data, new_values[0]))
                print 'updated rotation'
        print self.head_data
        print self.tail_data
        self.ar.reset(x=self.tail_data[:,0], y=self.tail_data[:,1], z=self.tail_data[:,2],
                        u=self.head_data[:,0], v=self.head_data[:,1], w=self.head_data[:,2], scale_mode='vector', scale_factor=1.0)
        self.dot_source.reset(x=self.tail_data[:,0], y=self.tail_data[:,1], z=self.tail_data[:,2],
        color=(1,1,0),opacity=0.5,scale_mode='vector', scale_factor=1.0)
        print 'finished'
    
        
        
if __name__ == "__main__":
    parser = ArgumentParser(description = 'Real-time display of TMS Coil position')
    parser.add_argument(type=str, dest='config')
    parser.add_argument(type=str, dest='mm_ip', nargs='?', default='127.0.0.1:7111')
    args = parser.parse_args()
    print("Using config file=%s, MM IP=%s" % (args.config, args.mm_ip))
    frame = MainWindow(args.config, args.mm_ip)
    frame.configure_traits()
    #t = threading.Thread(target =  frame.configure_traits())
    #t.daemon=True
    #t.start()