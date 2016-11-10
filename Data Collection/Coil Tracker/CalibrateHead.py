#Author: Matthew Ward mwar136@aucklanduni.ac.nz

import sys
import os
import logging
import time
import wx
import numpy as np
import quaternionarray as qa
import threading
from wx.lib.pubsub import pub

from IPython.config.loader import PyFileConfigLoader
from argparse import ArgumentParser

from PyDragonfly import Dragonfly_Module, MT_EXIT, CMessage, copy_to_msg, copy_from_msg
import Dragonfly_config as rc
from dragonfly_utils import respond_to_ping

import quaternion_average as mean
from TrackedPoint import TrackedPoint as TP

'''
CalibrateHead is a GUI that allows the generation of a .txt file and saving of anatomical locations of a participant to generate their LCS
'''

class CalibrateHead(wx.Frame):
    
    def __init__(self, config_file, server):
                     
        wx.Frame.__init__(self, None, 1, title="Calibrate Head", size=(450,200))

        czbutton = wx.Button(self, 2, 'Cz', size=(50,30), pos=(10,80))
        ltbutton = wx.Button(self, 3, 'Left Tragus', size=(100,30), pos=(80,80))
        rtbutton = wx.Button(self, 4, 'Right Tragus', size=(100,30), pos=(200, 80))
        nasionbutton = wx.Button(self, 5, 'Nasion', size=(100,30), pos=(320, 80))
        self.text = wx.TextCtrl(self, wx.ID_ANY, '', size=(200, 30), pos=(10, 30))

        dirbutton = wx.Button(self, 6, 'Choose', size=(80,20), pos=(280, 30))
        dirbutton.Bind(wx.EVT_BUTTON, self.choose_directory, id=6)
        
        createbutton = wx.Button(self, 7, 'Create', size=(50,20), pos=(220, 30))
        createbutton.Bind(wx.EVT_BUTTON, self.create_file, id=7)
        
        czbutton.Bind(wx.EVT_BUTTON, self.get_cz, id=2)
        ltbutton.Bind(wx.EVT_BUTTON, self.get_lt, id=3)
        rtbutton.Bind(wx.EVT_BUTTON, self.get_rt, id=4)
        nasionbutton.Bind(wx.EVT_BUTTON, self.get_nasion, id=5)
        
        self.calibrate = Acquirepositions(config_file, server)
        pub.subscribe(self.save_values, 'positions')
        self.Show(True)
        
    def choose_directory(self, event):
        dlg = wx.DirDialog(self, "Choose a directory:",
                           style=wx.DD_DEFAULT_STYLE|wx.DD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            self.location = dlg.GetPath()
        dlg.Destroy()
    
    
    def create_file(self, event):
        self.name = self.text.GetValue()
        textfile = open(self.location + '\\' + self.name + '.txt', 'w')
        textfile.write(self.name)
        textfile.close()
    
    def get_cz(self, event):
        self.calibrate.request('Cz')
        self.collect_positions = True
    
    def get_lt(self, event):
        self.calibrate.request('Left Tragus')
        self.collect_positions = True
    
    def get_rt(self, event):
        self.calibrate.request('Right Tragus')
        self.collect_positions = True
    
    def get_nasion(self, event):
        self.calibrate.request('Nasion')
        self.collect_positions = True
        
    def save_values(self, landmark, landmark_position, glasses_position, glasses_orientation):
        with open(self.location + '\\' + self.name + '.txt', 'a') as f:
            f.write('\n' + landmark)
            f.write('\nGlasses position = ')
            np.savetxt(f, glasses_position[None], fmt='%f' ,delimiter=',')
            f.write('Glasses orientation = ')
            np.savetxt(f, glasses_orientation[None], fmt='%f' ,delimiter=',')
            f.write('Landmark position = ')
            np.savetxt(f, landmark_position[None], fmt='%f' ,delimiter=',')
            f.close()
            sys.stdout.write('\n' + landmark)
        self.calibrated = True

           
'''
Acquirepositions is a background instance called by CalibrateHead to 
receive POLARIS_POSITIION ,messages and perform the necessary calculations of the pointer tool
'''
                        
class Acquirepositions(threading.Thread):

    def __init__(self, config_file, server):
        threading.Thread.__init__(self)
        self.glasses_pos = np.zeros((30,3))
        self.pointer_pos = np.zeros((30,3))
        self.glasses_ori = np.zeros((30,4))
        self.store_glasses = 0
        self.store_pointer = 0
        self.daemon = True
        self.collect_positions = False
        self.load_config(config_file)
        self.setup_dragonfly(server)
        self.start()
     
    def load_config(self, config_file):
        cfg = PyFileConfigLoader(config_file)
        cfg.load_config()
        self.config = cfg.config
        self.glasses = self.config.tools.index('ST568')
        self.pointer = self.config.tools.index('P717')
        self.pointer_Ti = np.array(self.config.tool_list[self.pointer].Ti)
        self.pointer_Qi = qa.norm(np.array(self.config.tool_list[self.pointer].Qi))
        self.pointer_Ni = np.array(self.config.tool_list[self.pointer].Ni)
        self.pointer_tp = TP(self.pointer_Qi, self.pointer_Ni, self.pointer_Ti)
        
        
    def setup_dragonfly(self, server):
        subscriptions = [MT_EXIT, \
                         rc.MT_PING, \
                         rc.MT_POLARIS_POSITION]
        self.mod = Dragonfly_Module(0, 0)
        self.mod.ConnectToMMM(server)
        for sub in subscriptions:
            self.mod.Subscribe(sub)
        self.mod.SendModuleReady()
        print "Connected to Dragonfly at ", server   
        
        
    def run(self):
        while True:
            msg = CMessage()
            rcv = self.mod.ReadMessage(msg, 0.001)
            if rcv == 1:
                msg_type = msg.GetHeader().msg_type
                dest_mod_id = msg.GetHeader().dest_mod_id
                if  msg_type == MT_EXIT:
                    if (dest_mod_id == 0) or (dest_mod_id == self.mod.GetModuleID()):
                        print 'Received MT_EXIT, disconnecting...'
                        self.mod.SendSignal(rc.MT_EXIT_ACK)
                        self.mod.DisconnectFromMMM()
                        break;
                elif msg_type == rc.MT_PING:
                    respond_to_ping(self.mod, msg, 'CalibrateHead')
                else:
                    self.process_message(msg) 
   
    def process_message(self, in_msg):
        # read a Dragonfly message
        msg_type = in_msg.GetHeader().msg_type
        dest_mod_id = in_msg.GetHeader().dest_mod_id
        if self.collect_positions == True:
            if(
               (self.store_pointer >= self.pointer_pos.shape[0]) & 
               (self.store_glasses >= self.glasses_pos.shape[0]) & 
               (self.store_glasses >= self.glasses_ori.shape[0])
            ):
                pub.sendMessage(
                                'positions', landmark=self.landmark, 
                                landmark_position=self.pointer_pos.mean(axis=0),
                                glasses_position=self.glasses_pos.mean(axis=0), 
                                glasses_orientation=mean.quaternion_mean(self.glasses_ori)
                )
                self.collect_positions = False
                self.store_glasses = 0
                self.store_pointer = 0
            
            else:
                if msg_type == rc.MT_POLARIS_POSITION:
                    in_mdf = rc.MDF_POLARIS_POSITION()
                    copy_from_msg(in_mdf, in_msg)
                    positions = np.array(in_mdf.xyz[:])
                    orientations = qa.norm(self.shuffle_q(np.array(in_mdf.ori[:])))
                    if in_mdf.tool_id == (self.glasses + 1):
                        self.glasses_pos[self.store_glasses, :] = positions
                        self.glasses_ori[self.store_glasses, :] = orientations
                        if np.any(np.isnan(positions)) == True:
                            sys.stdout.write('Glasses: nan present')
                            self.store_pointer -= 1
                        elif np.any(np.isnan(orientations)) == True:
                             sys.stdout.write('Glasses: nan present')
                             self.store_pointer -= 1
                        elif np.any(np.isinf(positions)) == True:
                            sys.stdout.write('\n Glasses: inf present')
                            self.store_pointer -= 1
                        elif np.any(np.isinf(orientations)) == True:
                            sys.stdout.write('\n Glasses: inf present')
                            self.store_pointer -= 1
                        elif np.all(positions) == positions[0]:
                            sys.stdout.write('\n Glasses: Zeros present')
                            self.store_pointer -= 1
                        elif np.all(orientations) == orientations[0]:
                            sys.stdout.write('\n Glasses: Zeros present')
                            self.store_pointer -= 1 
                        else:
                            self.store_glasses += 1
                            
                    elif in_mdf.tool_id == (self.pointer + 1):
                        pointer_pos, Qr = self.pointer_tp.get_pos(orientations, positions)
                        self.pointer_pos[self.store_pointer, :] = pointer_pos                      
                        if np.any(np.isnan(positions)) == True:
                             sys.stdout.write('Pointer: nan present')
                             self.store_glasses -= 1
                        elif np.any(np.isnan(orientations)) == True:
                             sys.stdout.write('Pointer: nan present')
                             self.store_glasses -= 1
                        elif np.any(np.isinf(positions)) == True:
                            sys.stdout.write('\n Pointer: inf present')
                            self.store_glasses -= 1
                        elif np.any(np.isinf(orientations)) == True:
                            sys.stdout.write('\n Pointer: inf present')
                            self.store_glasses -= 1
                        elif np.all(positions) == positions[0]:
                            sys.stdout.write('\n Pointer: Zeros present')
                            self.store_glasses -= 1
                        elif np.all(orientations) == orientations[0]:
                            sys.stdout.write('\n Pointer: Zeros present')
                            self.store_glasses -= 1 
                        else:
                            self.store_pointer += 1
                            
                    if (self.store_pointer < 0) or (self.store_glasses < 0):
                        self.store_pointer = 0
                        self.store_glasses = 0
        else:
            pass
        
        
    def request(self, landmark):
        self.collect_positions = True
        self.landmark = landmark
        
    def shuffle_q(self, q):
        return np.roll(q, -1, axis=0)
 
    
 
if __name__=="__main__":
    parser = ArgumentParser(description = 'Real-time display of TMS Coil position')
    parser.add_argument(type=str, dest='config')
    parser.add_argument(type=str, dest='mm_ip', nargs='?', default='127.0.0.1:7111')
    args = parser.parse_args()
    app = wx.PySimpleApp()
    frame = CalibrateHead(args.config, args.mm_ip)
    app.MainLoop()