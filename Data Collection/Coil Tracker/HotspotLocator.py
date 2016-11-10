#Author: Matthew Ward mwar136@aucklanduni.ac.nz

import sys
import os
import logging
import time
import numpy as np
from argparse import ArgumentParser
from IPython.config.loader import PyFileConfigLoader


import PyDragonfly
from PyDragonfly import CMessage, MT_EXIT, copy_to_msg, copy_from_msg
from dragonfly_utils import respond_to_ping

import Dragonfly_config as rc
import quaternionarray as qa
import amcmorl_py_tools.vecgeom as vg

from TrackedPoint import TrackedPoint as TP
import quaternion_average as mean

import threading
import wx
from PySide import QtGui, QtCore

'''
Calibrates the TMS coil and calculates the location of the hotspot in space
Requires 2 seconds to calibrate coil

Order of operations:
1. Loads config file
2. Begins logging
3. Setup dragonfly server
4. User begins calibration with enter- any other response is not accepted
5. Data storage array is created- needs information to be added
6. Run begins and starts message processing
7. Sample alignment is checked- If okay vector is calculated
8. Calibration becomes true so hotspot is now calculated in space using:

Ni = Calibration plate position at calibration time
Ti = Coil marker position at calibration time

Xi = Ni - Ti

Qi = Coil marker orientation at calibration time
Qk = Coil marker orientation at Tc+N
Qr = Rotation from calibration to Tc+N

Qr = Qk * Qi'

Tk = Position of coil at marker
Xk = Position of hotspot at Tc+N

Xk = (Qr*Xi)*Qr' + Tk
'''

class ChoiceDialog(wx.Frame):
    
    def __init__(self, parent, variable, id):
        wx.Frame.__init__(self, parent, id, 'Hotspot Locator', size=(320,100))
        new_calibration = wx.Button(self, 1, 'New Calibration', size=(100,20), pos=(30,20))
        load_calibration = wx.Button(self, 2, 'Load Calibration', size=(100,20), pos=(150,20))
        self.variable = variable
        new_calibration.Bind(wx.EVT_BUTTON, self.newCalibration)
        load_calibration.Bind(wx.EVT_BUTTON, self.loadCalibration)
        

    def newCalibration(self, event):
        self.variable.choice = 'New'
        self.Close()
        
    def loadCalibration(self, event):
        self.variable.choice = 'Load'
        self.Close()
        
class HotspotLocator (object):
    
    def __init__(self, config_file, mm_ip):
        self.plate_vector = np.array([1.6, 85.4, 2.4])
        self.plate_init_ori = np.array([0,0,0,-1])
        self.calibrating = False
        self.calibrated = False
        self.load_config(config_file)
        self.load_logging()
        self.setup_dragonfly(mm_ip)
        self.get_frequency()
        t = threading.Thread(target=self.run)
        app = wx.PySimpleApp()
        choice = ChoiceDialog(None, self, wx.ID_ANY)
        choice.Show(True)
        app.MainLoop()
        if self.choice == 'Load':
            self.load_old_calibration()
            self.calibrated = True
            t.start()
        else:    
            self.create_storage()
            t.start()
            self.user_start_calibrate()
        
    def load_old_calibration(self):
        app = wx.PySimpleApp()
        dialog = wx.FileDialog(None,  style=wx.FD_OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            file_path = dialog.GetPath()
        dialog.Destroy()
        with open(file_path, 'r') as calib:
            lines = calib.readlines()
            self.Qi = np.fromstring(lines[0], sep=',')
            Ni = np.fromstring(lines[1], sep=',')
            Ti = np.fromstring(lines[2], sep=',')
            self.cal_plate_vec = np.fromstring(lines[3], sep=',')
        self.tp = TP(self.Qi, Ni, Ti)
        
    def load_config(self, config_file):
        cfg = PyFileConfigLoader(config_file)
        cfg.load_config()
        self.config = cfg.config
        print "HotspotLocator: loading config"      
        #self.ntools = len(self.config.tool_list)
        self.plate = self.config.tools.index('CB609')
        self.marker = self.config.tools.index('CT315')
        self.glasses = self.config.tools.index('ST568')
        self.pointer = self.config.tools.index('P717')
         
    def load_logging(self):
        log_file = os.path.normpath(os.path.join(self.config.config_dir, 'coil_calibration.log'))
        print "log file: " + log_file
        logging.basicConfig(filename=log_file, level=logging.DEBUG)
        logging.info(' ')
        logging.info(' ')
        logging.debug("**** STARTING UP ****")
        logging.info("  %s  " % time.asctime())
        logging.info("*********************")
            
    def setup_dragonfly(self, mm_ip):
        self.mod = PyDragonfly.Dragonfly_Module(0, 0)
        self.mod.ConnectToMMM(mm_ip)
        self.mod.Subscribe(MT_EXIT)
        self.mod.Subscribe(rc.MT_PING)
        self.mod.Subscribe(rc.MT_POLARIS_POSITION)
        
        self.mod.SendModuleReady()
        print "HotspotLocator: connected to dragonfly"
    
    def get_frequency(self):
        # loop over receiving messages until we get a POLARIS_POSITION message
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
                    respond_to_ping(self.mod, msg, 'HotspotLocator')
                else:
                    msg_type = msg.GetHeader().msg_type
                    if msg_type == rc.MT_POLARIS_POSITION:
                        # handling input message
                        mdf = rc.MDF_POLARIS_POSITION()
                        copy_from_msg(mdf, msg)
                        self.fsamp = 1/mdf.sample_header.DeltaTime
                        if self.fsamp != 0:
                            print 'Got Frequency!'
                            break
        # (handle EXITS and PINGS appropriately)
        # return 1 / DeltaTime from first POLARIS_POSITION msg
        
        
    
    def user_start_calibrate(self):
        # get a POLARIS_POSITION message, read sample_header.DeltaTime to get
        # message frequency
        while True:
            x = raw_input("Press enter to calibrate...\n\r")
            if not x:
                break
            print '.......'
        sys.stdout.write('starting in:')
        sys.stdout.write('5\n')
        sys.stdout.flush()
        time.sleep(1)
        sys.stdout.write('4\n')
        sys.stdout.flush()
        time.sleep(1)
        sys.stdout.write('3\n')
        sys.stdout.flush()
        time.sleep(1)
        sys.stdout.write('2\n')
        sys.stdout.flush()
        time.sleep(1)
        sys.stdout.write('1\n')
        sys.stdout.flush()
        time.sleep(1)
        sys.stdout.write('Calibrating...')
        self.calibrating = True
    
    def create_storage (self):
        self.store_plate_pos = np.empty([5 * self.fsamp, 3])
        self.store_plate_ori = np.empty([5 * self.fsamp, 4])
        self.store_coil_pos = np.empty([5 * self.fsamp, 3])
        self.store_coil_ori = np.empty([5 * self.fsamp, 4])
        self.store_plate = 0
        self.store_coil = 0
        
        
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
                    respond_to_ping(self.mod, msg, 'HotspotLocator')
                else:
                    self.process_message(msg)
    
    def process_message(self, in_msg):
        msg_type = in_msg.GetHeader().msg_type
        if msg_type == rc.MT_POLARIS_POSITION:
            # handling input message
            in_mdf = rc.MDF_POLARIS_POSITION()
            copy_from_msg(in_mdf, in_msg)
            positions = np.array(in_mdf.xyz[:])
            orientations = qa.norm(self.shuffle_q(np.array(in_mdf.ori[:])))
                           
            #np.testing.assert_array_equal(positions[:,0], orientations[:,0], err_msg='Samples are not aligned')
            
            if self.calibrated:
                if in_mdf.tool_id == (self.marker + 1): 
                    # calculating output
                    hotspot_position, Qr = self.tp.get_pos(orientations, positions)
                    hotspot_vector_head = qa.rotate(Qr, self.cal_plate_vec)
                    if np.any(np.isnan(hotspot_position)) == True:
                         sys.stdout.write('x')
                         hotspot_vector_head[:] = np.NaN
                         hotspot_position[:] = np.NaN
                    elif np.any(np.isnan(orientations)) == True:
                         sys.stdout.write('x')
                         hotspot_vector_head[:] = np.NaN
                         hotspot_position[:] = np.NaN
                    elif np.all(orientations) == orientations[0]:
                        sys.stdout.write('x')
                        hotspot_vector_head[:] = np.NaN
                        hotspot_position[:] = np.NaN
                    elif np.any(np.isinf(orientations)) == True:
                        sys.stdout.write('x')
                        hotspot_vector_head[:] = np.NaN
                        hotspot_position[:] = np.NaN
                          #print '         *****nan present, check coil is within frame!*****'
                    
                        #creating output message
                    out_mdf = rc.MDF_HOTSPOT_POSITION()
                    out_mdf.xyz[:] = hotspot_position
                    out_mdf.ori[:] = np.append(hotspot_vector_head, 0)# Qk - coil active orientation
                    out_mdf.sample_header = in_mdf.sample_header
                    msg = CMessage(rc.MT_HOTSPOT_POSITION)
                    copy_to_msg(out_mdf, msg)
                    self.mod.SendMessage(msg)
                    sys.stdout.write("o")
                                    
            else:
                if self.calibrating:
                    if(
                       (self.store_plate >= self.store_plate_pos.shape[0]) & 
                       (self.store_plate >= self.store_plate_ori.shape[0]) & 
                       (self.store_coil >= self.store_coil_pos.shape[0]) & 
                       (self.store_coil >= self.store_coil_ori.shape[0])
                    ):
                        #If all the arrays are full calibration is complete
                        self.calibrating = False
                        self.make_calibration_vector()
                    else:
                        if in_mdf.tool_id == (self.marker + 1): 
                            self.store_coil_pos[self.store_coil, :] = positions
                            self.store_coil_ori[self.store_coil, :] = orientations
                            if np.any(np.isnan(positions)) == True:
                                sys.stdout.write('\n Coil: NaN present')
                                self.store_plate -=1
                            elif np.any(np.isnan(orientations)) == True:
                                sys.stdout.write('\n Coil: NaN present')
                                self.store_plate -=1
                            elif np.any(np.isinf(positions)) == True:
                                sys.stdout.write('\n Coil: inf present')
                                self.store_plate -=1
                            elif np.any(np.isinf(orientations)) == True:
                                sys.stdout.write('\n Coil: inf present')
                                self.store_plate -=1
                            elif np.all(positions) == positions[0]:
                                sys.stdout.write('\n Coil: Zeros present')
                                self.store_plate -=1
                            elif np.all(orientations) == orientations[0]:
                                sys.stdout.write('\n Coil: Zeros present')
                                self.store_plate -=1
                            else:
                                self.store_coil += 1
                                
                        if in_mdf.tool_id == (self.plate + 1):
                            self.store_plate_pos[self.store_plate, :] = positions
                            self.store_plate_ori[self.store_plate, :] = orientations
                            if np.any(np.isnan(positions)) == True:
                                sys.stdout.write('\n Plate: NaN present')
                                self.store_coil -=1
                            elif np.any(np.isnan(orientations)) == True:
                                sys.stdout.write('\n Plate: NaN present')
                                self.store_coil -=1
                            elif np.any(np.isinf(positions)) == True:
                                sys.stdout.write('\n Plate: inf present')
                                self.store_coil -=1
                            elif np.any(np.isinf(orientations)) == True:
                                sys.stdout.write('\n Plate: inf present')
                                self.store_coil -=1
                            elif np.all(positions) == positions[0]:
                                sys.stdout.write('\n Plate: Zeros present')
                                self.store_coil -=1
                            elif np.all(orientations) == orientations[0]:
                                sys.stdout.write('\n Plate: Zeros present')
                                self.store_coil -=1
                            else:
                                self.store_plate += 1
                                print orientations
                
                        if (self.store_plate < 0) or (self.store_coil < 0):
                            self.store_plate = 0
                            self.store_coil = 0
                else:
                    pass
        
        
        else:
            pass
    def make_calibration_vector(self):
        plate_ori = mean.quaternion_mean(self.store_plate_ori)
        Ni        = self.store_plate_pos.mean(axis=0)
        self.Qi   = mean.quaternion_mean(self.store_coil_ori)
        Ti        = self.store_coil_pos.mean(axis=0)
        self.tp = TP(self.Qi, Ni, Ti)
        plate_rot = qa.mult(plate_ori, qa.inv(self.plate_init_ori)).flatten()
        self.cal_plate_vec = qa.rotate(plate_rot, self.plate_vector).flatten()
        msg_str_pos = "%.5e, " * 3
        msg_str_ori = "%.5e, " * 4
        sys.stdout.write('\n Plate orientation:    ')
        sys.stdout.write(msg_str_ori % (plate_ori[0], plate_ori[1], plate_ori[2],\
                                    plate_ori[3]))
        sys.stdout.write('\n Plate position:       ')
        sys.stdout.write(msg_str_pos % (Ni[0], Ni[1], Ni[2]))
        sys.stdout.write('\n Coil orientation:     ')
        sys.stdout.write(msg_str_ori % (self.Qi[0], self.Qi[1], self.Qi[2], self.Qi[3]))
        sys.stdout.write('\n Coil position:        ')
        sys.stdout.write(msg_str_pos % (Ti[0], Ti[1], Ti[2]))
        sys.stdout.write("\n********** Calibration complete! ***********\n\r")
        sys.stdout.flush()
        app = wx.PySimpleApp()
        saveFileDialog = wx.FileDialog(None, 'Save File', os.getcwd(), "", "Text File (*.txt)|*.txt",
                                       wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        saveFileDialog.ShowModal()
        app.MainLoop()
        self.filename = saveFileDialog.GetPath()
        saveFileDialog.Destroy()
        print self.filename
        with open(self.filename, 'a') as f:
            np.savetxt(f, self.Qi[None], fmt='%f' ,delimiter=',')
            np.savetxt(f, Ni[None], fmt='%f' ,delimiter=',')
            np.savetxt(f, Ti[None], fmt='%f' ,delimiter=',')
            np.savetxt(f, self.cal_plate_vec[None], fmt='%f' ,delimiter=',')
            f.close()
        self.calibrated = True
    
    def shuffle_q(self, q):
        return np.roll(q, -1, axis=0)

if __name__ == "__main__":
    parser = ArgumentParser(description = 'Interface with Polaris hardware' \
        ' and emit HOTSPOT_POSITION messages')
    parser.add_argument(type=str, dest='config')
    parser.add_argument(type=str, dest='mm_ip', nargs='?', default='')
    args = parser.parse_args()
    print("Using config file=%s, MM IP=%s" % (args.config, args.mm_ip))
    pdf = threading.Thread(target=HotspotLocator, args=(args.config, args.mm_ip))
    pdf.start()
    print "Finishing up"