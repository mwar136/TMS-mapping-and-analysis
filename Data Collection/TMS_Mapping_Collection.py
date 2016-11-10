#Author: Matthew Ward mwar136@aucklanduni.ac.nz
'''
26/09/26 Changes:
Added dumps of all messages when start and stop buttons are pressed.
Message def auto copies to data directory
'''
import sys
import os
import logging
import time
import numpy as np
import wx
import shutil
import threading
from argparse import ArgumentParser
from IPython.config.loader import PyFileConfigLoader

import PyDragonfly
from PyDragonfly import CMessage, MT_EXIT, copy_to_msg, copy_from_msg
from dragonfly_utils import respond_to_ping

import Dragonfly_config as rc
import quaternionarray as qa
import amcmorl_py_tools.vecgeom as vg
from amcmorl_py_tools.vecgeom import transformations as tf

from TrackedPoint import TrackedPoint as TP

class TMS_Mapping_Collection(threading.Thread):
    def __init__(self, parent, config_file, server):
        threading.Thread.__init__(self)
        self.parent = parent
        self.serial_no = 2
        self.glasses_pos_buffer = np.zeros((1,3))
        self.glasses_orientation_buffer = np.zeros((1,4))
        self.coil_pos_buffer = np.zeros((1,3))
        self.coil_vector_buffer = np.zeros((1,3))
        self.trial_no = 0
        self.stimulus_no = 0
        self.status = False
        self.daemon = True
        self.testing = False
        self.calibration_loaded = False
        self.server = server
        self.load_config(config_file)
        self.load_logging()
        self.setup_dragonfly(server)
        self.initial_ping()
        self.find_data()
        '''
        self.calibration_matrix = self.make_head_axis(
                                                      self.Left_Tragus_pos, 
                                                      self.Right_Tragus_pos, 
                                                      self.Nasion_pos
        )
        self.glasses_to_cz = TP(
                                self.glasses_calib_ori,
                                self.cz_calib_pos,
                                self.glasses_calib_pos
        )
        '''
        self.calibration_loaded = True
        self.start()
            
    def find_data(self):
        app = wx.PySimpleApp()
        dialog = wx.FileDialog(None,  style=wx.FD_OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            self.file_path = dialog.GetPath() 
        dialog.Destroy()
        app = 0
        #self.file_path = "C:\\Users\\amcmorl\\Dave25-07-16.txt"
        self.open_file(self.file_path)
        self._copy_Mdef()
        
    def open_file(self, file_name):
        with open(file_name, 'r') as calib:
            Right_Tragus = []
            Left_Tragus = []
            Nasion = []
            Cz = []
            lines = calib.readlines()
            for i, line in enumerate(lines):
                if 'Right Tragus' in line:
                    Right_Tragus.extend((lines[i+1], lines[i+2], lines[i+3]))
                if 'Left Tragus' in line:
                    Left_Tragus.extend((lines[i+1], lines[i+2], lines[i+3]))
                if 'Nasion' in line:
                    Nasion.extend((lines[i+1], lines[i+2], lines[i+3]))
                if 'Cz' in line:
                    Cz.extend((lines[i+1], lines[i+2], lines[i+3]))
            '''
            self.Right_Tragus_pos = self._make_array(Right_Tragus, 2)
            self.Left_Tragus_pos = self._make_array(Left_Tragus, 2)
            self.Nasion_pos = self._make_array(Nasion, 2)
            self.cz_calib_pos = self._make_array(Cz, 2)
            self.glasses_calib_pos = self._make_array(Cz, 0)
            self.glasses_calib_ori =  self._make_array(Cz, 1)
            '''
            self.glasses_to_RT = TP(
                                    self._make_array(Right_Tragus, 1),
                                    self._make_array(Right_Tragus, 2),
                                    self._make_array(Right_Tragus, 0)
                                   )
            self.glasses_to_LT = TP(
                                    self._make_array(Left_Tragus, 1),
                                    self._make_array(Left_Tragus, 2),
                                    self._make_array(Left_Tragus, 0)
                                   )                       
            self.glasses_to_Nasion = TP(
                                    self._make_array(Nasion, 1),
                                    self._make_array(Nasion, 2),
                                    self._make_array(Nasion, 0)
                                   )
            self.glasses_to_cz = TP(
                                    self._make_array(Cz, 1),
                                    self._make_array(Cz, 2),
                                    self._make_array(Cz, 0)
                                   )
            
            
    def _make_array(self, arr, x): #x indexes array
        pos = arr[x].split()[-1]
        return np.array(pos.split(',')).astype(np.float)
         
    def make_head_axis(self, left_tragus, right_tragus, nasion):
        x = vg.unitvec(right_tragus-((left_tragus + right_tragus)/2)) #make right +x
        yk = vg.unitvec(nasion-((left_tragus + right_tragus)/2)) #make anterior +y
        zk = vg.unitvec(np.cross(x, yk))
        y =  vg.unitvec(np.cross(zk, x))
        z = np.cross(x,y) # make superior +z
        return tf.DCM2quat(np.array([x,
                                     y,
                                     z]))

    def load_config(self, config_file):
        cfg = PyFileConfigLoader(config_file)
        cfg.load_config()
        self.config = cfg.config
        
        # special casing for SAMPLE_GENERATED
        if (self.config.trigger == 'SAMPLE_GENERATED'):
            self.config.trigger_msg = rc.MT_SAMPLE_GENERATED
            self.config.trigger_mdf = rc.MDF_SAMPLE_GENERATED
        else:
            self.config.trigger_msg = \
                eval('rc.MT_' + self.config.trigger)
            self.config.trigger_mdf = \
                eval('rc.MDF_' + self.config.trigger)
        print "Triggering with", self.config.trigger
        print "TMS Mapping Collection: loading config"
        
        #self.ntools = len(self.config.tool_list)
        self.plate = self.config.tools.index('CB609')
        self.marker = self.config.tools.index('CT315')
        self.glasses = self.config.tools.index('ST568')
        self.pointer = self.config.tools.index('P717')
        self.metronome_config = self.config.metronome_config
    
    def setup_dragonfly(self, server):
        self.mod = PyDragonfly.Dragonfly_Module(0, 0)
        self.mod.ConnectToMMM(server)
        self.mod.Subscribe(MT_EXIT)
        self.mod.Subscribe(rc.MT_PING)
        self.mod.Subscribe(rc.MT_POLARIS_POSITION)
        self.mod.Subscribe(rc.MT_HOTSPOT_POSITION)
        self.mod.Subscribe(rc.MT_TMS_TRIGGER)
        
        self.mod.SendModuleReady()
        print "TMS_Mapping_Collection: connected to dragonfly"
        
    def load_logging(self):
        log_file = os.path.normpath(os.path.join(self.config.config_dir, 'TMS_mapping_collection.log'))
        print "log file: " + log_file
        logging.basicConfig(filename=log_file, level=logging.DEBUG)
        logging.info(' ')
        logging.info(' ')
        logging.debug("**** STARTING UP ****")
        logging.info("  %s  " % time.asctime())
        logging.info("*********************")
    
    def initial_ping(self):
        t = time.time()
        while time.time() < t + 10:
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
                    respond_to_ping(self.mod, msg, 'TMS_Mapping_Collection')
                    break
                else:
                    pass
    
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
                    respond_to_ping(self.mod, msg, 'TMS_Mapping_Collection')
                else:
                    self.process_message(msg)
    
    def process_message(self, in_msg):
        msg_type = in_msg.GetHeader().msg_type
        #print('? %d STATUS=%s TESTING=%s' % (msg_type, str(self.status), str(self.testing)))
        if self.status == False:
                pass          
        elif self.testing:
            if msg_type == rc.MT_POLARIS_POSITION:
                # handling input message
                in_mdf = rc.MDF_POLARIS_POSITION()
                copy_from_msg(in_mdf, in_msg)
                positions = np.array(in_mdf.xyz[:])
                orientations = qa.norm(self.shuffle_q(np.array(in_mdf.ori[:])))
                if self.Ptest_count == 250 & self.Gtest_count == 250:
                    self.testing = False
                    self.status = False
                    print 'done'
                else:
                    if in_mdf.tool_id == (self.pointer + 1): 
                        self.pointer_test[self.Ptest_count, :] = positions
                        if np.any(np.isnan(positions)) == True:
                            sys.stdout.write('\n Pointer: NaN present')
                            self.Gtest_count -=1
                        elif np.any(np.isnan(orientations)) == True:
                            sys.stdout.write('\n Pointer: NaN present')
                            self.Gtest_count -=1
                        elif np.any(np.isinf(positions)) == True:
                            sys.stdout.write('\n Pointer: inf present')
                            self.Gtest_count -=1
                        elif np.any(np.isinf(orientations)) == True:
                            sys.stdout.write('\n Pointer: inf present')
                            self.Gtest_count -=1
                        elif np.all(positions) == positions[0]:
                            sys.stdout.write('\n Pointer: Zeros present')
                            self.Gtest_count -=1
                        elif np.all(orientations) == orientations[0]:
                            sys.stdout.write('\n Pointer: Zeros present')
                            self.Gtest_count -=1
                        else:
                            self.Ptest_count += 1
                            
                    if in_mdf.tool_id == (self.glasses + 1):
                        self.headP_test[self.Gtest_count, :] = positions
                        self.headQ_test[self.Gtest_count, :] = orientations
                        if np.any(np.isnan(positions)) == True:
                            sys.stdout.write('\n Glasses: NaN present')
                            self.Ptest_count -=1
                        elif np.any(np.isnan(orientations)) == True:
                            sys.stdout.write('\n Glasses: NaN present')
                            self.Ptest_count -=1
                        elif np.any(np.isinf(positions)) == True:
                            sys.stdout.write('\n Glasses: inf present')
                            self.Ptest_count -=1
                        elif np.any(np.isinf(orientations)) == True:
                            sys.stdout.write('\n Glasses: inf present')
                            self.Ptest_count -=1
                        elif np.all(positions) == positions[0]:
                            sys.stdout.write('\n Glasses: Zeros present')
                            self.Ptest_count -=1
                        elif np.all(orientations) == orientations[0]:
                            sys.stdout.write('\n Glasses: Zeros present')
                            self.Ptest_count -=1
                        else:
                            self.Gtest_count += 1
                            print orientations
                
                    if (self.Ptest_count < 0) or (self.Gtest_count < 0):
                        self.Gtest_count = 0
 
        else:
            if msg_type == rc.MT_POLARIS_POSITION:
                # handling input message
                in_mdf = rc.MDF_POLARIS_POSITION()
                copy_from_msg(in_mdf, in_msg)
                if in_mdf.tool_id == (self.glasses + 1):
                    positions = np.array(in_mdf.xyz[:])
                    orientations = qa.norm(self.shuffle_q(np.array(in_mdf.ori[:])))
                    self.glasses_pos_buffer = positions
                    self.glasses_orientation_buffer = orientations
                else:
                    pass
            elif msg_type == rc.MT_HOTSPOT_POSITION:
                in_mdf = rc.MDF_HOTSPOT_POSITION()
                copy_from_msg(in_mdf, in_msg)
                positions = np.array(in_mdf.xyz[:])
                vector = np.array(in_mdf.ori[:3])
                self.coil_pos_buffer = positions
                self.coil_vector_buffer = vector
                
                    
            elif msg_type == rc.MT_TMS_TRIGGER:
                self.find_hotspot_to_cz(
                                        self.coil_pos_buffer, 
                                        self.glasses_orientation_buffer,
                                        self.glasses_pos_buffer,
                                        self.coil_vector_buffer
                )
                
    def shuffle_q(self, q):
        return np.roll(q, -1, axis=0)
        
    def find_hotspot_to_cz(self, HS_GCS, glasses_orientation, glasses_position, coil_vector): #Saves 0, 0, 0 array to file when object called
        check_HS_GCS = self.check(HS_GCS)
        check_glasses_orientation = self.check(glasses_orientation)
        check_glasses_position = self.check(glasses_position)
        check_coil_vector = self.check(coil_vector)
        if np.any(np.array([check_HS_GCS, check_glasses_orientation, check_glasses_position, check_coil_vector])):
            print np.array([check_HS_GCS, check_glasses_orientation, check_glasses_position, check_coil_vector])
            HS_LCS = (np.zeros(3)* np.NaN)[None]
            print HS_LCS
            vector_LCS = (np.zeros(3)* np.NaN)[None]
            with open(self.location, 'a') as location:
                np.savetxt(location, np.insert(np.array(HS_LCS), 0, self.trial_count)[None], fmt='%f', delimiter=',', newline='\r\n')
                location.close()
            with open(self.vector, 'a') as data:
                np.savetxt(data, np.insert(np.array(vector_LCS), 0, self.trial_count)[None], fmt='%f', delimiter=',', newline='\r\n')
                data.close()
            sys.stdout.write("No position\n")
        else:
            cz_pos, cz_rot = self.glasses_to_cz.get_pos(glasses_orientation, glasses_position)
            LT_pos, LT_rot = self.glasses_to_LT.get_pos(glasses_orientation, glasses_position)
            RT_pos, RT_rot = self.glasses_to_RT.get_pos(glasses_orientation, glasses_position)
            Nasion_pos, Nasion_rot = self.glasses_to_Nasion.get_pos(glasses_orientation, glasses_position)
            calibration_matrix = self.make_head_axis(
                                                          LT_pos, 
                                                          RT_pos, 
                                                          Nasion_pos)
            HS_LCS =  qa.rotate(calibration_matrix, (HS_GCS - cz_pos))
            vector_LCS =  qa.rotate(calibration_matrix, coil_vector)
            
            #HS_LCS =  qa.rotate(rot, qa.rotate(self.calibration_matrix, (HS_GCS - cz_pos)))
            #vector_LCS =  qa.rotate(rot, qa.rotate(self.calibration_matrix, (coil_vector)))
            with open(self.location, 'a') as location:
                np.savetxt(location, np.insert(np.array(HS_LCS), 0, self.trial_count)[np.newaxis], fmt='%f', delimiter=',', newline='\r\n')
                location.close()
            with open(self.vector, 'a') as data:
                np.savetxt(data, np.insert(np.array(vector_LCS), 0, self.trial_count)[None], fmt='%f', delimiter=',', newline='\r\n')
                data.close()
            self.stimulus_no += 1
            print self.stimulus_no
              
        in_mdf = rc.MDF_POLARIS_POSITION()
        out_mdf = rc.MDF_PLOT_POSITION()
        out_mdf.xyz[:] = HS_LCS[0]
        out_mdf.ori[:] = np.append(vector_LCS[0], 0)# Qk - coil active orientation
        out_mdf.sample_header = in_mdf.sample_header
        msg = CMessage(rc.MT_PLOT_POSITION)
        copy_to_msg(out_mdf, msg)
        self.mod.SendMessage(msg)
        self.trial_count += 1  
        sys.stdout.write("C")

    def collect(self):
        file_path = self.file_path.split('\\')[:-1]
        self.location = '\\'.join(file_path) + '\\coil_location' + str(self.trial_no ) + '.txt'
        self.vector = '\\'.join(file_path) + '\\coil_vector'+ str(self.trial_no ) + '.txt'
        self.status = True
        self.trial_count = 0
        self._flush_quicklogger(self.trial_no)
        self.send_metronome(1)
        
    def pause(self):
        self.status = False
        self.send_metronome(2)
        
    def restart(self):
        if self.status == False:
            self.status = True
            self.send_metronome(1)
        else:
            print 'Already running'
        
    def stop(self):
        self.status = False
        self._save_file()
        self.coil_pos_buffer[:] = np.NaN
        self.glasses_pos_buffer[:] = np.NaN
        self.glasses_orientation_buffer[:] = np.NaN
        self.coil_vector_buffer[:] = np.NaN
        self.trial_count = 0
        self.trial_no += 1
        self.send_metronome(0)
        self.stimulus_no = 0

        
    def send_metronome (self, state):
        mdf = rc.MDF_MNOME_STATE()
        self.serial_no += 1
        mdf.sample_header.SerialNo  = self.serial_no
        mdf.sample_header.Flags     = 0
        mdf.sample_header.DeltaTime = (1. / 5)
        mdf.state = state
        msg = CMessage(rc.MT_MNOME_STATE)
        copy_to_msg(mdf, msg)
        self.mod.SendMessage(msg)
        print "Sent message %d" %state
    
    def check(self, data):
        if np.any(np.isnan(data)) == True:
            return True
        elif np.any(np.isnan(data)) == True:
            return True
        elif np.any(np.isinf(data)) == True:
            return True
        elif np.any(np.isinf(data)) == True:
            return True
        elif np.all(data == data[0]):
            return True
        elif np.all(data == data[0]):
            return True
        else:
            return False
    
    def test_calibration(self, axis):
        if self.calibration_loaded:
            self.pointer_test = np.zeros((250, 3))
            self.headP_test = np.zeros((250, 3))
            self.headQ_test = np.zeros((250, 4))
            self.Ptest_count = 0
            self.Gtest_count = 0
            self.status = True
            self.testing = True
            test_result = np.zeros((250,3))
            print 'doing calculations'
            for i in np.arange(250):
                cz_pos, rot = self.glasses_to_cz.get_pos(self.headQ_test[i],self.headP_test[i])
                test_result[i] =  qa.rotate(rot, qa.rotate(self.calibration_matrix, (self.pointer_test[i] - cz_pos)))
            print 'complete'
            if axis == 'x':  
                diff = np.diff(test_result[~np.isnan(test_result).any(axis=1), 0],axis=0)
                if np.any(diff < 0):
                    return False
                else:
                    return True
            
            if axis == 'y':  
                diff = np.diff(test_result[~np.isnan(test_result).any(axis=1), 1],axis=0)
                if np.any(diff < 0):
                    return False
                else:
                    return True
                    
            if axis == 'z':  
                diff = np.diff(test_result[~np.isnan(test_result).any(axis=1), 2],axis=0)
                if np.any(diff < 0):
                    return False
                else:
                    return True
            
            self.status = False
            self.testing = False

    def _save_file(self):
        file_path = ('\\').join(self.file_path.split('\\')[:-1])
        subject = self.file_path.split('\\')[-1].split('.')[0]
        msg = CMessage(PyDragonfly.MT_SAVE_MESSAGE_LOG)
        mdf = PyDragonfly.MDF_SAVE_MESSAGE_LOG()
        ql_pathname  = os.path.join(file_path,
            '%s.%d.bin' % (subject, self.trial_no))
        mdf.pathname = str(ql_pathname)
        mdf.pathname_length = len(ql_pathname)
        msg = CMessage(PyDragonfly.MT_SAVE_MESSAGE_LOG)
        msg.SetData(mdf.this, 260)
            # NOT SURE IF THIS NUMBER IS RIGHT - SEEMS TO WORK
            # NEED TO WRAP DRAGONFLY FOR PYTHON BETTER!
        self.mod.SendMessage(msg)
        
    def _flush_quicklogger(self, id):
        file_path = ('\\').join(self.file_path.split('\\')[:-1])
        ql_pathname = os.path.join(file_path, 'ql_flush.%d.bin' % (id))
        print ql_pathname
        mdf = PyDragonfly.MDF_SAVE_MESSAGE_LOG()
        mdf.pathname = str(ql_pathname)
        mdf.pathname_length = len(ql_pathname)
        msg = CMessage(PyDragonfly.MT_SAVE_MESSAGE_LOG)
        msg.SetData(mdf.this, 260)
        self.mod.SendMessage(msg)
    
    def _copy_Mdef(self):
        copy_dir = ('\\').join(self.file_path.split('\\')[:-1])
        shutil.copy('C:\\Users\\amcmorl\\Documents\\BCI\\message_definitions\\Dragonfly_config.mat', copy_dir)
    
class MainWindow(wx.Frame):

    def __init__(self, id, config_file, server):
        wx.Frame.__init__(self, None, id, 'TMS mapping collection', size=(420,200))

        test_calibration_button = wx.Button(self, 1, 'Test Calibration', size=(100,20), pos=(10,120))
        test_calibration_button.Bind(wx.EVT_BUTTON, self.test_calibration, id=1)
        '''
        calibrate_coil_button = wx.Button(self, 2, 'Calibrate Coil', size=(80,20), pos=(110,80))
        calibrate_coil_button.Bind(wx.EVT_BUTTON, self.create_file, id=2)
        '''
        self.start_trial_button = wx.Button(self, 3, 'Start Trial', size=(80,60), pos=(10,40))
        self.pause_trial_button = wx.Button(self, 4, 'Pause Trial', size=(80,60), pos=(110,40))
        self.stop_trial_button = wx.Button(self, 5, 'Stop Trial', size=(80,60), pos=(210,40))
        self.restart_trial_button = wx.Button(self, 6, 'Restart Trial', size=(80,60), pos=(310,40))
                
        self.trial_no_text = wx.StaticText(self, 8, str(0), pos=(150, 120), size=(90,90))

        self.mapping_collection = TMS_Mapping_Collection(self, config_file, server)        
        self.start_trial_button.Bind(wx.EVT_BUTTON, self.start_trial, id=3)
        self.pause_trial_button.Bind(wx.EVT_BUTTON, self.pause, id=4)
        self.stop_trial_button.Bind(wx.EVT_BUTTON, self.stop, id=5)
        self.restart_trial_button.Bind(wx.EVT_BUTTON, self.restart , id=6)

 
        
    def start_trial(self, event):
        self.mapping_collection.collect()
        self.start_trial_button.Disable()
        self.restart_trial_button.Disable()
        self.stop_trial_button.Enable()
        self.pause_trial_button.Enable()
       
    
    def stop(self, event):
        self.mapping_collection.stop()
        self.start_trial_button.Enable()
    
         
    def pause(self, event):
        self.mapping_collection.pause()
        self.restart_trial_button.Enable()
        
    
    def restart(self, event):
        self.mapping_collection.restart()
        self.restart_trial_button.Disable()
       
    
    def test_calibration (self,event):

        if self.mapping_collection.calibration_loaded:
            app = wx.App()
            dlg = wx.MessageDialog(None, 'Move pointer from the left ear to the right ear', 'Check Calibration', 
            wx.OK|wx.ICON_INFORMATION)
            dlg.ShowModal()
            dlg.Destroy()
            result = self.mapping_collection.test_calibration('x')
            if result:
                sys.stdout.write('X axis is good')
            else:
                sys.stdout.write('X axis is bad')
            app = wx.PySimpleApp()        
            dlg = wx.MessageDialog(None, 'Move pointer from the vertex towards the nasion', 'Check Calibration', 
            wx.OK|wx.ICON_INFORMATION)
            dlg.ShowModal()
            if out == wx.OK:
                dlg.Destroy()
                result = self.mapping_collection.test_calibration('y')
                if result:
                    print 'Y axis is good'
                else:
                    print 'Y axis is bad'
            app = wx.PySimpleApp()
            dlg = wx.MessageDialog(None, 'Move pointer from the vertex to the right ear', 'Check Calibration', 
            wx.OK | wx.ICON_INFORMATION)
            out = dlg.ShowModal()
           
            if out == wx.OK:
                dlg.Destroy()
                result = self.mapping_collection.test_calibration('z')
                if result:
                    print 'Z axis is good'
                else:
                    print 'Z axis is bad'
        else:
            app = wx.PySimpleApp()
            dlg = wx.MessageDialog(None, 'Load calibration file first!', 'Error', 
            wx.OK | wx.ICON_ERROR)
            out = dlg.ShowModal()
            dlg.Destroy()

    
if __name__=="__main__":
    parser = ArgumentParser(description = 'Trial controller for TMS mapping')
    parser.add_argument(type=str, dest='config')
    parser.add_argument(type=str, dest='mm_ip', nargs='?', default='127.0.0.1:7111')
    args = parser.parse_args()
    app = wx.App()
    frame = MainWindow(wx.ID_ANY, args.config, args.mm_ip)
    frame.Show()
    app.MainLoop()