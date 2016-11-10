#Author: Matthew Ward mwar136@aucklanduni.ac.nz
'''
start from command line
Current bug:
Random dropping of single EMG recording for one trial- current fix has an output of trial count from when 
metronome signal is received and a second count when the data is save to a file. 
These should only be different by 1.

Currently the bug does not appear to be within the BackgroundProcess object

'''
import sys

import numpy as np

from ConfigParser import SafeConfigParser
from argparse import ArgumentParser

from PyDragonfly import Dragonfly_Module, MT_EXIT, CMessage, copy_to_msg, copy_from_msg
import Dragonfly_config as rc
from dragonfly_utils import respond_to_ping

from PySide import QtGui, QtCore
import pyqtgraph as pg

import pdb
import threading
import DAQ_IO_Triggers as dit
import time

class Config(object):
    pass

class FastDisplay(QtGui.QMainWindow):
    def __init__(self, config_file, server):
        #self.parent = parent
        #self.counter = 0
        QtGui.QMainWindow.__init__(self)
        self.load_config(config_file)
        self.current_tab = 'RMS'
        #self.setup_dragonfly(server)
        #self.respond_to_ping()
        self.filename = QtGui.QFileDialog.getSaveFileName(self, "Save file", "", ".csv")[0] #name initial file in participant folder, will iterate for each trial
        self.bgprocess = BackgroundProcess(self, config_file,server)
        #self.TMS_trigger = False
        self.init_gui()
        #self.count = 0
        self.file_number = 0
        self.trial_count = 0
        

    def load_config(self, config_file):
        cfg = SafeConfigParser()
        cfg.read(config_file)
        self.config = Config()
        #daq_config.minV  = cfg.getfloat('main', 'minV')
        #daq_config.maxV  = cfg.getfloat('main', 'maxV')
        self.config.nsamp = cfg.getint('main', 'nsamp_per_chan_per_second')
        self.config.nchan = cfg.getint('main', 'nchan')
        self.config.nemg = cfg.getint('main', 'nemg')
        self.config.nirq  = self.freq = cfg.getint('main', 'nirq_per_second')
        self.config.pre_trig = cfg.getfloat('main', 'pre_trigger')
        self.config.perchan = self.config.nsamp / self.config.nirq
        self.config.npt   = self.config.nsamp * self.config.nchan / self.config.nirq
        self.config.pre_trig_samp = self.config.pre_trig * self.config.nsamp
        
        sys.stdout.write("nsamp: " + str(self.config.nsamp))
        sys.stdout.write("nchan: " + str(self.config.nchan))
        sys.stdout.write("nirq: " + str(self.config.nirq))
        sys.stdout.write("perchan: " + str(self.config.perchan))
        sys.stdout.write("npt: " + str(self.config.npt))
        
        assert((self.config.nsamp * self.config.nchan) % self.config.nirq == 0)
        assert(self.config.nsamp % self.config.nirq == 0)
           
    def init_gui(self):
        
        tabs = QtGui.QTabWidget()
        
        rms_win = pg.GraphicsLayoutWidget()
        self.hotspot_win = pg.GraphicsLayoutWidget()
        collect_win = pg.GraphicsLayoutWidget()
        #win.resize(1000,600)
        
        self.create_menu()
        self.setMenuBar(self.menu)
        
        cols = 2
        rows = self.config.nemg / cols
        self.npt = 2500 #Should find a way to dynamically set this, impacts viewed time
        self.amplification = np.ones((self.config.nemg))*1000
        self.time_array = np.arange(self.config.pre_trig*(-1), (self.npt/2000.0)-self.config.pre_trig, 1.0/self.config.nsamp)
        
        save_box = QtGui.QCheckBox('Save?')
        self.saving_data = False
        save_box.stateChanged.connect(self.save_data)
        
        pick_chan = QtGui.QComboBox()
        pick_chan.currentIndexChanged.connect(lambda i: self.change_chan(i))
        #################################################################################
        ############ Creating RMS and Amplitude window boxes for Hotspot################
        #################################################################################
        
        self.hotspot_amp_start_win = QtGui.QDoubleSpinBox()
        self.hotspot_amp_start_win.setRange(self.time_array.min(), self.time_array.max())
        self.hotspot_amp_start_win.setValue(0.01)
        self.hotspot_amp_start_win.setDecimals(4)
        self.hotspot_amp_start_win.setSingleStep(0.001)
        hotspot_amp_start_win_label = QtGui.QLabel('P-P start')
        
        self.hotspot_amp_stop_win = QtGui.QDoubleSpinBox()
        self.hotspot_amp_stop_win.setRange(self.time_array.min(), self.time_array.max())
        self.hotspot_amp_stop_win.setValue(0.045)
        self.hotspot_amp_stop_win.setDecimals(4)
        self.hotspot_amp_stop_win.setSingleStep(0.001)
        hotspot_amp_stop_win_label = QtGui.QLabel('P-P stop')
        
        self.hotspot_rms_start_win = QtGui.QDoubleSpinBox()
        self.hotspot_rms_start_win.setRange(self.time_array.min(), self.time_array.max())
        self.hotspot_rms_start_win.setValue(-0.05)
        self.hotspot_rms_start_win.setDecimals(4)
        self.hotspot_rms_start_win.setSingleStep(0.001)
        hotspot_rms_start_win_label = QtGui.QLabel('RMS window start')
        
        self.hotspot_rms_stop_win = QtGui.QDoubleSpinBox()
        self.hotspot_rms_stop_win.setRange(self.time_array.min(), self.time_array.max())
        self.hotspot_rms_stop_win.setValue(0)
        self.hotspot_rms_stop_win.setDecimals(4)
        self.hotspot_rms_stop_win.setSingleStep(0.001)
        hotspot_rms_stop_win_label = QtGui.QLabel('RMS window stop')

        #################################################################################
        ############ Creating RMS and Amplitude window boxes for Collect################
        #################################################################################
        
        self.collect_amp_start_win = QtGui.QDoubleSpinBox()
        self.collect_amp_start_win.setRange(self.time_array.min(), self.time_array.max())
        self.collect_amp_start_win.setValue(0.01)
        self.collect_amp_start_win.setDecimals(4)
        self.collect_amp_start_win.setSingleStep(0.001)
        collect_amp_start_win_label = QtGui.QLabel('P-P start')
        
        self.collect_amp_stop_win = QtGui.QDoubleSpinBox()
        self.collect_amp_stop_win.setRange(self.time_array.min(), self.time_array.max())
        self.collect_amp_stop_win.setValue(0.045)
        self.collect_amp_stop_win.setDecimals(4)
        self.collect_amp_stop_win.setSingleStep(0.001)
        collect_amp_stop_win_label = QtGui.QLabel('P-P stop')
        
        self.collect_rms_start_win = QtGui.QDoubleSpinBox()
        self.collect_rms_start_win.setRange(self.time_array.min(), self.time_array.max())
        self.collect_rms_start_win.setValue(-0.05)
        self.collect_rms_start_win.setDecimals(4)
        self.collect_rms_start_win.setSingleStep(0.001)
        collect_rms_start_win_label = QtGui.QLabel('RMS window start')
        
        self.collect_rms_stop_win = QtGui.QDoubleSpinBox()
        self.collect_rms_stop_win.setRange(self.time_array.min(), self.time_array.max())
        self.collect_rms_stop_win.setValue(0)
        self.collect_rms_stop_win.setDecimals(4)
        self.collect_rms_stop_win.setSingleStep(0.001)
        collect_rms_stop_win_label = QtGui.QLabel('RMS window stop')
             
        rms_layout = QtGui.QGridLayout()
        hotspot_layout = QtGui.QGridLayout()
        collect_layout = QtGui.QGridLayout()
             
        rms_check = QtGui.QWidget()
        hotspot = QtGui.QWidget()
        collect = QtGui.QWidget()

        pg.setConfigOptions(antialias=True)
        
        self.rms_axes = np.empty((rows, cols), dtype=object)
        self.collect_axes = np.empty((rows, cols), dtype=object)
        self.rms = np.empty((rows, cols), dtype=object)
        self.collect_amp = np.empty((rows, cols), dtype=object)
        
        for i in xrange(rows):
            for j in xrange(cols):
                rms_ax = rms_win.addPlot(title="EMG%d" % (i * cols + j))
                #ax.disableAutoRange(axis=None)
                self.rms_axes[i,j]=rms_ax.plot(np.random.normal(1,1, size=1000))
                rms_ax.setLabel('left', text='Amplitude', units='V', color='#FFFFFF', size='10pt')
                rms_ax.setLabel('bottom', text='Time', unit= 's', color='#FFFFFF', size='10pt') 
                self.rms[i,j] = rms_win.addLabel(title = str(0.00), color='#FFFFFF', size='13pt')
            rms_win.nextRow()
        
        for i in xrange(self.config.nemg):
            pick_chan.addItem('Channel%d' % (i)) #creates combobox channels to select
        
        self.hotspot_ax = self.hotspot_win.addPlot(title="EMG")
        self.hotspot_plot =  self.hotspot_ax.plot(x=self.time_array, y=np.random.normal(1,1, size=self.npt))
        self.hotspot_ax.setLabel('left', text='Amplitude', units='V', color='#FFFFFF', size='10pt')
        self.hotspot_ax.setLabel('bottom', text='Time', unit= 's', color='#FFFFFF', size='10pt')
        self.hotspot_amp = self.hotspot_win.addLabel(title = str(0.00), color='#FFFFFF', size='13pt')
        self.chan_index = 0

        for i in  xrange(rows):
            for j in xrange(cols):
                collect_ax = collect_win.addPlot(title="EMG%d" % (i * cols + j))
                self.collect_axes[i,j] = collect_ax.plot(x=self.time_array, y=np.random.normal(1,1, size=self.npt))
                collect_ax.setLabel('left', text='Amplitude', units='V', color='#FFFFFF', size='10pt')
                collect_ax.setLabel('bottom', text='Time', unit= 's', color='#FFFFFF', size='10pt') 
                self.collect_amp[i,j] = collect_win.addLabel(title = str(0.00), color='#FFFFFF', size='13pt')
            collect_win.nextRow()

        #RMS tab Layout
        rms_layout.addWidget(rms_win, 0, 0)
        
        ### Hotspot Tab ###
        hotspot_layout.addWidget(self.hotspot_win, 2, 0, 1, 10)
        hotspot_layout.addWidget(pick_chan, 0, 8)
        hotspot_layout.addWidget(hotspot_amp_start_win_label, 0, 0)
        hotspot_layout.addWidget(self.hotspot_amp_start_win, 0, 1)
        hotspot_layout.addWidget(hotspot_amp_stop_win_label, 0, 2)
        hotspot_layout.addWidget(self.hotspot_amp_stop_win, 0, 3)
        hotspot_layout.addWidget(hotspot_rms_start_win_label , 0, 4)
        hotspot_layout.addWidget(self.hotspot_rms_start_win, 0, 5)
        hotspot_layout.addWidget(hotspot_rms_stop_win_label, 0, 6)
        hotspot_layout.addWidget(self.hotspot_rms_stop_win, 0, 7)
        ######################################################################
        
        ### Collect Tab ###
        collect_layout.addWidget(collect_win,1, 0, 1, 10)
        collect_layout.addWidget(save_box, self.collect_axes.shape[0]+2, 0)
        collect_layout.addWidget(collect_amp_start_win_label, 0, 0)
        collect_layout.addWidget(self.collect_amp_start_win, 0, 1)
        collect_layout.addWidget(collect_amp_stop_win_label , 0, 2)
        collect_layout.addWidget(self.collect_amp_stop_win , 0, 3)
        collect_layout.addWidget(collect_rms_start_win_label , 0, 4)
        collect_layout.addWidget(self.collect_rms_start_win , 0, 5)
        collect_layout.addWidget(collect_rms_stop_win_label , 0, 6)
        collect_layout.addWidget(self.collect_rms_stop_win , 0, 7)
        #######################################################################
        #Assign layout to tabs
        rms_check.setLayout(rms_layout)
        hotspot.setLayout(hotspot_layout)
        collect.setLayout(collect_layout)
        
        tabs.addTab(rms_check, 'RMS Check')
        tabs.addTab(hotspot, 'Hotspot and Threshold')
        tabs.addTab(collect, 'Collect Data')
        #pdb.set_trace()
        tabs.currentChanged.connect(self.currentTab)
        self.setWindowTitle('EMG/ TMS Fast display')
      
        self.setCentralWidget(tabs)

        timer = QtCore.QTimer(self)
        timer.connect(timer, QtCore.SIGNAL("timeout()"), self.timer_event)
        timer.start(0)

    def change_chan(self, i):
        self.chan_index = i

    def update_plot(self):
        for i in xrange(self.config.nemg):
            self.new_data = self.bgprocess.old_data
            if self.current_tab == 'RMS':
                if i<= self.config.nemg:
                    self.rms_axes.flat[i].setData(self.new_data[i]/self.amplification[i])
                    self.rms.flat[i].setText(text=str(np.sqrt(np.mean(np.square(self.new_data[i]/self.amplification[i]*1000))))[:6])
        if self.current_tab == 'Hotspot':
            if self.bgprocess.new_hotspot_data:
                self.hotspot_data = self.bgprocess.hotspot_data
                self.hotspot_plot.setData(x=self.time_array, y=self.hotspot_data[self.chan_index]/self.amplification[self.chan_index]) 
                self.hotspot_amp.setText(
                    text=('Amplitude: '
                        + str(np.amax(self.hotspot_data[self.chan_index, (self.hotspot_amp_start_win.value()+self.config.pre_trig)*self.config.nsamp: \
                        (self.hotspot_amp_stop_win.value()+self.config.pre_trig)*self.config.nsamp]
                                    / self.amplification[self.chan_index]*1000)
                            - np.amin(self.hotspot_data[self.chan_index, (self.hotspot_amp_start_win.value()+self.config.pre_trig)*self.config.nsamp: \
                            (self.hotspot_amp_stop_win.value()+self.config.pre_trig)*self.config.nsamp]
                                    / self.amplification[self.chan_index]*1000))[:6] 
                        + 'mV \r\n RMS: ' 
                        + str(np.sqrt(np.mean(np.square(self.hotspot_data[self.chan_index, (self.hotspot_rms_start_win.value()+self.config.pre_trig)*self.config.nsamp: \
                        (self.hotspot_rms_stop_win.value()+self.config.pre_trig)*self.config.nsamp]
                                / self.amplification[self.chan_index]*1000))))[:6] 
                        + 'mV')
                )
                self.bgprocess.new_hotspot_data = False
        if self.current_tab == 'Collect':
            if self.bgprocess.new_collect_data:
                
                self.collect_data = self.bgprocess.collect_data
                for i in xrange(self.config.nemg):
                
                    self.collect_axes.flat[i].setData(x=self.time_array, y=self.collect_data[i]/self.amplification[i]) 
                    self.collect_amp.flat[i].setText(
                        text=('Amplitude: ' 
                            + str(np.amax(self.collect_data[i, (self.collect_amp_start_win.value()++self.config.pre_trig)*self.config.nsamp: \
                            (self.collect_amp_stop_win.value()+self.config.pre_trig)*self.config.nsamp]
                                    / self.amplification[i]*1000)
                                - np.amin(self.collect_data[i, (self.collect_amp_start_win.value()+self.config.pre_trig)*self.config.nsamp: \
                                (self.collect_amp_stop_win.value()+self.config.pre_trig)*self.config.nsamp]
                                    / self.amplification[i]*1000))[:6] 
                            + 'mV \r\n ' + 'RMS: ' 
                            + str(np.sqrt(np.mean(np.square(self.collect_data[i, (self.collect_rms_start_win.value()+self.config.pre_trig)*self.config.nsamp: \
                            (self.collect_rms_stop_win.value()+self.config.pre_trig)*self.config.nsamp]
                                   / self.amplification[i]*1000))))[:6] 
                            + 'mV')
                    )
                    if self.saving_data:
                        save_array = np.insert(self.collect_data[i]/self.amplification[i], 0, i, axis=0)
                        
                        f = open(self.filename, 'a')
                        np.savetxt(f, np.insert(save_array, 0, [self.trial_count, self.bgprocess.trial_count])[np.newaxis], fmt='%f', delimiter=',', newline='\r\n')
                        f.close()
                self.trial_count += 1
            self.bgprocess.new_collect_data = False
  
    def timer_event(self):
        done = False
        #print self.bgprocess.trial_count
        #sys.stdout.write("~")
        sys.stdout.flush()
        self.update_plot()
        '''
        while not done:
            msg = CMessage()
            rcv = self.mod.ReadMessage(msg, 0)
            if rcv == 1:
                    self.process_message(msg)
            else:
                done = True
        '''
    def create_menu(self):
        self.menu = QtGui.QMenuBar()
        self.amp_menu = self.menu.addMenu('Amplification')
        self.chan_setting = np.empty((self.config.nemg), dtype=object)
        self.chan_groups = np.empty((self.config.nemg), dtype=object)
        for i in xrange(self.config.nemg):
            self.chan_groups[i], act1, act2 = self.create_action(i)
            chan = self.amp_menu.addMenu('Channel %d' %i)
            chan.addAction(act1)
            chan.addAction(act2)
            self.chan_setting[i] = chan

    def create_action(self,i):
        amp_group = QtGui.QActionGroup(self)
        five = QtGui.QAction(amp_group, checkable=True)
        five.setText('500')
        ten = QtGui.QAction(amp_group, checkable=True, checked=True)
        ten.setText('1000')
        amp_group.triggered.connect(lambda button: self.change_amp(button, i))
        return amp_group, five, ten
    
    def change_amp(self, button, i):
        if button.text() == '500':
            if self.amplification[i] == 500:
                pass
            else:
                self.amplification[i] = 500
                print 'changed to 500'
        if button.text() == '1000':
            if self.amplification[i] == 1000:
                pass
            else:
                self.amplification[i] = 1000
                print 'changed to 1000'
                
    def save_data(self, btn):
        if btn == 2:
            self.saving_data = True
            self.filename = '.'.join(((self.filename.split('.'))[0], str(self.file_number), self.filename.split('.')[-1]))
            print self.filename
        else:
             self.saving_data = False
             self.file_number += 1
             self.trial_count = 0
    
    def currentTab(self, i):
        if i == 0:
            self.current_tab = 'RMS'
        elif i == 1:
            self.current_tab = 'Hotspot'
        elif i == 2:
            self.current_tab = 'Collect'

    def _update_when_triggered(self):
        while self.current_tab == 'Hotspot':
            self.bgprocess.external_trigger()
            time.sleep(2)
            
            
            
            
class BackgroundProcess(threading.Thread):
    
    def __init__(self, parent, config_file, server):
        threading.Thread.__init__(self)
        self.load_config(config_file)
        self.parent = parent
        self.collect_count = 0
        self.hotspot_count = 0
        self.trial_count = 0
        self.TMS_trigger = False
        self.new_hotspot_data = False  
        self.new_collect_data = False
        self.ext_trig = False
        self.daemon = True
        self.setup_dragonfly(server)
        self.ext_trig = dit.DAQOut()
        self.setup_buffers()
        self.start()
     
    def load_config(self, config_file):
        cfg = SafeConfigParser()
        cfg.read(config_file)
        self.config = Config()
        #daq_config.minV  = cfg.getfloat('main', 'minV')
        #daq_config.maxV  = cfg.getfloat('main', 'maxV')
        self.config.nsamp = cfg.getint('main', 'nsamp_per_chan_per_second')
        self.config.nchan = cfg.getint('main', 'nchan')
        self.config.nemg = cfg.getint('main', 'nemg')
        self.config.nirq  = self.freq = cfg.getint('main', 'nirq_per_second')
        self.config.pre_trig = cfg.getfloat('main', 'pre_trigger')
        self.config.trig_chan = cfg.getfloat('main', 'trig_chan')
        self.config.perchan = self.config.nsamp / self.config.nirq
        self.config.npt   = self.config.nsamp * self.config.nchan / self.config.nirq
        self.config.pre_trig_samp = self.config.pre_trig * self.config.nsamp
        assert((self.config.nsamp * self.config.nchan) % self.config.nirq == 0)
        assert(self.config.nsamp % self.config.nirq == 0)
        
    
    def setup_dragonfly(self, server):
        subscriptions = [MT_EXIT, \
                         rc.MT_PING, \
                         rc.MT_DAQ_DATA, \
                         rc.MT_SAMPLE_GENERATED, \
                         rc.MT_TMS_TRIGGER, \
                         rc.MT_MNOME_STATE]
        self.mod = Dragonfly_Module(0, 0)
        self.mod.ConnectToMMM(server)
        for sub in subscriptions:
            self.mod.Subscribe(sub)
        self.mod.SendModuleReady()
        print "Connected to Dragonfly at ", server
        
    def setup_buffers(self):
        self.cols = 2
        self.rows = self.config.nemg / self.cols
        self.npt = 2500 
        self.old_data = np.zeros((self.config.nchan, self.npt+300)) #  bigger to account for variable trig location in buffer
        self.new_data = np.zeros((self.config.nchan, self.npt+300))
        self.collect_data = np.zeros((self.config.nemg, self.npt))
        self.hotspot_data = np.zeros((self.config.nemg, self.npt))
        
        
    def run(self):
        while True:
            msg = CMessage()
            rcv = self.mod.ReadMessage(msg, 0)
            if rcv == 1:
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
                    respond_to_ping(self.mod, msg, 'fast_display_rms')    
                else:
                    self.process_message(msg)
                

    def process_message(self, msg):
        msg_type = msg.GetHeader().msg_type
        dest_mod_id = msg.GetHeader().dest_mod_id
        if msg_type == rc.MT_TMS_TRIGGER:
            self.trial_count += 1
            self.ext_trig.run()
            self.TMS_trigger = True
        if msg_type == rc.MT_MNOME_STATE:
            mdf = rc.MDF_MNOME_STATE()
            copy_from_msg(mdf, msg)
            if mdf.state == 0:
                self.trial_count = 0
                print 'reset count'
            else:
                pass
        else:
            # if it is a NiDAQ message from channels 0-7, plot the data
            #self.counter += 1
            if msg_type == rc.MT_DAQ_DATA:
                #sys.stdout.write("*")
                #sys.stdout.flush()
                mdf = rc.MDF_DAQ_DATA()
                copy_from_msg(mdf, msg)
                # add data to data buffers (necessary, or just use graphics buffers?)
                # update plots to new data buffers
                buf = mdf.buffer
 
                self.new_data[:,:-self.config.perchan] = self.old_data[:,self.config.perchan:]
                for i in xrange(self.config.nchan):
                    #if i == 0:
                    #    print mdf.buffer[perchan * i:perchan * (i + 1)].size
                    self.new_data[i, -self.config.perchan:] = buf[i:self.config.nchan * self.config.perchan:self.config.nchan]
                self.old_data[:] = self.new_data[:]
                
        if self.parent.current_tab ==  'Collect':
            if self.TMS_trigger:
                if self.config.pre_trig_samp <= np.argmax(self.old_data[self.config.trig_chan, :] >= 3) <= (self.config.pre_trig_samp)+200:
                    self.trig_index = np.argmax(self.old_data[self.config.trig_chan, :] >= 3)
                    self.collect_data = self.old_data[:self.config.nemg, self.trig_index - self.config.pre_trig_samp:self.trig_index + self.npt - self.config.pre_trig_samp]
                    self.old_data = self.old_data * 0
                    self.new_data = self.new_data * 0
                    self.new_collect_data = True
                    self.TMS_trigger = False
                    #print self.trial_count
                
        if self.parent.current_tab == 'Hotspot':
            if self.config.pre_trig_samp <= np.argmax(self.old_data[self.config.trig_chan, :] >= 3) <= self.config.pre_trig_samp+200:
                self.trig_index = np.argmax(self.old_data[self.config.trig_chan, :] >= 3)
                self.hotspot_data = self.old_data[:self.config.nemg, self.trig_index - self.config.pre_trig_samp:self.trig_index + self.npt - self.config.pre_trig_samp]
                self.new_hotspot_data = True
                self.old_data = self.old_data * 0
                self.new_data = self.new_data * 0
    

    
if __name__ == "__main__":
    parser = ArgumentParser(description = 'Real-time display of 8-channel EMG')
    parser.add_argument(type=str, dest='config')
    parser.add_argument(type=str, dest='mm_ip', nargs='?', default='127.0.0.1:7111')
    args = parser.parse_args()
    print("Using config file=%s, MM IP=%s" % (args.config, args.mm_ip))
    app = QtGui.QApplication(sys.argv)
    fd = FastDisplay(args.config, args.mm_ip)
    fd.show()
    exit_status = app.exec_()
    sys.exit(exit_status)
    print "Got here"