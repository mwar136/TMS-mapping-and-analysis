#Author: Matthew Ward mwar136@aucklanduni.ac.nz

import PyDAQmx as pdq
import threading
import numpy as np
import threading
import time

class DAQIn(pdq.Task):

    def __init__(self, parent, nsamp, nsamp_per_irq):
        pdq.Task.__init__(self)
        #self.CreateAIVoltageChan("Dev1/ai8", "", pdq.DAQmx_Val_RSE, 
        #    -5, 5, pdq.DAQmx_Val_Volts, None)
        self.CreateCICountEdgesChan("/Dev1/ctr1", "", pdq.DAQmx_Val_Rising, 0, pdq.DAQmx_Val_CountUp)
        #self.CfgDigEdgeStartTrig("/Dev1/PFI1", pdq.DAQmx_Val_Rising,)
        self.CfgSampClkTiming("/Dev1/PFI1", 10, pdq.DAQmx_Val_Rising, 
            pdq.DAQmx_Val_HWTimedSinglePoint, 1)
        #self.ActualCallback = pdq.DAQmxDoneEventCallbackPtr(self.actual_callback)
        self.AutoRegisterSignalEvent(pdq.DAQmx_Val_SampleClock, 0)
    
    def DoneCallback(self, status):
        print 'task done'
        self.StopTask()
        self.ClearTask()
       
    def SignalCallback(self):
        self.parent_callback()
        self.StopTask()
        self.ClearTask()
        
    def register_callback(self, fn):
        self.parent_callback = fn

class DAQOut(pdq.Task):

    def __init__(self):
        pdq.Task.__init__(self)
        self.CreateDOChan("Dev1/port0/line4", "", pdq.DAQmx_Val_ChanForAllLines)
        self.low = np.zeros((1), dtype=np.uint8)
        self.high = np.ones((1), dtype=np.uint8)
        self.WriteDigitalLines(1, True, 0, pdq.DAQmx_Val_GroupByChannel,
                                   self.low, None, None)
       
    def run (self):
        self.StartTask()
        self.WriteDigitalLines(1, True, 0, pdq.DAQmx_Val_GroupByChannel,
                                   self.high, None, None)
        time.sleep(.05)
        self.WriteDigitalLines(1, False, 0, pdq.DAQmx_Val_GroupByChannel,
                                   self.low, None, None)
        self.StopTask()
        
        