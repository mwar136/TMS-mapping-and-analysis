#Author: Matthew Ward mwar136@aucklanduni.ac.nz

import numpy as np
import Dragonfly_config as rc
from argparse import ArgumentParser
from ConfigParser import SafeConfigParser
from PyDragonfly import Dragonfly_Module, CMessage, copy_to_msg, copy_from_msg, MT_EXIT
from dragonfly_utils import respond_to_ping
import winsound
import sys
import threading

class Metronome(object):
    
    def __init__(self, config_file, mm_ip):
        self.load_config(config_file)
        self.count = 0
        self.pause_state = True
        self.setup_Dragonfly(mm_ip)
        self.calc_rates()
        self.run()
        
    def load_config(self, config_file):
        self.config = SafeConfigParser()
        self.config.read(config_file)
        self.pretrigger_time = self.config.getfloat('metronome', 'pretrigger time')
        self.metronome_period = self.config.getfloat('metronome', 'metronome period')
        self.in_msg_type = 'DAQ_DATA' # trigger msg
        self.in_msg_num = eval('rc.MT_%s' % (self.in_msg_type.upper()))
        print self.in_msg_num, 'config load complete'
        
    def calc_rates(self):
        self.in_msg_freq = 1 / self.chk_msg()
        self.metronome_count = self.metronome_period * self.in_msg_freq
        if self.pretrigger_time > 0: #negative pre-trigger fire after metronome
            self.trigger_out_count = self.metronome_count - self.pretrigger_time * self.in_msg_freq
        else:    
            self.trigger_out_count = self.metronome_count + self.pretrigger_time * self.in_msg_freq
        print 'Got frequency! %d' %self.in_msg_freq
        print  self.metronome_count, self.trigger_out_count

    def chk_msg(self):
        while True:
            in_msg = CMessage()
            rcv = self.mod.ReadMessage(in_msg, 0.1)
            if rcv == 1:
                msg_type = in_msg.GetHeader().msg_type
                if  msg_type == MT_EXIT:
                    if (dest_mod_id == 0) or (dest_mod_id == self.mod.GetModuleID()):
                        print 'Received MT_EXIT, disconnecting...'
                        self.mod.SendSignal(rc.MT_EXIT_ACK)
                        self.mod.DisconnectFromMMM()
                        break
                elif msg_type == rc.MT_PING:
                    respond_to_ping(self.mod, in_msg, 'Metronome')
                elif msg_type == self.in_msg_num:
                    in_mdf = eval('rc.MDF_%s()' % (self.in_msg_type.upper()))
                    copy_from_msg(in_mdf, in_msg)
                    return in_mdf.sample_header.DeltaTime
                   
                
    def setup_Dragonfly(self, mm_ip):
        self.mod = Dragonfly_Module(0, 0)
        self.mod.ConnectToMMM(mm_ip)
        self.mod.Subscribe(MT_EXIT)
        self.mod.Subscribe(rc.MT_PING)
        self.mod.Subscribe(self.in_msg_num)
        self.mod.Subscribe(rc.MT_MNOME_STATE)
        self.mod.SendModuleReady()
        print "Connected to Dragonfly at", mm_ip
        
    def run(self):
        while True:
            in_msg = CMessage()
            rcv = self.mod.ReadMessage(in_msg, 0.1)
            if rcv == 1:
                msg_type = in_msg.GetHeader().msg_type
                if  msg_type == MT_EXIT:
                    if (dest_mod_id == 0) or (dest_mod_id == self.mod.GetModuleID()):
                        print 'Received MT_EXIT, disconnecting...'
                        self.mod.SendSignal(rc.MT_EXIT_ACK)
                        self.mod.DisconnectFromMMM()
                        break
                elif msg_type == rc.MT_PING:
                    respond_to_ping(self.mod, in_msg, 'Metronome')
                elif msg_type == rc.MT_MNOME_STATE:
                    print 'got message'
                    in_mdf = rc.MDF_MNOME_STATE()
                    copy_from_msg(in_mdf, in_msg)
                    if in_mdf.state == 0:
                        print 'got stop'
                        self.pause_state = True
                        self.count = 0
                    elif in_mdf.state == 1:
                        print 'got start'
                        self.pause_state = False
                        self.count = 0
                    elif in_mdf.state == 2:
                        print 'got pause'
                        self.pause_state = True
                        self.count = 0
                elif msg_type == self.in_msg_num:
                    if self.pause_state:
                        pass
                    else:
                        self.count += 1
                        if self.pretrigger_time > 0:
                            if self.count == self.metronome_count:
                                in_mdf = eval('rc.MDF_%s()' % (self.in_msg_type.upper()))
                                copy_from_msg(in_mdf, in_msg)
                                out_mdf = rc.MDF_TMS_TRIGGER()
                                out_mdf.sample_header = in_mdf.sample_header
                                out_msg = CMessage(rc.MT_TMS_TRIGGER)
                                copy_to_msg(out_mdf, out_msg)
                                self.mod.SendMessage(out_msg)
                                self.count = 0 - int(np.random.uniform(0, 1.5, 1)[0] * self.in_msg_freq)
                                    
                            if self.count == self.trigger_out_count:
                                sound_thread  = threading.Thread(target=self.play_sound)
                                sound_thread.start()

                        else:
                            if self.count == self.trigger_out_count:
                                in_mdf = eval('rc.MDF_%s()' % (self.in_msg_type.upper()))
                                copy_from_msg(in_mdf, in_msg)
                                out_mdf = rc.MDF_TMS_TRIGGER()
                                out_mdf.sample_header = in_mdf.sample_header
                                out_msg = CMessage(rc.MT_TMS_TRIGGER)
                                copy_to_msg(out_mdf, out_msg)
                                self.mod.SendMessage(out_msg)
                                    
                            if self.count == self.metronome_count:
                                self.count = 0 - int(np.random.uniform(0, 1.5, 1)[0] * self.in_msg_freq)
                                sound_thread  = threading.Thread(target=self.play_sound)
                                sound_thread.start()
                                
    
    def play_sound(self):
        winsound.Beep(1500, 1000)
    
if __name__ == "__main__":
    parser = ArgumentParser(description = "Generate metronome beep at regular interval")
    parser.add_argument('config', metavar='config_file', type=str)
    parser.add_argument(type=str, dest='mm_ip', nargs='?', default='')
    args = parser.parse_args()
    print("Using config file=%s, MM IP=%s" % (args.config, args.mm_ip))
    mw = Metronome(args.config, args.mm_ip)