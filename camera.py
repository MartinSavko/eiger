#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy.misc import imsave

import time
import os
import sys
import numpy as np
import PyTango
import logging
import gevent

from goniometer import goniometer
import redis
from pymba import *

class camera(object):
    def __init__(self, 
                 camera_type='prosilica',
                 y_pixels_in_detector=1024, 
                 x_pixels_in_detector=1360,
                 channels=3,
                 default_exposure_time=0.05,
                 default_gain=8.,
                 pixel_format='RGB8Packed',
                 tango_address='i11-ma-cx1/ex/imag.1',
                 tango_beamposition_address='i11-ma-cx1/ex/md2-beamposition',
                 use_redis=True):

        self.y_pixels_in_detector = y_pixels_in_detector
        self.x_pixels_in_detector = x_pixels_in_detector
        self.channels = channels
        self.default_exposure_time = default_exposure_time
        self.current_exposure_time = None
        self.default_gain = default_gain
        self.current_gain = None
        self.pixel_format=pixel_format
        self.goniometer = goniometer()
        self.use_redis = use_redis
        if use_redis == True:
            self.camera = None
            self.redis = redis.StrictRedis()
        else:
            self.camera = PyTango.DeviceProxy(tango_address)
            
        self.beamposition = PyTango.DeviceProxy(tango_beamposition_address)
        self.camera_type = camera_type
        self.shape = (y_pixels_in_detector, x_pixels_in_detector, channels)
        
        #self.focus_offsets = \
           #{1: -0.0819,
            #2: -0.0903,
            #3: -0.1020,
            #4: -0.1092,
            #5: -0.1098,
            #6: -0.1165,
            #7: -0.1185,
            #8: -0.1230,
            #9: -0.1213,
            #10: -0.1230}
        # After changing zoom 2019-03-07
        self.focus_offsets = \
           {1: -0.137,
            2: -0.104,
            3: -0.114,
            4: -0.114,
            5: -0.127,
            6: -0.132,
            7: -0.125,
            8: -0.138,
            9: -0.139,
            10: -0.128}
       
        self.zoom_motor_positions = \
           {1: 34500.0,
            2: 31165.0,
            3: 27185.0,
            4: 23205.0,
            5: 19225.0,
            6: 15245.0,
            7: 11265.0,
            8: 7285.0,
            9: 3305.0,
            10: 0.0 }
        
        self.backlight = \
           {1: 10.0,
            2: 10.0,
            3: 11.0,
            4: 13.0,
            5: 15.0,
            6: 21.0,
            7: 29.0,
            8: 41.0,
            9: 50.0,
            10: 61.0}
        
        self.frontlight = \
           {1: 10.0,
            2: 10.0,
            3: 11.0,
            4: 13.0,
            5: 15.0,
            6: 21.0,
            7: 29.0,
            8: 41.0,
            9: 50.0,
            10: 61.0}
        
        self.gain = \
           {1: 8.,
            2: 8.,
            3: 8.,
            4: 8.,
            5: 8.,
            6: 8.,
            7: 8.,
            8: 8.,
            9: 8.,
            10: 8.}
        
        # calibrations zoom X2
        #self.calibrations = \
           #{1: np.array([ 0.0008041425, 0.0008059998]),
            #2: np.array([ 0.0006467445, 0.0006472503]),
            #3: np.array([ 0.0004944528, 0.0004928871]),
            #4: np.array([ 0.0003771583, 0.0003756801]),
            #5: np.array([ 0.0002871857, 0.0002864572]),
            #6: np.array([ 0.0002194856, 0.0002190059]),
            #7: np.array([ 0.0001671063, 0.0001670309]),
            #8: np.array([ 0.0001261696, 0.0001275330]),
            #9: np.array([ 0.0000966609, 0.0000974697]),
            #10: np.array([0.0000790621, 0.0000784906])}
       
        self.calibrations = \
           {1: np.array([0.00160829, 0.001612  ]),
            2: np.array([0.00129349, 0.0012945 ]),
            3: np.array([0.00098891, 0.00098577]),
            4: np.array([0.00075432, 0.00075136]),
            5: np.array([0.00057437, 0.00057291]),
            6: np.array([0.00043897, 0.00043801]),
            7: np.array([0.00033421, 0.00033406]),
            8: np.array([0.00025234, 0.00025507]),
            9: np.array([0.00019332, 0.00019494]),
            10: np.array([0.00015812, 0.00015698])}
       
       
        self.magnifications = np.array([np.mean(self.calibrations[1]/self.calibrations[k]) for k in range(1, 11)])      
        self.master = False
        
    def get_point(self):
        return self.get_image()
        
    def get_image(self, color=True):
        if color:
            return self.get_rgbimage()
        else:
            return self.get_bwimage()
    
    def get_image_id(self):
        if self.use_redis:
            image_id = self.redis.get('last_image_id')
        else:
            image_id = self.camera.imagecounter
        return image_id
        
    def get_rgbimage(self):
        if self.use_redis:
            data = self.redis.get('last_image_data')
            rgbimage = np.ndarray(buffer=data, dtype=np.uint8, shape=(1024, 1360, 3))
        else:
            rgbimage = self.camera.rgbimage.reshape((self.shape[0], self.shape[1], 3))
        return rgbimage
    
    def get_bwimage(self):
        rgbimage = self.get_rgbimage()
        return rgbimage.mean(axis=2)
    
    def save_image(self, imagename, color=True):
        if color:
            image_id, image = self.get_image_id(), self.get_rgbimage()
        else:
            image_id, image = self.get_image_id(), self.get_bwimage()
        imsave(imagename, image)
        return imagename, image, image_id
        
    def get_zoom(self):
        return self.goniometer.md2.coaxialcamerazoomvalue
    
    def set_zoom(self, value, wait=True):
        if value is not None:
            value = int(value)
            self.set_gain(self.gain[value])
            self.goniometer.md2.backlightlevel = self.backlight[value]
            self.goniometer.set_position({'Zoom': self.zoom_motor_positions[value], 'AlignmentX': self.focus_offsets[value]}, wait=wait)
            self.goniometer.md2.coaxialcamerazoomvalue = value
        
    def get_calibration(self):
        return np.array([self.get_vertical_calibration(), self.get_horizontal_calibration()])
        
    def get_vertical_calibration(self):
        return self.goniometer.md2.coaxcamscaley
        
    def get_horizontal_calibration(self):
        return self.goniometer.md2.coaxcamscalex

    def set_exposure(self, exposure=0.05):
        if not (exposure >= 3.e-6 and exposure<3):
            print('specified exposure time is out of the supported range (3e-6, 3)')
            return -1
        if not self.use_redis:
            self.camera.exposure = exposure
        if self.master:
            self.camera.ExposureTimeAbs = exposure * 1.e6
        self.redis.set('camera_exposure_time', exposure)
        self.current_exposure_time = exposure
        
    def get_exposure(self):
        if not self.use_redis:
            exposure = self.camera.exposure
        if self.master:
            exposure = self.camera.ExposureTimeAbs/1.e6
        else:
            exposure = float(self.redis.get('camera_exposure_time'))
    
    def set_exposure_time(self, exposure_time):
        self.set_exposure(exposure_time)
    
    def get_exposure_time(self):
        if not self.use_redis:
            return self.get_exposure()
    
    def get_gain(self):
        if not self.use_redis:
            gain = self.camera.gain
        elif self.master:
            gain = self.camera.GainRaw
        else:
            gain = float(self.redis.get('camera_gain'))
        return gain
    
    def set_gain(self, gain):
        if not (gain >= 0 and gain <=24):
            print('specified gain value out of the supported range (0, 24)')
            return -1
        if not self.use_redis:
            self.camera.gain = gain
        elif self.master:
            self.camera.GainRaw = int(gain)
        self.redis.set('camera_gain', gain)                
        self.current_gain = gain
        
    def get_beam_position_vertical(self):
        return self.beamposition.read_attribute('Zoom%d_Z' % self.get_zoom()).value
    
    def get_beam_position_horizontal(self):
        return self.beamposition.read_attribute('Zoom%d_X' % self.get_zoom()).value
    
    def set_frontlightlevel(self, frontlightlevel):
        self.goniometer.md2.frontlightlevel = frontlightlevel
        
    def get_frontlightlevel(self):
        return self.goniometer.md2.frontlightlevel
    
    def set_backlightlevel(self, backlightlevel):
        self.goniometer.md2.backlightlevel = backlightlevel
        
    def get_backlightlevel(self):
        return self.goniometer.md2.backlightlevel
    
    def get_width(self):
        return self.x_pixels_in_detector
    
    def get_height(self):
        return self.y_pixels_in_detector
    
    def get_image_dimensions(self):
        return [self.get_width(), self.get_height()]
    
    def run_camera(self):
        self.master = True
        
        vimba = Vimba()
        system = vimba.getSystem()
        vimba.startup()
        
        if system.GeVTLIsPresent:
            system.runFeatureCommand("GeVDiscoveryAllOnce")
            gevent.sleep(3)
        
        cameraIds = vimba.getCameraIds()
        print('cameraIds %s' % cameraIds)
        self.camera = vimba.getCamera('DEV_000F3102FD4E')
        self.camera.openCamera()
        self.camera.PixelFormat = self.pixel_format
        
        self.frame0 = self.camera.getFrame()    # creates a frame
        self.frame0.announceFrame()
        
        self.image_dimensions = (self.frame0.width, self.frame0.height)
        
        self.set_exposure(self.default_exposure_time)
        self.set_gain(self.default_gain)
        
        self.current_gain = self.get_gain()
        self.current_exposure_time = self.get_exposure_time()
        
        
        self.camera.startCapture()
        
        self.camera.runFeatureCommand("AcquisitionStart")
        
        k = 0
        last_frame_id = None
        _start = time.time()
        while self.master:
            self.frame0.waitFrameCapture()
            try:
                self.frame0.queueFrameCapture()
            except:
                print('camera: frame dropped')
                continue
            
            #img = self.frame0.getImage()
            if self.frame0._frame.frameID != last_frame_id:
                k+=1
                data = self.frame0.getBufferByteData()
                img = np.ndarray(buffer=data, 
                                 dtype=np.uint8, 
                                 shape=(self.frame0.height, self.frame0.width, self.frame0.pixel_bytes))
                
                self.redis.set('last_image_data', img.ravel().tostring())
                self.redis.set('last_image_timestamp', str(time.time()))
                self.redis.set('last_image_id', self.frame0._frame.frameID)
                self.redis.set('last_image_frame_timestamp', str(self.frame0._frame.timestamp))
                requested_gain = float(self.redis.get('camera_gain'))
                if requested_gain != self.current_gain:
                    self.set_gain(requested_gain)
                requested_exposure_time = float(self.redis.get('camera_exposure_time'))
                if requested_exposure_time != self.current_exposure_time:
                    self.set_exposure(requested_exposure_time)
                    
            if k%10 == 0:
                print('camera last frame id %d fps %.3f ' % (self.frame0._frame.frameID, k/(time.time() - _start)))
                _start = time.time()
                k = 0
            gevent.sleep(0.01)
            
        self.camera.runFeatureCommand("AcquisitionStop")
        self.close_camera()
    
    def close_camera(self):
        self.master = False
        
        with Vimba() as vimba:
            self.camera.flushCaptureQueue()
            self.camera.endCapture()
            self.camera.revokeAllFrames()
            vimba.shutdown()
    
    def start_camera(self):
        return

if __name__ == '__main__':
    cam = camera()
    cam.run_camera()