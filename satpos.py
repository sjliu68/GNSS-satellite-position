# -*- coding: utf-8 -*-
"""
Created on Tue May  1 09:48:38 2018

@author: Shengjie Liu, stop68@foxmail.com

Calculate satellite position from Rinex navigation file

Test on Rinex v2.10, v3.02 with GPS navigation file

Requirements:
    numpy
    argparse
    
Example:
    python satpos.py --file=rinex302.18N
    
  Perform time correction with --timeCor=True
    python satpos.py --file=rinex210.18N --timeCor=True
  
  Use Householder's iteration instead of Newton's
    python satpos.py --file=rinex210.18N --iteration=Householder

"""

import numpy as np
import argparse

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--file', type=str, default = None)
parser.add_argument('--timeCor', type=bool, default = False)
parser.add_argument('--iteration', type=str, default = 'Newton')
args = parser.parse_args()

GM = 3.986005*np.power(10.0,14)
c = 2.99792458*np.power(10.0,8)
omegae_dot = 7.2921151467*np.power(10.0,-5)

def readRinexN(file):
    data = {}
    with open(file,'rt') as f:
        foundend = 0
        idx = 0
        for line in f:
            a = line.split(' ')
            if 'END' in a:
                foundend = 1
                continue
            if foundend==0:
                continue
            b = [x for x in a if x!='']
            try:
                b.remove('\n')
            except:
                pass
            if len(b)!=4:
                idx += 1
                data[str(idx)] = b
            else:
                for each in b:
                    data[str(idx)].append(each)
                    
    data2 = np.zeros([len(data),38])
    outercount = 0
    for (k,v) in data.items(): 
        count = 0
        for each in v:
            if 'D' in each:
                tmp = each.split('D')
                tmp = float(tmp[0])*np.power(10.0,float(tmp[1]))
            else:
                tmp = float(each)
            data2[outercount,count] = tmp
            count += 1
        outercount += 1
    return data, data2

def readRinexN302(file):
    data = {}
    with open(file,'rt') as f:
        foundend = 0
        idx = 0
        for line in f:
            a = line.split(' ')
            if 'END' in a:
                foundend = 1
                continue
            if foundend==0:
                continue
            b = [x for x in a if x!='']
            try:
                b.remove('\n')
            except:
                pass
            if len(b)!=4:
                idx += 1
                data[str(idx)] = b
            else:
                for each in b:
                    data[str(idx)].append(each)
                    
    data2 = np.zeros([len(data),38])
    outercount = 0
    for (k,v) in data.items(): 
        count = 0
        for each in v:
            if 'D' in each:
                tmp = each.split('D')
                tmp = float(tmp[0])*np.power(10.0,float(tmp[1]))
            else:
                if 'G' in each:
                    tmp = float(each.split('G')[1])
                elif 'C' in each:
                    tmp = float(each.split('C')[1])
                else:
                    tmp = float(each)
            data2[outercount,count] = tmp
            count += 1
        outercount += 1
    return data, data2

def calSatPos(data, timeCor=False, iteration='Newton'):
    sats = np.zeros([data.shape[0],5])
    for j in range(data.shape[0]):
        
        ## load variables
        A = np.power(data[j,17],2) 
        toe = data[j,18] # Time of Ephemeris
        tsv = data[j,18]
        tk = tsv - toe
        
        n0 = np.sqrt(GM/np.power(A,3)) #
        dn = data[j,12]
        n = n0 + dn
        m0 = data[j,13]
        M = m0+n*tk
        
        af0 = data[j,7]
        af1 = data[j,8]
        w = data[j,24]
        cuc = data[j,14] 
        cus = data[j,16]
        crc = data[j,23]
        crs = data[j,11]
        i0 = data[j,22]
        idot = data[j,26]
        omg0 = data[j,20]
        odot = data[j,25] 
        e = data[j,15] # Eccentricity
        
        ## time correction
        if timeCor == True:
            NRnext = 0
            NR = 1
            m = 1
            while np.abs(NRnext-NR)>np.power(10.0,-16):
                NR = NRnext
                f = NR-e*np.sin(NR)-M
                f1 = 1-e*np.cos(NR)
                f2 = e*np.sin(NR)
                if iteration=='Householder':
                    NRnext = NR - f/(f1-(f2*f/(2*f1)))
                else:
                    NRnext = NR - f/f1
                m += 1
            
            E = NRnext
            
            F = -2*np.sqrt(GM)/np.power(c,2)
            delta_tr = F*e*np.sqrt(A)*np.sin(E)
            delta_tsv = af0+af1*(tsv-toe)+delta_tr
            t = tsv-delta_tsv
            tk = t-toe
            M = m0+n*tk
        
        NRnext = 0
        NR = 1
        m = 1
        while np.abs(NRnext-NR)>np.power(10.0,-16):
            NR = NRnext
            f = NR-e*np.sin(NR)-M
            f1 = 1-e*np.cos(NR)
            f2 = e*np.sin(NR)
            if iteration=='Householder':
                NRnext = NR - f/(f1-(f2*f/(2*f1)))
            else:
                NRnext = NR - f/f1
            m += 1
    
        E = NRnext
        v = np.arctan2(np.sqrt(1-np.power(e,2))*np.sin(E),np.cos(E)-e)
        phi = v + w
        u = phi + cuc*np.cos(2*phi) + cus*np.sin(2*phi)
        r = A*(1-e*np.cos(E)) + crc*np.cos(2*phi) + crs*np.sin(2*phi)
        i = i0 + idot*tk
        Omega = omg0 + (odot-omegae_dot)*tk - omegae_dot*toe
        x1 = np.cos(u)*r
        y1 = np.sin(u)*r
        
        sats[j,0] = data[j,0]
        sats[j,1] = x1*np.cos(Omega) - y1*np.cos(i)*np.sin(Omega)
        sats[j,2] = x1*np.sin(Omega) + y1*np.cos(i)*np.cos(Omega)
        sats[j,3] = y1*np.sin(i)
    return sats

if __name__ == "__main__":
    print('\n--- Calculate satellite position ---\nN file:',args.file)
    print('Time correction =',args.timeCor,
          '\nIteration strategy =',args.iteration,'\n')
    
    rawdata,data = readRinexN302(args.file)
    satp = calSatPos(data,timeCor=args.timeCor,iteration=args.iteration)
    for each in satp:
        print('Sat:',format(np.uint8(each[0]),'2d'),each[1:-1])
    np.savetxt('satellite position.csv',satp,delimiter=',')
    print('--- Save file as "satellite position.csv" in the current directory ---')
    
    
    
    
        
