# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 11:36:20 2016
Orbit plotting
@author: adam
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob as glob
import tables

G = 2.959122082855911E-4

def read_aei(filename):
    df = pd.read_csv(filename, names=['time', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz'], delim_whitespace=True, skiprows=4)
    return df

def secondary_orbit(mp, ms, d):
    v = (G*(mp+ms)/d)**.5
    return v

def planet_orbit(mp, ms, d, a, e, theta):
    theta = np.deg2rad(theta)
    v=(G*(mp+ms)/a*(1-e)/(1+e))**.5
    vx=-np.sin(theta)*v
    vy=np.cos(theta)*v
    #print v, vx, vy, (vx**2+vy**2)**.5
    p= (1+e)*a # at apocenter
    x=p*np.cos(theta)
    y=p*np.sin(theta)
    xp = d*ms/(mp+ms)
    x = x+ xp
    vp = secondary_orbit(mp, ms, d)*ms/(mp+ms)
    vy=vy+vp
    #print 'PL'
    #print x, y, 0.0
    #print vx, vy, 0.0
    #print "0.0 0.0 0.0"
    z = 0.0
    vz = 0.0
    #return v
    return x, y, z, vx, vy, vz

def mod_planet_orbit(mp, ms, d, a, e, theta):
    theta = np.deg2rad(theta)
    v=v=mod_mean_mo(mp,ms,d,a)*a*((1-e)/(1+e))**.5
    vx=-np.sin(theta)*v
    vy=np.cos(theta)*v
    vp = secondary_orbit(mp, ms, d)*ms/(mp+ms)
    p= (1+e)*a
    x=p*np.cos(theta)
    y=p*np.sin(theta)
    xp = d*ms/(mp+ms)
    x = x+ xp
    vy=vy+vp
    z = 0.0
    vz = 0.0
    #print 'PL'
    #print x, y, 0.0
    #print vx, vy, 0.0
    #print "0.0 0.0 0.0"
    return x, y, z, vx, vy, vz

def mean_mo(mp,ms,a):
    return ( G*(mp+ms)/a**3 )**.5

def mod_mean_mo(mp,ms,d,a):
    n2 = G*( 45*d**4*mp*ms*(mp**2-mp*ms+ms**2) + 48*d**2*mp*ms*(mp+ms)**2*a**2 + 64*(mp+ms)**4*a**4)/(64 * (mp+ms)**3 * a**7)
    return n2**.5

def period(mp, ms, a): # in days
    T = np.pi*2*(a**3/(G*(mp+ms)))**.5
    return T

def sma(mp,ms,P):
    a = (G*(mp+ms)*P**2/(4*np.pi**2))**(1./3)
    return a

def dist_resonance(mp, ms, d, r):
    a = d*r**(2./3.)
    return a

def solve(mp, ms, d, a,n):
    n1 = mod_mean_mo(mp,ms,d,a)
    a=int(a)
    for x in xrange(16):
        r=1
        while (r < n):
            a = a + 10**-x
            r =  n1/mod_mean_mo(mp,ms,d,a)
        a=a-10**-x
    return a

def elements(df):
    #tic()
    p1el = np.zeros((len(df),6))
    #prog = 5
    for n in xrange(len(df)):
        p1el[n,:] = xv2el_swifter(G,df.x.values[n],df.y.values[n],df.z.values[n],df.vx.values[n],df.vy.values[n],df.vz.values[n])
        #if int(float(n)/len(df)*100) == prog:
        #    print prog, '%'
        #    prog = prog+5
    #toc()
    el = pd.DataFrame(p1el,columns=['a', 'ecc','i', 'pomega','capom','capm'], index=df.index)
    df_total = pd.concat([df, el], axis=1)
    return df_total

def cal_elements(folder,coord):
    tic()
    hdfs = glob.glob(folder+'/*.hdf')
    for hdf in hdfs:
        df = pd.read_hdf(hdf,coord)
        el = elements(df)
        el.to_hdf(hdf,coord)
    toc()

def process_all_binary(folder):
    aei2hdf(folder)
    binary_bary(folder)
    cal_elements(folder,'bary')
    jacobi(folder)
    cal_elements(folder,'jacobi')

def process_all_single(folder):
    aei2hdf(folder)
    cal_elements(folder,'central')

def read_param(paramfile):
    with open(paramfile) as param:
        for line in param:
            if 'mass' in line:
                mass = line[24:]
    return float(mass)

def get_fate(hdf):
    hdf = tables.open_file(hdf, mode='r')
    return hdf.root._v_attrs.fate
    hdf.close()

def xv2el_swifter(mu,x,y,z,vx,vy,vz):
  # data has m,x,y,z,u,v,w
  x_vec=[x,y,z]
  v_vec=[vx,vy,vz]
  # Adapted from swifter routine
  TINY=4.E-15
  ellipse=-1
  parabola=0
  hyperbola=1

  a = 0.0; ecc = 0.0; inc = 0.0; capom = 0.0; omega = 0.0; capm = 0.0
  
  r = np.sqrt(np.dot(x_vec, x_vec))
  v2 = np.dot(v_vec, v_vec)
  hx = x_vec[1]*v_vec[2] - x_vec[2]*v_vec[1]
  hy = x_vec[2]*v_vec[0] - x_vec[0]*v_vec[2]
  hz = x_vec[0]*v_vec[1] - x_vec[1]*v_vec[0]
  h2 = hx*hx + hy*hy + hz*hz
  h = np.sqrt(h2)
  if (h2 == 0.0): return
  rdotv = np.dot(x_vec,v_vec)
  energy = 0.5*v2 - mu/r
  fac = hz/h
  if (fac < -1.0): inc = np.pi
  elif (fac < 1.0): inc = np.math.acos(fac)
  fac = np.sqrt(hx*hx + hy*hy)/h
  if (fac < TINY): 
    u = np.math.atan2(x_vec[1], x_vec[0])
    if (hz < 0.0): u = -u
  else:
    capom = np.math.atan2(hx, -hy)
    if(np.sin(inc)==0.): u=0. # RAS-- to get rid of error
    else: u = np.math.atan2(x_vec[2]/np.sin(inc), x_vec[0]*np.cos(capom) + x_vec[1]*np.sin(capom))
  if (capom < 0.0): capom = capom + 2.*np.pi
  if (u < 0.0): u = u + 2.*np.pi
  if (abs(energy*r/mu) < np.sqrt(TINY)): iorbit_type = parabola
  else: 
    a = -0.5*mu/energy
    if (a < 0.0): 
      fac = -h2/(mu*a)
      if (fac > TINY): iorbit_type = hyperbola
      else: iorbit_type = parabola
    else: iorbit_type = ellipse
  if (iorbit_type == ellipse):
    fac = 1.0 - h2/(mu*a)
    if (fac > TINY):
      ecc = np.sqrt(fac)
      cape = 0.0
      face = (a - r)/(a*ecc)
      if (face < -1.0): cape = np.pi
      elif (face < 1.0): cape = np.math.acos(face)
      if (rdotv < 0.0): cape = 2.*np.pi - cape
      fac = 1.0 - ecc*np.cos(cape)
      cw = (np.cos(cape) - ecc)/fac
      sw = np.sqrt(1.0 - ecc*ecc)*np.sin(cape)/fac
      w = np.math.atan2(sw, cw)
      if (w < 0.0): w = w + 2.*np.pi
    else:
      cape = u
      w = u
    capm = cape - ecc*np.sin(cape)
  elif (iorbit_type==parabola):
    a = 0.5*h2/mu
    ecc = 1.0
    w = 0.0
    fac = 2.0*a/r - 1.0
    if (fac < -1.0): w = np.pi
    elif (fac < 1.0): w = np.math.acos(fac)
    if (rdotv < 0.0): w = 2.*np.pi - w
    tmpf = np.tan(0.5*w)
    capm = tmpf*(1.0 + tmpf*tmpf/3.0)
  elif (iorbit_type==hyperbola):
    ecc = np.sqrt(1.0 + fac)
    tmpf = (a - r)/(a*ecc)
    if (tmpf < 1.0): tmpf = 1.0
    capf = np.log(tmpf + np.sqrt(tmpf*tmpf - 1.0))
    if (rdotv < 0.0): capf = -capf
    fac = ecc*np.cosh(capf) - 1.0
    cw = (ecc - np.cosh(capf))/fac
    sw = np.sqrt(ecc*ecc - 1.0)*np.sinh(capf)/fac
    w = np.math.atan2(sw, cw)
    if (w < 0.0): w = w + 2.*np.pi
    capm = ecc*np.sinh(capf) - capf
  omega = u - w
  if (omega < 0.0): omega = omega + 2.*np.pi
  return [a, ecc, np.degrees(inc), np.degrees(omega), np.degrees(capom), np.degrees(capm)]

def tic():
    #Homemade version of matlab tic and toc functions
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time
    if 'startTime_for_tictoc' in globals():
        print "Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds."
    else:
        print "Toc: start time not set"

def single(folder):
    tic()
    aeis = glob.glob(folder+'*.aei')
    for aei in aeis:
        print aei
        p = read_aei(aei)
        p = elements(p)
        p.to_hdf(aei[:-4]+'.hdf','central')
    toc()

def aei2hdf(folder):
    tic()
    aeis = glob.glob(folder+'*.aei')
    for aei in aeis:
        print aei
        p = read_aei(aei)
        p.to_hdf(aei[:-4]+'.hdf','central')
    s1 = pd.DataFrame()
    s1['x'] =  p.x*0.0
    s1['y'] =  p.y*0.0
    s1['z'] =  p.z*0.0
    s1['vx'] = p.vx*0.0
    s1['vy'] = p.vy*0.0
    s1['vz'] = p.vz*0.0
    s1['time'] = p.time
    s1['mass'] = read_param(folder+'/param.in')
    s1.to_hdf(folder+"STAR1.hdf",'central')
    toc()
    
def fates(folder):
    # if the fate is itself, it means it was hit by others but was the more massive body
    hdfs = glob.glob(folder+'*.hdf')
    lines = [line.rstrip('\n') for line in open(folder+'info.out')]
    for hdf in hdfs:
        planet = hdf[len(folder):-4]
        fate = 'survived'
        for line in lines:
            if (line.find(planet) > 0):
                if line.find("was hit by") > 0:
                    fate = line[1:6]
                if line.find("ejected") > 0:
                    fate = "ejected"
        print planet, fate
        fates = pd.DataFrame([fate])
        fates.to_hdf(hdf,'fate')
        hdf = tables.open_file(hdf, mode='a')
        hdf.root._f_setattr('fate',fate)
        hdf.close()

def binary_bary(folder):
    tic()
    print 'Reading Secondary orbit'
    s2 = pd.read_hdf(folder+'/STAR2.hdf','central')
    s1mass = float(read_param(folder+'/param.in'))
    print 'Calculating binary barycenter'
    xb = s2.x*s2.mass[0]/(s1mass+s2.mass[0])
    yb = s2.y*s2.mass[0]/(s1mass+s2.mass[0])
    zb = s2.z*s2.mass[0]/(s1mass+s2.mass[0])
    ub = s2.vx*s2.mass[0]/(s1mass+s2.mass[0])
    vb = s2.vy*s2.mass[0]/(s1mass+s2.mass[0])
    wb = s2.vz*s2.mass[0]/(s1mass+s2.mass[0])
    hdfs = glob.glob(folder+'*.hdf')
    for hdf in hdfs:
        print hdf
        p = pd.read_hdf(hdf,'central')
        pb = pd.DataFrame()
        pb['x'] = p.x - xb
        pb['y'] = p.y - yb
        pb['z'] = p.z - zb
        pb['vx'] = p.vx - ub
        pb['vy'] = p.vy - vb
        pb['vz'] = p.vz - wb
        pb['time'] = p.time
        #pb = elements(pb)
        pb.to_hdf(hdf,'binarybary')
    toc()

def jacobi(folder):
    tic()
    s2 = pd.read_hdf(folder+'STAR2.hdf','central')
    m1 = float(read_param(folder+'param.in'))
    mtot = s2.mass
    sj = pd.DataFrame()
    sj['x'] = s2.x
    sj['y'] = s2.y
    sj['z'] = s2.z
    sj['vx'] = s2.vx
    sj['vy'] = s2.vy
    sj['vz'] = s2.vz
    sj['time'] = s2.time
    sj['mass'] = s2.mass
    sj.to_hdf(folder+'STAR2.hdf','jacobi')
    mx = s2.mass * s2.x
    my = s2.mass * s2.y
    mz = s2.mass * s2.z
    mu = s2.mass * s2.vx
    mv = s2.mass * s2.vy
    mw = s2.mass * s2.vz
    
    hdfs = glob.glob(folder+'PL*.hdf')
    for hdf in hdfs:  # excluding secondary (and last planet)
        p = pd.read_hdf(hdf,'central')
        temp = 1.0 / (mtot + m1)
        mtot = mtot + p.mass
        pj = pd.DataFrame()
        pj['x'] =  p.x  -  temp * mx
        pj['y'] =  p.y  -  temp * my
        pj['z'] =  p.z  -  temp * mz
        pj['vx'] = p.vx  -  temp * mu
        pj['vy'] = p.vy  -  temp * mv
        pj['vz'] = p.vz  -  temp * mw
        mx = mx  +  p.mass * p.x
        my = my  +  p.mass * p.y
        mz = mz  +  p.mass * p.z
        mu = mu  +  p.mass * p.vx
        mv = mv  +  p.mass * p.vy
        mw = mw  +  p.mass * p.vz
        pj['time'] = p.time
        pj['mass'] = p.mass
        pj.to_hdf(hdf,'jacobi')
    
    # write jacobi for central body
    s1 = pd.DataFrame()
    s1['x'] =  p.x*0.0
    s1['y'] =  p.y*0.0
    s1['z'] =  p.z*0.0
    s1['vx'] = p.vx*0.0
    s1['vy'] = p.vy*0.0
    s1['vz'] = p.vz*0.0
    s1['time'] = p.time
    s1['mass'] = read_param(folder+'/param.in')
    s1.to_hdf(folder+"STAR1.hdf",'jacobi')
    
#    if len(aeis) > 2: # yes the last one is really processed differently just to avoid the mx calculations
#        p = read_aei(aeis[-2])
#        temp = 1.0 / (mtot + m1)
#        mtot = mtot + p.mass
#        pj = pd.DataFrame()
#        pj['x'] =  p.x  -  temp * mx
#        pj['y'] =  p.y  -  temp * my
#        pj['z'] =  p.z  -  temp * mz
#        pj['vx'] = p.vx  -  temp * mu
#        pj['vy'] = p.vy  -  temp * mv
#        pj['vz'] = p.vz  -  temp * mw
#        pj['time'] = p.time
#        pj['mass'] = p.mass
#        pj.to_hdf(aei[-2][:-4]+'.hdf','jacobi')
#        print aeis[-2][:-4]
        
    toc()

def bary(folder):
    tic()
    hdfs = glob.glob(folder+'*.hdf')
    #p1 = pd.read_hdf(folder+'PL','central')
    s1mass = float(read_param(folder+'param.in'))
    tot_mass = s1mass
    for hdf in hdfs:
        print hdf
        p = pd.read_hdf(hdf,'central')
        tot_mass += p.mass[0]
    #print tot_mass
    
    xb = 0.0
    yb = 0.0
    zb = 0.0
    ub = 0.0
    vb = 0.0
    wb = 0.0
    
    for hdf in hdfs:
        print hdf
        p = pd.read_hdf(hdf,'central')
        xb += p.x*p.mass/tot_mass
        yb += p.y*p.mass/tot_mass
        zb += p.z*p.mass/tot_mass
        ub += p.vx*p.mass/tot_mass
        vb += p.vy*p.mass/tot_mass
        wb += p.vz*p.mass/tot_mass
    
    for hdf in hdfs:
        print hdf
        p = pd.read_hdf(hdf,'central')
        pb = pd.DataFrame()
        pb['x'] = p.x - xb
        pb['y'] = p.y - yb
        pb['z'] = p.z - zb
        pb['vx'] = p.vx - ub
        pb['vy'] = p.vy - vb
        pb['vz'] = p.vz - wb
        pb['time'] = p.time
        pb = elements(pb)
        pb.to_hdf(hdf,'totalbary')
        
    #print hdf
    #p = pd.DataFrame()
    #p['x'] =  - xb
    #p['y'] =  - yb
    #p['z'] =  - zb
    #p['vx'] = - ub
    #p['vy'] = - vb
    #p['vz'] = - wb
    #p['time'] = pb.time
    #p = elements(p)
    #p.to_hdf(folder+"STAR1.hdf",'totalbary')
    toc()


class quick:
    def __init__(self, directory):
        # 2. to refer to the inner class, you must use self.Bank
        # 3. no need to use an inner class here
        self.directory = directory
        lines = [line.rstrip('\n') for line in open('param_default.in')]
        param = np.array([])
        for line in lines:
            if line.find("=") >0:
                param = np.append(param,line[line.find("=")+2:])
        self.param = param
        self.params()
        
        lines = [line.rstrip('\n') for line in open('mercury_default.inc')]
        merc_in= np.array([])
        for line in lines:
            if line.find('parameter') > 0:
                if line.find('PI') < 0:
                    if line.find("=") >0:
                        merc_in = np.append(merc_in,line[line.find("=")+2:-1])
        self.merc_in = merc_in
        self.max_num_bod= float(merc_in[0])
        self.max_num_close=float(merc_in[1])
        self.max_num_messages=float(merc_in[2])
        self.huge = merc_in[3]
        self.max_num_files=float(merc_in[4])
        self.K2 = merc_in[5]
        self.AU2cm = merc_in[6]
        self.mass_sun2grams = merc_in[7]
        self.isbinary = merc_in[8]
        self.primary_name = merc_in[9]
        self.ce_binary = merc_in[10]
        
        self.secondary_mass = .5
        self.binary_separation = .5


        self.planet_num = 0
        self.small_num = 0
        #planet.planet_num += 1
        #self.name = name
        #planet_num = 0
        self.planet_data = np.array(["#","Name","m","a","e","x","y","z","vx","vy","vz"])
        self.small_data = np.array(["#","Name","m","a","e","x","y","z","vx","vy","vz"])
        #self.names = np.array([])
        #self.a = np.array([])
        #self.e = np.array([])
        #self.x = np.array([])
        #self.y = np.array([])
        #self.z = np.array([])

    def add_planet(self, m, a, e=0, kep_or_mod="Mod", name="PL", theta=0):
        #self.names = np.append(self.names,name)
        #self.a = np.append(self.a,a)
        #self.e = np.append(self.e,e)
        self.planet_num += 1
        if name=="PL":
            name = name+str(self.planet_num)
        newrow = [self.planet_num,name,m,a,e]
        if (kep_or_mod[0] == "k") or (kep_or_mod[0] == "K"):
            newrow = np.append(newrow,planet_orbit(self.primary_mass,self.secondary_mass,self.binary_separation,a,e,theta))
        
        if (kep_or_mod[0] == "m") or (kep_or_mod[0] == "M"):
            newrow = np.append(newrow,mod_planet_orbit(self.primary_mass,self.secondary_mass,self.binary_separation,a,e,theta))
        self.planet_data = np.vstack([self.planet_data, newrow])
        
    def add_small(self, a, m=0, e=0, kep_or_mod="Mod", name="SM", theta=0):
        #self.names = np.append(self.names,name)
        #self.a = np.append(self.a,a)
        #self.e = np.append(self.e,e)
        self.small_num += 1
        if name=="SM":
            name = name+str(self.small_num)
        newrow = [self.small_num,name,m,a,e]
        if (kep_or_mod[0] == "k") or (kep_or_mod[0] == "K"):
            newrow = np.append(newrow,planet_orbit(self.primary_mass,self.secondary_mass,self.binary_separation,a,e,theta))
        if (kep_or_mod[0] == "m") or (kep_or_mod[0] == "M"):
            newrow = np.append(newrow,mod_planet_orbit(self.primary_mass,self.secondary_mass,self.binary_separation,a,e,theta))
        self.small_data = np.vstack([self.small_data, newrow])
    
    def add_resonance(self, m, n, e=0, kep_or_mod="Mod", name="PL", theta=0):
        # adds planet at resonance to last added planet
        a = solve(self.primary_mass,self.secondary_mass,self.binary_separation,float(self.planet_data[self.planet_num][3]),n)
        self.add_planet(m, a, e, kep_or_mod, name, theta)    

    def save_param(self,name):#name is string
        self.collect_param()
        np.save(name,self.param)
    
    def load_param(self,name):
        self.param = np.load(name+'.npy')
        self.params()
        
    def params(self):
        self.algorithm = self.param[0]
        self.start_time = self.param[1]
        self.stop_time = self.param[2]
        self.output_interval = self.param[3]
        self.timestep = self.param[4]
        self.accuracy_parameter =self.param[5]
        self.stop_after_close_encounter = self.param[6]
        self.allow_collisions = self.param[7]
        self.include_collisional_fragmentation = self.param[8]
        self.express_time_in_days_or_years = self.param[9]
        self.express_time_relative_to_integration_start_time = self.param[10]
        self.output_precision = self.param[11]
        self.include_relativity = self.param[12]
        self.include_user_defined_force = self.param[13]
        self.ejection_distance = self.param[14]
        self.radius_of_central_body = self.param[15]
        self.primary_mass = float(self.param[16])
        self.central_J2 = self.param[17]
        self.central_J4 = self.param[18]
        self.central_J6 = self.param[19]
        self.Hybrid_integrator_changeover = self.param[20]
        self.data_dumps = self.param[21]
        self.periodic_effects =self.param[22]
        
    def collect_param(self):
        self.param[0] = self.algorithm
        self.param[1] = self.start_time
        self.param[2] = self.stop_time
        self.param[3] = self.output_interval
        self.param[4] = self.timestep
        self.param[5] = self.accuracy_parameter
        self.param[6] = self.stop_after_close_encounter
        self.param[7] = self.allow_collisions
        self.param[8] = self.include_collisional_fragmentation
        self.param[9] = self.express_time_in_days_or_years
        self.param[10] = self.express_time_relative_to_integration_start_time
        self.param[11] = self.output_precision
        self.param[12] = self.include_relativity
        self.param[13] = self.include_user_defined_force
        self.param[14] = self.ejection_distance
        self.param[15] = self.radius_of_central_body
        self.param[16] = self.primary_mass
        self.param[17] = self.central_J2
        self.param[18] = self.central_J4
        self.param[19] = self.central_J6
        self.param[20] = self.Hybrid_integrator_changeover
        self.param[21] = self.data_dumps
        self.param[22] = self.periodic_effects
    
    def collect_merc(self):
        self.merc_in = [self.max_num_bod,self.max_num_close,self.max_num_messages,self.huge,
                        self.max_num_files,self.K2,self.AU2cm,self.mass_sun2grams,
                        self.isbinary, self.primary_name,self.ce_binary]
    
    def build(self):
        #make param.in
        if self.isbinary == '.FALSE.':
            self.secondary_mass = 0.0
        
        self.collect_param()
        param_count=0
        lines = [line.rstrip('\n') for line in open(self.directory+'param.in')]
        p = open(self.directory+"param.in",'w')
        for line in lines:
            if line.find("=") <0:
                p.write(line+'\n')
            if line.find("=") >0:
                p.write(line[:line.find("=")+1]+" "+self.param[param_count]+'\n')
                param_count +=1
        p.close()
        
        #make mercury.inc
        self.collect_merc()
        merc_count = 0
        lines = [line.rstrip('\n') for line in open(self.directory+'mercury.inc')]
        m = open(self.directory+"mercury.inc",'w')
        for line in lines:
            if line.find('parameter') > 0:
                if line.find('PI') < 0:
                    if line.find("=") >0:
                        m.write(line[:line.find("=")+1]+" "+str(self.merc_in[merc_count])+')\n')
                        merc_count +=1
            else:
                m.write(line+'\n')
        m.close()
        
        #make big.in
        lines = [line.rstrip('\n') for line in open(self.directory+'big.in')]
        b = open(self.directory+"big.in",'w')
        for line in lines[:6]:
            b.write(line+'\n')
        b.write(" STAR2 m="+str(self.secondary_mass)+'\n')
        b.write(' '+str(self.binary_separation)+" 0.0 0.0"+'\n')
        b.write(" 0.0 "+str(secondary_orbit(self.primary_mass, self.secondary_mass, self.binary_separation))+" 0.0"+'\n')
        b.write(" 0.0 0.0 0.0"+'\n')
        for planet in xrange(self.planet_num):
            #b.write(self.planet_data[planet])
            b.write(' '+self.planet_data[planet+1][1]+" m= "+self.planet_data[planet+1][2]+'\n')
            b.write(' '+self.planet_data[planet+1][5]+' '+self.planet_data[planet+1][6]+' '+self.planet_data[planet+1][7]+'\n')
            b.write(' '+self.planet_data[planet+1][8]+' '+self.planet_data[planet+1][9]+' '+self.planet_data[planet+1][10]+'\n')
            b.write(" 0.0 0.0 0.0"+'\n')
        b.close()
        
        #make small.in
        lines = [line.rstrip('\n') for line in open(self.directory+'small.in')]
        s = open(self.directory+"small.in",'w')
        for line in lines[0:5]:
            s.write(line+'\n')
        for small in xrange(self.small_num):
            #b.write(self.small_data[small])
            s.write(' '+self.small_data[small+1][1]+" m= "+self.small_data[small+1][2]+'\n')
            s.write(' '+self.small_data[small+1][5]+' '+self.small_data[small+1][6]+' '+self.small_data[small+1][7]+'\n')
            s.write(' '+self.small_data[small+1][8]+' '+self.small_data[small+1][9]+' '+self.small_data[small+1][10]+'\n')
            s.write(" 0.0 0.0 0.0"+'\n')
        s.close()
        
    #def run(self):
        #runmer = sh.Command(folder+'mercury6')
        
        
    #def run(self):
        #str1 = 'gfortran -o '+self.directory+'mercury6 '+self.directory+'mercury6_ras.for'
        #os.system('gfortran -o mercury66 mercury6_ras.for')
        #os.system(str1)
        #str1 = 'rm '+self.directory+'*.out'
        #print str1
        #os.system(str1)
        #str1 = 'rm '+self.directory+'*.dmp'
        #print str1
        #os.system(str1)
        #str1 = 'rm '+self.directory+'*.aei'
        #print str1
        #os.system(str1)
        #str1 = self.directory+"mercury6"
        #print str1
        #os.system(str1)
        #str1 = self.directory+'element6'
        #print str1
        #os.system(str1)
        #str1 = 'rm '+self.directory+'*.clo'
        #print str1
        #os.system(str1)
        #str1 = self.directory+'close6'
        #print str1
        #os.system(str1)
        #op.process_all_binary(self.directory)
        