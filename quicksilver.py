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
#import tables

G = 2.959122082855911E-4

def read_aei(filename):
    df = pd.read_csv(filename, names=['time', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz'], delim_whitespace=True, skiprows=4)
    df = df.astype('float64')
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

def mod_planet_orbit(mp, ms, d, a, e, theta, peri):
    #planet starts out at pericenter
    theta = np.deg2rad(theta)
    if peri:
        v=mod_mean_mo(mp,ms,d,a)*a*((1+e)/(1-e))**.5
    else:
        v=mod_mean_mo(mp,ms,d,a)*a*((1-e)/(1+e))**.5
    vx=-np.sin(theta)*v
    vy=np.cos(theta)*v
    vp = secondary_orbit(mp, ms, d)*ms/(mp+ms)
    if peri:
        p= (1-e)*a
    else:
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

def mod_epic(mp,ms,d,a):
    k2 = (G*(-135*d**4*ms*mp*(ms**2-ms*mp+mp**2) - 48*d**2*ms*mp*(ms+mp)**2*a**2 + 64*(ms+mp)**4*a**4))/(64.*(ms+mp)**3*a**7)
    return k2**.5

def period(mp, ms, a, years=False): # in days
    T = np.pi*2*(a**3/(G*(mp+ms)))**.5
    if years:
        T = T/period(1.0,0.0,1.0)
    return T

def mod_period(mp,ms,d,a,years=False):
    p = (mod_mean_mo(mp,ms,d,a)**-1)*2*np.pi
    if years:
        p = p/period(1.0,0.0,1.0)
    return p

def sma(mp,ms,P):
    a = (G*(mp+ms)*P**2/(4*np.pi**2))**(1./3)
    return a

def mod_a_from_n(mp, ms, d, n):
    a=sma(mp,ms,2*np.pi/n)/2
    for x in xrange(16):
        r=0.0
        count = 0
        while (r < 1):
            count +=1
            a = a + 10**-x
            r =  n/mod_mean_mo(mp,ms,d,a)
            if count > 100:
                break
        a=a-10**-x
    return a

def dist_resonance(mp, ms, d, r):
    a = d*r**(2./3.)
    return a

def solve(mp, ms, d, a,ratio):
    n1 = mod_mean_mo(mp,ms,d,a)
    n2 = n1/ratio
    a2 = mod_a_from_n(mp,ms,d,n2)
    return a2

def solve_with_second(mp, ms, d,ratio):
    return mod_a_from_n(mp,ms,d,mean_mo(mp,ms,d)/ratio)

def MCO_EL2X(mu,a,e,i,p,n,l):
    # Calculates Cartesian coordinates and velocities given Keplerian orbital
    # elements (for elliptical, parabolic or hyperbolic orbits).
    #
    # Based on a routine from Levison and Duncan's SWIFT integrator.
    #
    #  mu = grav const * (central + secondary mass)
    #  q = perihelion distance
    #  e = eccentricity
    #  i = inclination                 )
    #  p = longitude of perihelion !!! )   in
    #  n = longitude of ascending node ) radians
    #  l = mean anomaly                )
    #
    #  x,y,z = Cartesian positions  ( units the same as a )
    #  u,v,w =     "     velocities ( units the same as sqrt(mu/a) )

    # Change from longitude of perihelion to argument of perihelion
    g = p - n
    #
    # Rotation factors
    si, ci = np.sin(i), np.cos(i)
    sg, cg = np.sin(g), np.cos(g)
    sn, cn = np.sin(n), np.cos(n)
    z1 = cg * cn
    z2 = cg * sn
    z3 = sg * cn
    z4 = sg * sn
    d11 =  z1 - z4*ci
    d12 =  z2 + z3*ci
    d13 = sg * si
    d21 = -z3 - z2*ci
    d22 = -z4 + z1*ci
    d23 = cg * si
    #
    # Semi-major axis
          #a = q / (1.0 - e)
    #
    # Ellipse
    romes = np.sqrt(1.0 - e*e)
    temp = MCO_KEP(e,l)
    se, ce = np.sin(temp), np.cos(temp)
    z1 = a * (ce - e)
    z2 = a * romes * se
    temp = np.sqrt(mu/a) / (1.0 - e*ce)
    z3 = -se * temp
    z4 = romes * ce * temp

    x = d11*z1 + d21*z2
    y = d12*z1 + d22*z2
    z = d13*z1 + d23*z2
    u = d11*z3 + d21*z4
    v = d12*z3 + d22*z4
    w = d13*z3 + d23*z4

    return x,y,z,u,v,w
def MCO_KEP(e,oldl):
    #Solves Kepler's equation for eccentricities less than one.
    #Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330.
    #e = eccentricity
    #l = mean anomaly      (radians)
    #u = eccentric anomaly (   "   )

    pi = np.pi
    twopi = 2.0 * pi
    piby2 = 0.5 * pi

    #Reduce mean anomaly to lie in the range 0 < l < pi
    if (oldl >= 0):
        l = oldl % twopi
    else:
        l = oldl % twopi + twopi
    sign = 1.0
    if (l > pi) :
        l = twopi - l
        sign = -1.0

    ome = 1.0 - e

    if (l >= .45) or (e <  .55) :

    # Regions A,B or C in Nijenhuis
    # -----------------------------
    #
    # Rough starting value for eccentric anomaly
        if (l < ome) :
            u1 = ome
        else:
            if (l > (pi-1.0-e)) :
                u1 = (l+e*pi)/(1.0+e)
            else:
                u1 = l + e

    # Improved value using Halley's method
        flag = u1 > piby2
        if (flag) :
            x = pi - u1
        else:
            x = u1
        
        x2 = x*x
        sn = x*(1.0 + x2*(-.16605 + x2*.00761) )
        dsn = 1.0 + x2*(-.49815 + x2*.03805)
        if (flag): dsn = -dsn
        f2 = e*sn
        f0 = u1 - f2 - l
        f1 = 1.0 - e*dsn
        u2 = u1 - f0/(f1 - .50*f0*f2/f1)
    else:

    # Region D in Nijenhuis
    # ---------------------
    #
    # Rough starting value for eccentric anomaly
        z1 = 4.0*e + .50
        p = ome / z1
        q = .50 * l / z1
        p2 = p*p
        z2 = np.exp( np.log( np.sqrt( p2*p + q*q ) + q )/1.5 )
        u1 = 2.0*q / ( z2 + p + p2/z2 )

    # Improved value using Newton's method
        z2 = u1*u1
        z3 = z2*z2
        u2 = u1 - .0750*u1*z3 / (ome + z1*z2 + .3750*z3)
        u2 = l + e*u2*( 3.0 - 4.0*u2*u2 )

    # Accurate value using 3rd-order version of Newton's method
    # N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy!

    # First get accurate values for u2 - sin(u2) and 1 - cos(u2)
    bigg = (u2 > piby2)
    if (bigg) :
        z3 = pi - u2
    else:
        z3 = u2
      
    big = (z3 > (.50*piby2))
    if (big) :
        x = piby2 - z3
    else:
        x = z3

    x2 = x*x
    ss = 1.0
    cc = 1.0

    ss = x*x2/6.*(1. - x2/20.*(1. - x2/42.*(1. - x2/72.*(1. - x2/110.*(1. - x2/156.*(1. - x2/210.*(1. - x2/272.)))))))
    cc = x2/2.*(1. -x2/12.*(1. -x2/30.*(1. -x2/56.*(1. -x2/90.*(1. -x2/132.*(1.-x2/182.*(1.-x2/240.*(1.-x2/306.))))))))

    if (big) :
        z1 = cc + z3 - 1.0
        z2 = ss + z3 + 1.0 - piby2
    else:
        z1 = ss
        z2 = cc
      

    if (bigg) :
        z1 = 2.0*u2 + z1 - pi
        z2 = 2.0 - z2
      

    f0 = l - u2*ome - e*z1
    f1 = ome + e*z2
    f2 = .50*e*(u2-z1)
    f3 = e/6.0*(1.0-z2)
    z1 = f0/f1
    z2 = f0/(f2*z1+f1)
    mco_kep = sign*( u2 + f0/((f3*z1+f2)*z2+f1) )
    
    return mco_kep
    

def elements(df,m):
    #tic()
    p1el = np.zeros((len(df),6))
    #prog = 5
    for n in xrange(len(df)):
        p1el[n,:] = xv2el_swifter(G*m.values[n],df.x.values[n],df.y.values[n],df.z.values[n],df.vx.values[n],df.vy.values[n],df.vz.values[n])
        #if int(float(n)/len(df)*100) == prog:
        #    print prog, '%'
        #    prog = prog+5
    #toc()
    el = pd.DataFrame(p1el,columns=['a', 'ecc','inc', 'pomega','capom','capm'], index=df.index)
    df_total = pd.concat([df, el], axis=1)
    return df_total

def cal_elements(folder,coord):
    # Coordinate systems:
    # central body, binary barycenter, jacobi, total barycenter
    tic()
    hdfs = glob.glob(folder+'/*.hdf')
    primary = pd.read_hdf(folder+'/STAR1.hdf','central')
    central_mass = primary.mass
    if coord == 'binarybary':
        sec = pd.read_hdf(folder+'/STAR2.hdf','binarybary')
        secondary_mass = sec.mass
        central_mass = central_mass + secondary_mass 
    for hdf in hdfs:
        df = pd.read_hdf(hdf,coord)
        el = elements(df,central_mass)
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
    s1 = pd.read_hdf(folder+'/STAR1.hdf','central')
    hdfs = glob.glob(folder+'*.hdf')
    for hdf in hdfs:
        s = pd.read_hdf(hdf,'central')
        s['a'], s['ecc'], s['inc'], s['pomega'], s['capom'], s['capm'] = xv2el_array_bound(G*s1[:len(s)].mass,s['x'],s['y'],s['z'],s['vx'],s['vy'],s['vz'])
        s.to_hdf(hdf[:-4]+'.hdf','central')

def read_param(paramfile):
    with open(paramfile) as param:
        for line in param:
            if 'mass' in line:
                mass = line[23:]
    return float(mass)

def read_clo(closefile):
    df = pd.read_csv(closefile, names=['time', 'object', 'dmin', 'a1', 'e1', 'i1', 'a2', 'e2', 'i2'], delim_whitespace=True, skiprows=4)
    return df

def rebound_process(folder):
    """Requires a 'info.txt' file with the number of bodies in the integration. 
    Also modifications to reb_output_ascii in output.c to output time and mass."""
    names2=["time","mass","x","y","z","vx","vy","vz"]
    num_bod = pd.read_csv(folder+"info.txt").NBod[0]
    xyz = pd.read_csv(folder+"output.txt",delim_whitespace=True,header=None,names=names2)
    star_count = 0
    for body in xrange(num_bod):
        print body,
        p = xyz.iloc[body::num_bod, :].reset_index(drop=True)
        p.time = p.time/np.pi/2.0
        if body == 0:
            central_mass = p.mass.median()
        if p.mass.median() > 0.007*central_mass:
            p.to_hdf(folder+"STAR"+str(body+1)+'.hdf','bary')
            star_count += 1
        else:
            p.to_hdf(folder+"PL"+str(body+1-star_count)+'.hdf','bary')

#def get_fate(hdf):
#    hdf = tables.open_file(hdf, mode='r')
#    return hdf.root._v_attrs.fate
#    hdf.close()

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

def xv2el_array_bound(mu,x,y,z,vx,vy,vz):

    mu = np.array(mu)
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    vx = np.array(vx)
    vy = np.array(vy)
    vz = np.array(vz)
    
    num = len(x)
    
    TINY=np.ones(num)*4.E-15
    #ellipse= np.ones(num)
    #parabola=np.zeros(num)
    #hyperbola= np.ones(num)
    
    a = np.zeros(num)
    ecc = np.zeros(num)
    inc = np.zeros(num)
    capom = np.zeros(num)
    omega = np.zeros(num)
    capm = np.zeros(num)
    
    u = np.zeros(num)
    w = np.zeros(num)
    cape = np.zeros(num)
    face = np.zeros(num)
    capf = np.zeros(num)
    tmpf = np.zeros(num)
    
    r = (x**2 + y**2 + z**2)**.5
    v2 = vx**2 + vy**2 + vz**2
    hx = y*vz - z*vy
    hy = z*vx - x*vz
    hz = x*vy - y*vx
    h2 = hx*hx + hy*hy + hz*hz
    h = np.sqrt(h2)
    #if np.min(h2) == 0.0: return
    rdotv = x*vx + y*vy + z*vz
    energy = 0.5*v2 - mu/r
    fac = hz/h
    inc[(fac <  -1.0)] = np.ones(num)[(fac < -1.0)]*np.pi
    inc[(fac >= -1.0) & (fac < 1.0)] = np.arccos(fac[(fac >= -1.0) & (fac < 1.0)])
    fac = np.sqrt(hx*hx + hy*hy)/h
    u = np.arctan2(y,x)
    u[(hz < 0.0)] = -u[(hz < 0.0)]
    capom[(fac >= TINY)] = np.arctan2(hx[(fac >= TINY)], -hy[(fac >= TINY)])
    u[(fac >= TINY) & (np.sin(inc)==0.)] = np.zeros(num)[(fac >= TINY) & (np.sin(inc)==0.)]
    tmp = np.arctan2(z/np.sin(inc), x*np.cos(capom) + y*np.sin(capom))
    u[(fac >= TINY) & (np.sin(inc)!=0.)] = tmp[(fac >= TINY) & (np.sin(inc)!=0.)]
    capom[(capom < 0.0)] = capom[(capom < 0.0)] + 2.*np.pi
    u[(u < 0.0)] = u[(u < 0.0)] + 2.*np.pi
    
    a = -0.5*mu/energy
    fac = 1.0 - h2/(mu*a)
    
    cape[(fac <= TINY)] = u[(fac <= TINY)]
    w[(fac <= TINY)] = u[(fac <= TINY)]
    
    ecc[(fac > TINY)] = np.sqrt(fac[(fac > TINY)])
    cape[(fac > TINY)] = np.zeros(num)[(fac > TINY)]
    tmp = (a - r)/(a*ecc)
    face[(fac > TINY)] = tmp[(fac > TINY)]
    cape[(fac > TINY) & (face <  -1.0)] = np.pi*np.ones(num)[(fac > TINY) & (face < -1.0)]
    cape[(fac > TINY) & (face >= -1.0) & (face < 1.0)] = np.arccos(face[(fac > TINY) & (face >= -1.0) & (face < 1.0)])
    cape[(fac > TINY) & (rdotv < 0.0)] = 2.*np.pi - cape[(fac > TINY) & (rdotv < 0.0)]
    fac2 = 1.0 - ecc*np.cos(cape)
    cw = (np.cos(cape) - ecc)/fac2
    sw = np.sqrt(1.0 - ecc*ecc)*np.sin(cape)/fac2
    w2 = np.arctan2(sw, cw)
    w[(fac > TINY)] = w2[(fac > TINY)]
    w[(w < 0.0)] = w[(w < 0.0)] + 2.*np.pi
    capm = cape - ecc*np.sin(cape)
    omega = u - w
    omega[(omega < 0.0)] = omega[(omega < 0.0)] + 2.*np.pi

    return a, ecc, np.degrees(inc), np.degrees(omega), np.degrees(capom), np.degrees(capm)

def mod_elements(s1,s2,df):
    theta = np.deg2rad(-df.capom)
    x1 = df.x * np.cos(theta) - df.y *np.sin(theta)
    y1 = df.x * np.sin(theta) + df.y *np.cos(theta)
    z1 = df.z
    vx1 = df.vx * np.cos(theta) - df.vy *np.sin(theta)
    vy1 = df.vx * np.sin(theta) + df.vy *np.cos(theta)
    vz1 = df.vz

    beta = np.deg2rad(-df.inc)
    df['x0'] = x1
    df['y0'] = y1 * np.cos(beta) - z1 *np.sin(beta)
    df['z0'] = y1 * np.sin(beta) + z1 *np.cos(beta)
    df['vx0'] = vx1
    df['vy0'] = vy1 * np.cos(beta) - vz1 *np.sin(beta)
    df['vz0'] = vy1 * np.sin(beta) + vz1 *np.cos(beta)
    
    df['r'] = (df.x**2 + df.y**2 + df.z**2)**.5
    df['r0'] = (df['x0']**2 + df['y0']**2 + df['z0']**2)**.5
    df['w'] = np.arctan2(df['y0'],df['x0'])
    df['w'] = np.degrees(df['w'])
    df['w'] = df['w'] % 360

    if df.a.median() <0:
        print("negative a")
    windo = period(s1.mass.iloc[0],s2.mass.iloc[0],df.a.median())/period(1,0,1)*1.25
    #dt = df.time.iloc[2]-df.time.iloc[1]
    #dt = df.time.diff()[df.time.diff()>0].median()
    dt = df.time.diff().mean()
    windo = int(windo/dt)
    if windo < 1:
        print("Output too infrequent")
    else:
        if windo < 100:
            print("Warning: Fewer than 100 outputs per orbit. Increase sampling.")
        mask = (df['r0'] == df['r0'].rolling(windo, center = True).min())
        df['mask_peri']= mask
    
        df['geo_w'] = pd.Series.copy(df['w'])
        df['geo_w'][~mask] = np.nan
        www = df['geo_w'].copy()[mask]
        wt = df.time[mask]
        mask2 = www-np.roll(www,-1) > 280
        mask3 = www-np.roll(www,-1) < -280
        wt1 = wt[mask2]
        wt2 = wt[mask3]
        for n in xrange(len(wt1)):#fixes linear interpolation problems with angles
            df['geo_w'][df.time > wt1.values[n]] = df['geo_w'][df.time > wt1.values[n]]+360
        for n in xrange(len(wt2)):
            df['geo_w'][df.time > wt2.values[n]] = df['geo_w'][df.time > wt2.values[n]]-360
    
        df['geo_w'] = df['geo_w'].interpolate(method='linear')
        df['geo_w'] = df['geo_w'] % 360
        df['geo_a'] = (df['r0'].rolling(2*windo, center = True).min() + df['r0'].rolling(2*windo, center = True).max())/2.0
        df['geo_e'] = (df['r0'].rolling(2*windo, center = True).max() - df['r0'].rolling(2*windo, center = True).min())/(df['r0'].rolling(2*windo, center = True).min() + df['r0'].rolling(2*windo, center = True).max())
            
        df["mod_M"] = np.empty((len(df)))
        df.mod_M= np.NAN
        df.mod_M[df.mask_peri] =np.arange(0,360*len(df.mod_M[df.mask_peri]),360)
        df['mod_M'] = df['mod_M'].interpolate(method='linear')
        df['mod_M'] = df['mod_M'] % 360
    
        df['true_anom'] = (df['w']-df['geo_w'])%360
    
    return df

def high_freq_pro(s1,s2,p1):
    n = len(p1.time.diff(1)[p1.time.diff(1)>p1.time.diff(1).median()*100])
    nhf = len(p1.time)/(n+1)
    df = pd.DataFrame()
    for step in xrange(n+1):
        p11 = p1[1+step*nhf:1+(step+1)*nhf]
        s11 = s1[1+step*nhf:1+(step+1)*nhf]
        s22 = s2[1+step*nhf:1+(step+1)*nhf]
        p11 = mod_elements(s11,s22,p11)
        df = pd.concat([df,p11])
    return df
        

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

def aei2hdf(folder):
    tic()
    aeis = glob.glob(folder+'*.aei')
    if len(aeis) == 0:
        print "ERROR: No .aei files found"
    else:
        maxtime = 0
        for aei in aeis:
            print aei
            p = read_aei(aei)
            if len(p) > maxtime:
                pmax = p#.dropna()
                maxtime = len(pmax)
            p.to_hdf(aei[:-4]+'.hdf','central')
        s1 = pd.DataFrame()
        s1['x'] =  pmax.x*0.0
        s1['y'] =  pmax.y*0.0
        s1['z'] =  pmax.z*0.0
        s1['vx'] = pmax.vx*0.0
        s1['vy'] = pmax.vy*0.0
        s1['vz'] = pmax.vz*0.0
        s1['time'] = pmax.time
        s1['mass'] = read_param(folder+'/param.in')
        s1.to_hdf(folder+"STAR1.hdf",'central')
    toc()
    
def fates(folder):
    # if the fate is itself, it means it was hit by others but was the more massive body
    hdfs = glob.glob(folder+'*.hdf')
    lines = [line.rstrip('\n') for line in open(folder+'info.out')]
    for hdf in hdfs:
        planet = hdf[len(folder):-4]
        planet = planet+" "
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
        #hdf = tables.open_file(hdf, mode='a')
        #hdf.root._f_setattr('fate',fate)
        #hdf.close()

def close(folder):
    clos = glob.glob(folder+'*.clo')
    for clo in clos:
        print clo
        p = read_clo(clo)
        p.to_hdf(clo[:-4]+'.hdf','close')

def binary_bary(folder, el = True, unbound=False):
    tic()
    print 'Reading Secondary orbit'
    s1 = pd.read_hdf(folder+'/STAR1.hdf','central')
    s2 = pd.read_hdf(folder+'/STAR2.hdf','central')
    print 'Calculating binary barycenter'
    xb = s2.x*s2.mass/(s1.mass+s2.mass)
    yb = s2.y*s2.mass/(s1.mass+s2.mass)
    zb = s2.z*s2.mass/(s1.mass+s2.mass)
    ub = s2.vx*s2.mass/(s1.mass+s2.mass)
    vb = s2.vy*s2.mass/(s1.mass+s2.mass)
    wb = s2.vz*s2.mass/(s1.mass+s2.mass)
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
        pb['mass'] = p.mass
        #pb = elements(pb)
        if el:
            if unbound:
                pb = elements(pb,s1.mass+s2.mass)
            else:
                pb['a'], pb['ecc'], pb['inc'], pb['pomega'], pb['capom'], pb['capm'] = xv2el_array_bound(G*(s1.mass+s2.mass),pb['x'],pb['y'],pb['z'],pb['vx'],pb['vy'],pb['vz'])
        pb.to_hdf(hdf,'binarybary')
    toc()

def jacobi(folder, el = True, unbound=False):
    tic()
    
    hdfs = glob.glob(folder+'STAR*.hdf')
    s = pd.read_hdf(folder+"STAR1.hdf",'central')
    mtot = s.mass
    mx = s.mass * s.x
    my = s.mass * s.y
    mz = s.mass * s.z
    mu = s.mass * s.vx
    mv = s.mass * s.vy
    mw = s.mass * s.vz
    
    # write jacobi for central body
    s1 = pd.DataFrame()
    s1['x'] =  s.x*0.0
    s1['y'] =  s.y*0.0
    s1['z'] =  s.z*0.0
    s1['vx'] = s.vx*0.0
    s1['vy'] = s.vy*0.0
    s1['vz'] = s.vz*0.0
    s1['time'] = s.time
    s1['mass'] = s.mass
    if el:
        if unbound:
            s1 = elements(s1,mtot)
        else:
            s1['a'], s1['ecc'], s1['inc'], s1['pomega'], s1['capom'], s1['capm'] = xv2el_array_bound(G*mtot,s['x'],s['y'],s['z'],s['vx'],s['vy'],s['vz'])
    s1.to_hdf(folder+"STAR1.hdf",'jacobi')
    
    if len(hdfs) > 1:
        print 'Binary present'
        s2 = pd.read_hdf(folder+"STAR2.hdf",'central')
        sj = pd.DataFrame()
        sj['x'] = s2.x
        sj['y'] = s2.y
        sj['z'] = s2.z
        sj['vx'] = s2.vx
        sj['vy'] = s2.vy
        sj['vz'] = s2.vz
        sj['time'] = s2.time
        sj['mass'] = s2.mass
        mtot += s2.mass
        if el:
            if unbound:
                sj = elements(sj,mtot)
            else:
                sj['a'], sj['ecc'], sj['inc'], sj['pomega'], sj['capom'], sj['capm'] = xv2el_array_bound(G*mtot,sj['x'],sj['y'],sj['z'],sj['vx'],sj['vy'],sj['vz'])
        sj.to_hdf(folder+'STAR2.hdf','jacobi')
        mx = mx + s2.mass * s2.x
        my = my + s2.mass * s2.y
        mz = mz + s2.mass * s2.z
        mu = mu + s2.mass * s2.vx
        mv = mv + s2.mass * s2.vy
        mw = mw + s2.mass * s2.vz
        
    
    hdfs = glob.glob(folder+'PL*.hdf')
    hdfs.sort()
    for hdf in hdfs:  # excluding secondary (and last planet)
        print hdf
        p = pd.read_hdf(hdf,'central')
        temp = 1.0 / mtot
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
        mtot = mtot + p.mass
        if el:
            if unbound:
                pj = elements(pj,mtot)
            else:
                pj['a'], pj['ecc'], pj['inc'], pj['pomega'], pj['capom'], pj['capm'] = xv2el_array_bound(G*mtot,pj['x'],pj['y'],pj['z'],pj['vx'],pj['vy'],pj['vz'])
        pj.to_hdf(hdf,'jacobi')

    hdfs = glob.glob(folder+'SM*.hdf')
    hdfs.sort()
    for hdf in hdfs:  # excluding secondary (and last planet)
        print hdf
        p = pd.read_hdf(hdf,'central')
        temp = 1.0 / mtot
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
        #mtot = mtot + p.mass
        if el:
            if unbound:
                pj = elements(pj,mtot)
            else:
                pj['a'], pj['ecc'], pj['inc'], pj['pomega'], pj['capom'], pj['capm'] = xv2el_array_bound(G*mtot,pj['x'],pj['y'],pj['z'],pj['vx'],pj['vy'],pj['vz'])
        pj.to_hdf(hdf,'jacobi')
    
    
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

def bary(folder, el=True, unbound=False):
    tic()
    hdfs = glob.glob(folder+'*.hdf')
    #p1 = pd.read_hdf(folder+'PL','central')
    #s1mass = float(read_param(folder+'param.in'))
    tot_mass = 0.0
    #for hdf in hdfs:
    #    print hdf
    #    p = pd.read_hdf(hdf,'central')
    #    tot_mass += p.mass
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
        xb += p.x*p.mass
        yb += p.y*p.mass
        zb += p.z*p.mass
        ub += p.vx*p.mass
        vb += p.vy*p.mass
        wb += p.vz*p.mass
        tot_mass += p.mass
    
    xb = xb/tot_mass
    yb = yb/tot_mass
    zb = zb/tot_mass
    ub = ub/tot_mass
    vb = vb/tot_mass
    wb = wb/tot_mass
    
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
        pb['mass'] = p.mass
        #pb = elements(pb)
        if el:
            if unbound:
                pb = elements(pb,tot_mass)
            else:
                pb['a'], pb['ecc'], pb['inc'], pb['pomega'], pb['capom'], pb['capm'] = xv2el_array_bound(G*tot_mass,pb['x'],pb['y'],pb['z'],pb['vx'],pb['vy'],pb['vz'])
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

def bary_to_central(folder, el = False, unbound=False):
    s1 = pd.read_hdf(folder+"STAR1.hdf",'bary')
    xc = s1.x
    yc = s1.y
    zc = s1.z
    uc = s1.vx
    vc = s1.vy
    wc = s1.vz
    hdfs = glob.glob(folder+'*.hdf')
    for hdf in hdfs:
        print hdf
        p = pd.read_hdf(hdf,'bary')
        pb = pd.DataFrame()
        pb['x'] = p.x - xc
        pb['y'] = p.y - yc
        pb['z'] = p.z - zc
        pb['vx'] = p.vx - uc
        pb['vy'] = p.vy - vc
        pb['vz'] = p.vz - wc
        pb['time'] = p.time
        pb['mass'] = p.mass
        if el:
            if unbound:
                pb = elements(pb,s1.mass)
            else:
                pb['a'], pb['ecc'], pb['inc'], pb['pomega'], pb['capom'], pb['capm'] = xv2el_array_bound(G*(s1.mass),pb['x'],pb['y'],pb['z'],pb['vx'],pb['vy'],pb['vz'])
        pb.to_hdf(hdf,'central')

class quick:
    def __init__(self, directory):
        self.directory = directory
        lines = [line.rstrip('\n') for line in open(self.directory+'param.in')]
        param = np.array([])
        for line in lines:
            if line.find("=") >0:
                param = np.append(param,line[line.find("=")+1:])
        self.param = param
        self.param_set()
        
        lines = [line.rstrip('\n') for line in open(self.directory+'mercury.inc')]
        merc_in= np.array([])
        for line in lines:
            if line.find('parameter') > 0:
                if line.find('PI') < 0:
                    if line.find("=") >0:
                        merc_in = np.append(merc_in,line[line.find("=")+1:-1])
        self.merc_in = merc_in
        self.merc_set()
        
        self.secondary_mass = .5
        self.binary_separation = .5
        self.binary_ecc = 0.0


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

    def add_big(self, m, a, e=0, kep_or_mod="Mod", name="PL", theta=0, peri=False):
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
            newrow = np.append(newrow,mod_planet_orbit(self.primary_mass,self.secondary_mass,self.binary_separation,a,e,theta,peri))
        self.planet_data = np.vstack([self.planet_data, newrow])
        
    def add_small(self, a, m=0, e=0, kep_or_mod="Mod", name="SM", theta=0, peri=False):
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
            newrow = np.append(newrow,mod_planet_orbit(self.primary_mass,self.secondary_mass,self.binary_separation,a,e,theta,peri))
        self.small_data = np.vstack([self.small_data, newrow])
    
    def add_resonance(self, m, n, e=0, kep_or_mod="Mod", name="PL", theta=0):
        # adds planet at resonance to last added planet
        a = solve(self.primary_mass,self.secondary_mass,self.binary_separation,float(self.planet_data[self.planet_num][3]),n)
        if n > 1:
            self.add_big(m, a, e, kep_or_mod, name, theta, peri=False)
        else:
            self.add_big(m, a, e, kep_or_mod, name, theta, peri=True)    

    def save_settings(self,name):#name is string
        self.collect_param()
        self.collect_merc()
        settings = np.append(self.param,self.merc_in)
        settings = np.append(settings,[self.secondary_mass,self.binary_separation,self.binary_ecc])
        np.save(name,settings)
    
    def load_settings(self,name):
        settings = np.load(name+'.npy')
        self.param = settings[:23]
        self.merc_in = settings[23:-3]
        self.secondary_mass = float(settings[-3])
        self.binary_separation = float(settings[-2])
        self.binary_ecc = float(settings[-1])
        self.param_set()
        self.merc_set()
    
    def param_set(self):
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

    def merc_set(self):
        self.max_num_bod= float(self.merc_in[0])
        self.max_num_close=float(self.merc_in[1])
        self.max_num_messages=float(self.merc_in[2])
        self.huge = self.merc_in[3]
        self.max_num_files=float(self.merc_in[4])
        self.K2 = self.merc_in[5]
        self.AU2cm = self.merc_in[6]
        self.mass_sun2grams = self.merc_in[7]
        self.isbinary = self.merc_in[8]
        self.primary_name = self.merc_in[9]
        self.ce_binary = self.merc_in[10]
    
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
                p.write(line[:line.find("=")+1]+self.param[param_count]+'\n')
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
                        m.write(line[:line.find("=")+1]+str(self.merc_in[merc_count])+')\n')
                        merc_count +=1
                else:
                    m.write(line+'\n')
            else:
                m.write(line+'\n')
        m.close()
        
        #make big.in
        lines = [line.rstrip('\n') for line in open(self.directory+'big.in')]
        b = open(self.directory+"big.in",'w')
        for line in lines[:6]:
            b.write(line+'\n')
        b.write(" STAR2 m="+str(self.secondary_mass)+'\n')
        b.write(' '+str(self.binary_separation*(1+self.binary_ecc))+" 0.0 0.0"+'\n')#start at binary apocenter
        b.write(" 0.0 "+str(secondary_orbit(self.primary_mass, self.secondary_mass, self.binary_separation)*((1-self.binary_ecc)/(1+self.binary_ecc))**.5)+" 0.0"+'\n')
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
        
# Plotting
def plot_elements(array_df, same_plot=False):
    n = len(array_df)
    if same_plot: n=1
    plt.figure(figsize=(12,2.5*n))
    m=0
    for df in array_df:
        plt.subplot(n,6,1+m)
        plt.plot(df.a)
        plt.subplot(n,6,2+m)
        plt.plot(df.ecc)
        plt.ylim([0,1])
        plt.subplot(n,6,3+m)
        plt.plot(df.inc)
        plt.ylim([0,360])
        plt.subplot(n,6,4+m)
        plt.plot(df.capm)
        plt.ylim([0,360])
        plt.subplot(n,6,5+m)
        plt.plot(df.capom)
        plt.ylim([0,360])
        plt.subplot(n,6,6+m)
        plt.plot(df.pomega)
        plt.ylim([0,360])
        if not same_plot:m += 6

def plot_mod_elements(array_df, same_plot=False):
    n = len(array_df)
    if same_plot: n=1
    plt.figure(figsize=(12,2.5*n))
    m=0
    for df in array_df:
        plt.subplot(n,6,1+m)
        plt.plot(df.geo_a)
        plt.ylim([0,df.geo_a.max()*1.05])
        plt.subplot(n,6,2+m)
        plt.plot(df.geo_e)
        plt.ylim([0,1])
        plt.subplot(n,6,3+m)
        plt.plot(df.inc)
        plt.ylim([0,360])
        plt.subplot(n,6,4+m)
        plt.plot(df.mod_M)
        plt.ylim([0,360])
        plt.subplot(n,6,5+m)
        plt.plot(df.capom)
        plt.ylim([0,360])
        plt.subplot(n,6,6+m)
        plt.plot(df.geo_w)
        plt.ylim([0,360])
        if not same_plot: m += 6