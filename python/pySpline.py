#!/usr/local/bin/python
'''
pySpline

Contains an relatively thin interface to the cmlib spline functions

Copyright (c) 2009 by G. Kenway
All rights reserved. Not to be used for commercial purposes.
Revision: 1.0   $Date: 24/05/2009$


Developers:
-----------
- Gaetan Kenway (GKK)

History
-------
	v. 1.0 - Initial Class Creation (GKK, 2009)
'''

__version__ = '$Revision: $'


# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, string

# =============================================================================
# External Python modules
# =============================================================================
from numpy import pi, sin, cos, linspace

# =============================================================================
# Extension modules
# =============================================================================
import pyspline as spline
import pyspline as spline_cs

# =============================================================================
# pySpline class
# =============================================================================
class spline():
	
    '''
    Spline Object Class 
    '''
    
    def __init__(self, *args, **kwargs):
		
	


	
		# Default Variable Set
		def_data = {
			"Name":'Lifting Segment',
			"Type":"external",
			"Area":1000.0,
			"Span":50.0,
			"Taper":0.4,
			"SweepLE":25.0,
			"Dihedral":5.0,
			"xc_offset":0.0,
			"root_Incidence":1.0,
			"root_Thickness":0.14,
			"root_Airfoil_type":"naca5",
			"root_Airfoil_ID":"230xx",
			"tip_Incidence":-2.0,
			"tip_Thickness":0.12,
			"tip_Airfoil_type":"naca5",
			"tip_Airfoil_ID":"230xx",
			"xrLE":0.0, "yrLE":0.0, "zrLE":0.0,
			"xRot":0.0, "yRot":0.0, "zRot":0.0,
			#'xpRot':0.0, 'ypRot':0.0, 'zpRot':0.0,
			"surface_SW_segments":5, "surface_CW_segments":20,
			"surface_SW_blending":'Equal',
			'_components':{},
			'_parent':{},
		}
		
		# Initialize Abstract Class
		GeoObject.__init__(self, def_data, *args, **kwargs)
		
		
	def _on_setitem(self,i,gobj):
		
		'''
		Sets a Component Object to Component Objects List (Component Dependent Part)
		
		Documentation last updated:  April. 14, 2008 - Ruben E. Perez
		'''
		
		# 
		if not isinstance(gobj,Effector): #or isinstance(gobj,HighLiftDevice) or isinstance(gobj,FlowControlDevice):
			try:
				gobj = Effector(gobj)
			except:
				raise TypeError("Input is not a Valid Geometry Object instance\n" \
				"Lifting Segment Components can only be Effectors\n")
			#end
		#end
		
		# 
		gobj._parent = self
		
		# 
		gobj._on_update()
		
		return gobj
		
		
	def _on_update(self):
		
		'''
		Update Geometric Attibutes (Component Dependent Part)
		
		Documentation last updated:  April. 7, 2008 - Ruben E. Perez
		'''
		
		# Inputs
		Components = self._components
		ncomp = len(Components)
		
		
		# Geometric Attributes
		self._Geometry()
		
		
		# Parametric Surface
		self._Surface()
		
		
		# Components _on_update()
		for comp in xrange(ncomp):
			Components[comp]._on_update()
		#end
		
		
	def _Geometry(self):
		
		'''
		Calculate Lifting Segment Geometric Properties
		
		Documentation last updated:  April. 7, 2008 - Ruben E. Perez
		'''
		
		# Inputs
		Segment_Type = self.Type
		Area = self.Area
		Span = self.Span
		Taper = self.Taper
		Sweep_LE = self.SweepLE*(pi/180)
		Dihedral = self.Dihedral*(pi/180)
		i_root = self.root_Incidence*(pi/180.0)
		tc_root = self.root_Thickness
		i_tip = self.tip_Incidence*(pi/180.0)
		tc_tip = self.tip_Thickness
		xrLE = self.xrLE
		yrLE = self.yrLE
		zrLE = self.zrLE
		xRot = self.xRot*(pi/180.0)
		yRot = self.yRot*(pi/180.0)
		zRot = self.zRot*(pi/180.0)
		
		
		# Aspect Ratio
		AR = Span**2/Area
		
		# Root Chord
		C_root = (2.0*Area)/(abs(Span)*(1.0+Taper))
		
		# Tip Chord
		C_tip = Taper*C_root
		
		# Mean Aerodynamic Chord
		MAC = (2.0/3.0)*((C_root+C_tip)-((C_root*C_tip)/(C_root+C_tip)))
		
		# Mean Geometric Chord
		MGC = Area/Span
		
		# Sweep Angles
		Sweep_C4 = atan(tan(Sweep_LE)-((4.0*(1.0/4.0))/(2.0*AR))*((1.0-Taper)/(1.0+Taper)))
		Sweep_C3 = atan(tan(Sweep_LE)-((4.0*(1.0/3.0))/(2.0*AR))*((1.0-Taper)/(1.0+Taper)))
		Sweep_C2 = atan(tan(Sweep_LE)-((4.0*(1.0/2.0))/(2.0*AR))*((1.0-Taper)/(1.0+Taper)))
		Sweep_TE = atan(tan(Sweep_LE)-((4.0*(1.0/1.0))/(2.0*AR))*((1.0-Taper)/(1.0+Taper)))
		
		# Wetted Area
		if (Segment_Type == "external","wingtip","winglet"):
			wetted_Area = Area*(1.977+0.52*(tc_root)*((1.0+(tc_tip/tc_root)*Taper)/(1.0+Taper)))
		else:
			wetted_Area = None
		#end
		
		# Volume
		#if (Segment_Type == "internal","external"):
		#	Volume = 
		#else:
		#	Volume = None
		##end
		
		# Planform Coordinates
		xrC4 = xrLE + (0.25*C_root)*cos(i_root)
		xrTE = xrC4 + (0.75*C_root)*cos(i_root)
		xtC4 = xrC4 + Span*tan(Sweep_C4)
		xtLE = xtC4 - (0.25*C_tip)*cos(i_tip)
		xtTE = xtC4 + (0.75*C_tip)*cos(i_tip)
		yrC4 = yrLE
		yrTE = yrLE
		ytC4 = yrLE + Span*cos(Dihedral)
		ytLE = ytC4
		ytTE = ytC4
		zrC4 = zrLE - (0.25*C_root)*sin(i_root)
		zrTE = zrC4 - (0.75*C_root)*sin(i_root)
		ztC4 = zrC4 + Span*sin(Dihedral)		
		ztLE = ztC4 + (0.25*C_tip)*sin(i_tip)
		ztTE = ztC4 - (0.75*C_tip)*sin(i_tip)
		
		#yMAC = yrLE + (1/3)*((1+2*Taper)/(1+Taper))*(Span/2)
		#xMAC = xrLE + (yMAC-yrLE)*tan(SweepLE)
		#zMAC = zrLE + (yMAC-yrLE)*tan(Dihed)
		
		
		# Output
		self.AR = AR
		self.root_Chord = C_root
		self.tip_Chord = C_tip
		self.MAC = MAC
		self.MGC = MGC
		self.SweepC4 = Sweep_C4*(180.0/pi)
		self.SweepC3 = Sweep_C3*(180.0/pi)
		self.SweepC2 = Sweep_C2*(180.0/pi)
		self.SweepTE = Sweep_TE*(180.0/pi)
		self.wetted_Area = wetted_Area
		#self.Volume = Volume
		self.xrC4 = xrC4
		self.xrTE = xrTE
		self.xtC4 = xtC4
		self.xtLE = xtLE
		self.xtTE = xtTE
		self.yrC4 = yrC4
		self.yrTE = yrTE
		self.ytC4 = ytC4
		self.ytLE = ytLE
		self.ytTE = ytTE
		self.zrC4 = zrC4
		self.zrTE = zrTE
		self.ztC4 = ztC4
		self.ztLE = ztLE
		self.ztTE = ztTE
		
		
	def _Surface(self):
		
		'''
		Calculate Lifting Segment Parametric Surface
		
		Documentation last updated:  March. 28, 2009 - Ruben E. Perez
		'''
		
		# Inputs
		Area = self.Area
		Span = self.Span
		Taper = self.Taper
		Sweep_LE = self.SweepLE*(pi/180)
		Sweep_C4 = self.SweepC4*(pi/180)
		Dihedral = self.Dihedral*(pi/180)
		C_root = self.root_Chord
		i_root = self.root_Incidence*(pi/180.0)
		tc_root = self.root_Thickness
		i_tip = self.tip_Incidence*(pi/180.0)
		tc_tip = self.tip_Thickness
		Airfoil_type_root = self.root_Airfoil_type
		Airfoil_ID_root = self.root_Airfoil_ID
		Airfoil_type_tip = self.tip_Airfoil_type
		Airfoil_ID_tip = self.tip_Airfoil_ID
		xrLE = self.xrLE
		yrLE = self.yrLE
		zrLE = self.zrLE
		#xRot = self.xRot*(pi/180.0)
		#yRot = self.yRot*(pi/180.0)
		#zRot = self.zRot*(pi/180.0)
		#xpRot = self.xpRot
		#ypRot = self.ypRot
		#zpRot = self.zpRot
		
		# Discretization Parameters
		nSW = self.surface_SW_segments	# Spanwise Discretization Parameter
		nCW = self.surface_CW_segments	# Chordwise Discretization Parameter
		bSW = self.surface_SW_blending  # Spanwise Airfoil Blending Bias 
		
		# Root Airfoil Shape
		if (Airfoil_type_root == 'datafile'):
			Airfoil_ID_r = Airfoil_ID_root
		else:
			xxID_root = string.find(Airfoil_ID_root,'xx')
			tc_r = int(round(tc_root*100))
			if (tc_r < 10):
				tcID_r = str('0') + str(tc_r)
			else:
				tcID_r = str(tc_r)
			#end
			Airfoil_ID_r = Airfoil_ID_root[:xxID_root] + tcID_r + Airfoil_ID_root[xxID_root+2:]
		#end

		Airfoil_root = Airfoil(Airfoil_type_root,Airfoil_ID_r,nCW)
		
		# Tip Airfoil Shape
		if (Airfoil_type_tip == 'datafile'):
			Airfoil_ID_t = Airfoil_ID_tip
		else:
			xxID_tip = string.find(Airfoil_ID_tip,'xx')
			tc_t = int(round(tc_tip*100))
			if (tc_t < 10):
				tcID_t = str('0') + str(tc_t)
			else:
				tcID_t = str(tc_t)
			#end
			Airfoil_ID_t = Airfoil_ID_tip[:xxID_tip] + tcID_t + Airfoil_ID_tip[xxID_tip+2:]
		#end
		Airfoil_tip = Airfoil(Airfoil_type_tip,Airfoil_ID_t,nCW)
		
		
		# Airfoil Coordinates
		np = int(nCW/2+0.5)*2 + 1
		if (Dihedral <= 30*(pi/180.0)):
			tau = 0.0
		elif (Dihedral > 60*(pi/180.0)):
			tau = pi/2.0
		elif (Dihedral > 30*(pi/180.0) and Dihedral <= 60*(pi/180.0)):
			tau = pi/4.0
		#end
		
		xf_root = Airfoil_root.shapepoints[:,0]
		yf_root = Airfoil_root.shapepoints[:,1]*sin(-tau)
		zf_root = Airfoil_root.shapepoints[:,1]*cos(-tau)
		xf_tip = Airfoil_tip.shapepoints[:,0]
		yf_tip = Airfoil_tip.shapepoints[:,1]*sin(-tau)
		zf_tip = Airfoil_tip.shapepoints[:,1]*cos(-tau)
		
		
		# Surface Coordinates
		xs = numpy.zeros([nSW+1,np], float)
		ys = numpy.zeros([nSW+1,np], float)
		zs = numpy.zeros([nSW+1,np], float)
		for i in xrange(0,nSW+1):	#for i in xrange(nSW,-1,-1):
			
			# eta_Airfoil
			xf_eta = numpy.zeros(np,float)
			yf_eta = numpy.zeros(np,float)
			zf_eta = numpy.zeros(np,float)
			for p in xrange(0,np):
				if (bSW == 'Equal'):
					eta = (float(i)/nSW)
					#eta = 0.5 - 0.5*cos((pi)*(float(i)/nSW))
				elif (bSW == 'Root'):
					eta = 1.0 - cos((pi/2.0)*(float(i)/nSW))
				elif (bSW == 'Tip'):
					eta = cos((pi/2.0) - (pi/2.0)*(float(i)/nSW))
				#end
				xf_eta[p] = (xf_tip[p] - xf_root[p])*eta + xf_root[p]
				yf_eta[p] = (yf_tip[p] - yf_root[p])*eta + yf_root[p]
				zf_eta[p] = (zf_tip[p] - zf_root[p])*eta + zf_root[p]
			#end
			
			# eta_Incidence
			i_eta = ((i_tip - i_root)/float(nSW))*float(i) + i_root
			xf_inc =  (xf_eta - 0.25)*cos(i_eta) + zf_eta*sin(i_eta)
			yf_inc = yf_eta
			zf_inc = -(xf_eta - 0.25)*sin(i_eta) + zf_eta*cos(i_eta)
			
			# Taper, Sweep, and Dihedral Correction
			for j in xrange(0,np):
				bi = (Span/float(nSW))*float(i)
				xs[i][j] = xf_inc[j]*C_root*(1.0-(bi/Span)*(1.0-Taper)) + (bi*tan(Sweep_C4))
				ys[i][j] = yf_inc[j]*C_root*(1.0-(bi/Span)*(1.0-Taper)) + (bi*cos(Dihedral))
				zs[i][j] = zf_inc[j]*C_root*(1.0-(bi/Span)*(1.0-Taper)) + (bi*sin(Dihedral))
			#end
			
		#end
		
		# Relocate Surface Coordinates
		xs = xs + xrLE + (0.25*C_root)*cos(i_root)
		ys = ys + yrLE
		zs = zs + zrLE - (0.25*C_root)*sin(i_root)
		
		# Rotate Surface Coordinates
		
		
		
		# Surface Geometric Properties (Wetted_Area and Volume)
		#if (Segment_Type == "external","wingtip","winglet"):
		#	wetted_Area =
		#else:
		#	wetted_Area = None
		##end
		#if (Segment_Type == "internal","external"):
		#	Volume = 
		#else:
		#	Volume = None
		##end
		
		
		# Output
		#self.root_Airfoil_ID = Airfoil_ID_root
		self.root_Airfoil = Airfoil_root
		#self.tip_Airfoil_ID = Airfoil_ID_tip
		self.tip_Airfoil = Airfoil_tip
		self.Surface_x = xs
		self.Surface_y = ys
		self.Surface_z = zs
		#self.Wetted_Area = wetted_Area
		#self.Volume = Volume
		
		
	def _on_plotPlanform(self,fig,update=False,*args,**kwargs):
		
		'''
		Plot Planform (Component Dependent Part)
		
		Keyword Arguments:
		------------------
		fig: OBJECT -> Parent Plot Figure
		sym: BOOLEAN -> Symmetric Surface Flag
		
		Documentation last updated:  April. 7, 2008 - Ruben E. Perez
		'''
		
		# Inputs
		xrLE = self.xrLE
		xrTE = self.xrTE
		xtLE = self.xtLE
		xtTE = self.xtTE
		yrLE = self.yrLE
		yrTE = self.yrTE
		ytLE = self.ytLE
		ytTE = self.ytTE
		zrLE = self.zrLE
		zrTE = self.zrTE
		ztLE = self.ztLE
		ztTE = self.ztTE
		
		# Parent Inputs
		if not isinstance(self._parent,GeoObject):
			sym = False
		else:
			sym = self._parent.Symmetry
		#end
		
		# Components Inputs
		ncomp = len(self._components)
		
		
		# Planform Coordinates
		x = numpy.array([xrLE, xtLE, xtTE, xrTE, xrLE])
		y = numpy.array([yrLE, ytLE, ytTE, yrTE, yrLE])
		z = numpy.array([zrLE, ztLE, ztTE, zrTE, zrLE])
		
		# 
		if not update:
			
			# Append Planform Figure Object(s)
			fig.append(PlotCurve(x,y,z,*args,**kwargs))
			if sym:
				fig.append(PlotCurve(x,-y,z,*args,**kwargs))
			#end
			
		else:
			
			# Update Surface Figure Object(s)
			fig[fig.counter].update(x,y,z,*args,**kwargs)
			fig.counter += 1
			if sym:
				fig[fig.counter].update(x,-y,z,*args,**kwargs)
				fig.counter += 1
			#end
			
		#end
		
		# Plot Components
		for comp in xrange(ncomp):
			fig = self._components[comp].plotPlanform(fig,sym,update,*args,**kwargs)
		#end
		
		return fig
		
		
	def _on_plotSurface(self,fig,update=False,*args,**kwargs):
		
		'''
		Plot Surface (Component Dependent Part)
		
		Keyword Arguments:
		------------------
		fig: OBJECT -> Parent Plot Figure
		sym: BOOLEAN -> Symmetric Surface Flag
		
		Documentation last updated:  April. 7, 2008 - Ruben E. Perez
		'''
		
		# Inputs
		x = self.Surface_x
		y = self.Surface_y
		z = self.Surface_z
		
		# Parent Inputs
		if not isinstance(self._parent,GeoObject):
			sym = False
		else:
			sym = self._parent.Symmetry
		#end
		
		# 
		if not update:
			
			# Append Surface Figure Object(s)
			fig.append(PlotSurface(x,y,z,*args,**kwargs))
			if sym:
				fig.append(PlotSurface(x,-y,z,*args,**kwargs))
			#end
			
		else:
			
			# Update Surface Figure Object(s)
			fig[fig.counter].update(x,y,z,*args,**kwargs)
			fig.counter += 1
			if sym:
				fig[fig.counter].update(x,-y,z,*args,**kwargs)
				fig.counter += 1
			#end
			
		#end
		
		return fig
		
		
	def writeOut(self,outfile,format=None):
		
		'''
		Write Formatted Output to File
		
		Keyword Arguments:
		------------------
		outfile: OBJECT -> Output File I/O Object.
		format: STRING -> Output Format [None,TeX].  Default = None
		
		Documentation last updated:  April. 7, 2008 - Ruben E. Perez
		'''
		
		# Inputs
		Components = self._components
		ncomp = len(Components)
		
		
		# Write Title Header
		outfile.write('%s - Lifting Segment\n' %(self.name))
		
		# Write Geometric Properties
		#outfile.write('\n')
		
		
		# Write Control Effectors
		outfile.write('- Control Effectors\n')
		for comp in xrange(ncomp):
			if isinstance(Components[comp],ControlEffector):
				Components[comp].writeOut(outfile,format)
			#end
		#end
		
		# Write High Lift Devices
		outfile.write('- High Lift Devices\n')
		for comp in xrange(ncomp):
			if isinstance(Components[comp],HighLiftDevice):
				Components[comp].writeOut(outfile,format)
			#end
		#end
		
		
	def writeDAT(self,datfile):
		
		'''
		Write to DAT File
		
		Keyword Arguments:
		------------------
		datfile: OBJECT   ->  DAT File I/O Object.
		
		Documentation last updated:  January. 29, 2009 - Ruben E. Perez
		'''
		
		# Inputs
		xs = self.Surface_x
		ys = self.Surface_y
		zs = self.Surface_z
		
		# Parent Inputs
		if not isinstance(self._parent,GeoObject):
			sym = False
		else:
			sym = self._parent.Symmetry
		#end
		
		
		# Surface to Patches
		fs,vs = surf2patch(xs,ys,zs,type='quad')
		nvs = vs.shape[0]
		nfs = fs.shape[0]	
		
		# Write DAT File
		datfile.write('ZONE N=%s, E=%s, F=FEPOINT, ET=QUADRILATERAL\n' %(nvs,nfs))
		for i in xrange(nvs):
			datfile.write('%s %s %s\n' %(vs[i,0],vs[i,1],vs[i,2]))
		#end
		datfile.write('\n')
		for j in xrange(nfs):
			datfile.write('%s %s %s %s\n' %(fs[j,0]+1,fs[j,1]+1,fs[j,2]+1,fs[j,3]+1))
		#end
		datfile.write('\n')
		
		# 
		if sym:
			
			# Write DAT File
			datfile.write('ZONE N=%s, E=%s, F=FEPOINT, ET=QUADRILATERAL\n' %(nvs,nfs))
			for i in xrange(nvs):
				datfile.write('%s %s %s\n' %(vs[i,0],-vs[i,1],vs[i,2]))
			#end
			datfile.write('\n')
			for j in xrange(nfs):
				datfile.write('%s %s %s %s\n' %(fs[j,0]+1,fs[j,1]+1,fs[j,2]+1,fs[j,3]+1))
			#end
			datfile.write('\n')
			
		#end
		
		
	def writeGEO(self,geofile,nSW=None,nCW=None):
		
		'''
		Write to GEO File
		
		Keyword Arguments:
		------------------
		geofile: OBJECT   ->  GEO File I/O Object.
		nSW: INTEGER: -> Spanwise Discretization.  Default = None (use Surface Discretization)
		nCW: INTEGER: -> Chordwise Discretization.  Default = None (use Surface Discretization)
		
		Documentation last updated:  October. 20, 2007 - Ruben E. Perez
		'''
		
		# Inputs
		Span = self.Span
		Taper = self.Taper
		Sweep_LE = self.SweepLE*(pi/180.0)
		Dihedral = self.Dihedral*(pi/180.0)
		C_root = self.root_Chord
		i_root = self.root_Incidence
		tc_root = self.root_Thickness
		i_tip = self.tip_Incidence
		tc_tip = self.tip_Thickness
		Airfoil_type_root = self.root_Airfoil_type
		Airfoil_ID_root = self.root_Airfoil_ID
		Airfoil_root = self.root_Airfoil
		Airfoil_type_tip = self.tip_Airfoil_type
		Airfoil_ID_tip = self.tip_Airfoil_ID
		Airfoil_tip = self.tip_Airfoil
		xrLE = self.xrLE
		yrLE = self.yrLE
		zrLE = self.zrLE
		zrC4 = self.zrC4
		
		# Discretization Parameters
		if (nSW == None):
			nSW = self.surface_SW_segments
		#end
		if (nCW == None):
			nCW = self.surface_CW_segments
		#end

			# Root Airfoil Shape
		if (Airfoil_type_root == 'datafile'):
			Airfoil_ID_r = Airfoil_ID_root
		else:
			xxID_root = string.find(Airfoil_ID_root,'xx')
			tc_r = int(round(tc_root*100))
			if (tc_r < 10):
				tcID_r = str('0') + str(tc_r)
			else:
				tcID_r = str(tc_r)
			#end
			Airfoil_ID_r = Airfoil_ID_root[:xxID_root] + tcID_r + Airfoil_ID_root[xxID_root+2:]
		#end

		Airfoil_root = Airfoil(Airfoil_type_root,Airfoil_ID_r,nCW)
		
		# Tip Airfoil Shape
		if (Airfoil_type_tip == 'datafile'):
			Airfoil_ID_t = Airfoil_ID_tip
		else:
			xxID_tip = string.find(Airfoil_ID_tip,'xx')
			tc_t = int(round(tc_tip*100))
			if (tc_t < 10):
				tcID_t = str('0') + str(tc_t)
			else:
				tcID_t = str(tc_t)
			#end
			Airfoil_ID_t = Airfoil_ID_tip[:xxID_tip] + tcID_t + Airfoil_ID_tip[xxID_tip+2:]
		#end
		Airfoil_tip = Airfoil(Airfoil_type_tip,Airfoil_ID_t,nCW)

		
		# Source Airfoil Shapes
		if (nCW != None):
			Airfoil_root = Airfoil(Airfoil_type_root,Airfoil_ID_r,nCW)
			Airfoil_tip = Airfoil(Airfoil_type_tip,Airfoil_ID_t,nCW)
		#end
		afx_root = Airfoil_root.shapepoints[:,0]
		afy_root = Airfoil_root.shapepoints[:,1]
		afx_tip = Airfoil_tip.shapepoints[:,0]
		afy_tip = Airfoil_tip.shapepoints[:,1]
		np = int(nCW/2+0.5)*2 + 1
		npU = int(np/2) + 1
		npL = npU
		
		# Airfoils Shapes Roll
		tau = Dihedral*(180.0/pi)
		
		if self._parent.Name=='Vertical Tail':
			tau = -tau
			print tau
		#endif
		#if (Dihedral <= 30*(pi/180.0)):
		#	tau = 0.0
		#elif (Dihedral > 60*(pi/180.0)):
		#	#tau = pi/2.0
		#	tau = 90.0
		#elif (Dihedral > 30*(pi/180.0) and Dihedral <= 60*(pi/180.0)):
		#	#tau = pi/4.0
		#	tau = 45.0
		##end
		
		# Undefined Parameters
		THICK = 1.0
		NEW = 1.0
		ROUND = 1.0
		YSYM = 0.0
		
		# Flap & Slats (Not Implemented Yet!)
		DSLAT = 0.0
		XSLAT = 0.0
		SWSLAT = 0.0
		DFLAP = 0.0
		XFLAP = 0.0
		SWFLAP = 0.0

		start = 1
		#print self.ListAttributes(),self.Type
		#stop
		if self.Type=='internal':
			start=0
			#print 'surfacesw',self.surface_SW_segments
			#self.surface_SW_segments+=1
		#endif
		
		# Write GEO File
		for i in xrange(start,nSW+1):
		#for i in xrange(0,nSW+1):
			
			# Spanwise Station Properties
			b_eta = (Span/float(nSW))*float(i)
			tc_eta = ((tc_tip - tc_root)/float(nSW))*float(i) + tc_root
			i_eta = ((i_tip - i_root)/float(nSW))*float(i) + i_root
			C_eta = C_root*(1.0-(b_eta/Span)*(1.0-Taper))
			xLE_eta = xrLE + b_eta*tan(Sweep_LE)
			yLE_eta = yrLE + b_eta*cos(Dihedral)
			zC4_eta = zrC4 + b_eta*sin(Dihedral)
			zLE_eta = zC4_eta + (0.25*C_eta)*sin(i_eta*(pi/180.0))
			#zLE_eta = zrLE + b_eta*sin(Dihedral)
			
			# Spanwise Station Airfoil
			afx_eta = numpy.zeros(np,float)
			afy_eta = numpy.zeros(np,float)
			for j in xrange(0,np):
				afx_eta[j] = ((afx_tip[j] - afx_root[j])/nSW)*float(i) + afx_root[j]
				afy_eta[j] = ((afy_tip[j] - afy_root[j])/nSW)*float(i) + afy_root[j]
			#end
			
				
			# Write Station Header
			geofile.write('        Z       XLE       YLE     CHORD     THICK     TWIST       NEW     ROUND  ROLL\n')
			geofile.write('%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f \n' %(yLE_eta,xLE_eta,zLE_eta,C_eta,THICK,i_eta,NEW,ROUND,-tau))

## 			if i == 3 or i==4:
## 				# Flap & Slats (Not Implemented Yet!)
## 				DSLAT = 0.0
## 				XSLAT = 0.0
## 				SWSLAT = 0.0
## 				DFLAP = 15.0#0.0
## 				XFLAP = xLE_eta+ 0.8*C_eta#0.0
## 				SWFLAP = 0.0
## 			#endif
			geofile.write(' YSYM   NU   NL      DSLAT        XSLAT     SWSLAT      DFLAP        XFLAP     SWFLAP\n')
			geofile.write(' %4.0f %4.0f %4.0f %10.5f %12.5f %10.5f %10.5f %12.5f %10.5f \n' %(YSYM,npU,npL,DSLAT,XSLAT,SWSLAT,DFLAP,XFLAP,SWFLAP))
			
			# Write Upper Surface Coordinates
			geofile.write(' Upper surface\n')
			for j in xrange(0,npU):
				geofile.write('%12.8f %12.8f\n' %(afx_eta[j],afy_eta[j]))
			#end
			
			# Write Lower Surface Coordinates
			geofile.write(' Lower surface\n')
			for j in xrange(np-1,npL-2,-1):
				geofile.write('%12.8f %12.8f\n' %(afx_eta[j],afy_eta[j]))
			#end
			
		#end
	


# =============================================================================
# LiftingSurface Class
# =============================================================================
class LiftingSurface(GeoObject):
	
	'''
	Lifting Surface Geometry Object Class - Inherited from Geometry Abstract Class
	'''
	
	def __init__(self, *args, **kwargs):
		
		'''
		LiftingSurface Object Class Initialization
		
		Keyword Arguments:
		------------------
		inp_data: DICTIONARY   ->  Names of the design variables as keys with their values as entries.
		
		Input Attributes:
		-----------------
		self.Name -> STRING: Name of the Lifting Surface.  Default = None
		self.Symmetry -> BOOLEAN: Indicates whether or not the Lifting Surface is symmetric.  Default = True
		self.xrLE -> SCALAR: LE Root x-coordinate position.  Default = 0.0
		self.yrLE -> SCALAR: LE Root y-coordinate position.  Default = 0.0
		self.zrLE -> SCALAR: LE Root z-coordinate position.  Default = 0.0
		self.xRot -> SCALAR: x-coordinate rotation.  Default = 0.0
		self.yRot -> SCALAR: y-coordinate rotation.  Default = 0.0
		self.zRot -> SCALAR: z-coordinate rotation.  Default = 0.0
		
		self._components -> DICTIONARY ->  Lifting Segment(s) Data
		
		Documentation last updated:  April. 7, 2008 - Ruben E. Perez
		'''
		
		# Default Variable Set
		def_data = {
			'Name':'Lifting Surface',
			'Symmetry':False,
			'xrLE':0.0, 'yrLE':0.0, 'zrLE':0.0,
			'xRot':0.0, 'yRot':0.0, 'zRot':0.0,
			#'xpRot':0.0, 'ypRot':0.0, 'zpRot':0.0,
			'_components':{},
			'_parent':{},
		}
		
		# Initialize Abstract Class
		GeoObject.__init__(self, def_data, *args, **kwargs)
		
		
	def _on_setitem(self,i,gobj):
		
		'''
		Sets a Component Object to Component Objects List (Component Dependent Part)
		
		Documentation last updated:  April. 7, 2008 - Ruben E. Perez
		'''
		
		# 
		if not isinstance(gobj,LiftingSegment):
			try:
				gobj = LiftingSegment(gobj)
			except:
				raise TypeError("Input is not a Valid Geometry Object instance\n" \
				"Lifting Surface Components can only be Lifting Segments\n")
			#end
		#end
		
		# 
		gobj._parent = self
		
		# 
		if (i == 0):
			gobj.xrLE = self.xrLE
			gobj.yrLE = self.yrLE
			gobj.zrLE = self.zrLE
			gobj._on_update()
		else:
			gobj.xrLE = self._components[i-1].xtLE + gobj.xc_offset*self._components[i-1].tip_Chord
			gobj.yrLE = self._components[i-1].ytLE 
			gobj.zrLE = self._components[i-1].ztLE
			gobj._on_update()
		#end
		
		return gobj
		
		
	def _on_update(self):
		
		'''
		Update Geometric Attibutes (Component Dependent Part)
		
		Documentation last updated:  April. 7, 2008 - Ruben E. Perez
		'''
		
		# Inputs
		Components = self._components
		ncomp = len(Components)
		
		# 
		if ncomp == 0:
			#print 'Skipped Update - Lifting Surface has no Segments'
			return
		#end
		
		
		# Components _on_update()
		for comp in xrange(ncomp):
			if (comp == 0):
				Components[comp].xrLE = self.xrLE
				Components[comp].yrLE = self.yrLE
				Components[comp].zrLE = self.zrLE
			else:
				Components[comp].xrLE = Components[comp-1].xtLE + Components[comp].xc_offset*Components[comp-1].tip_Chord
				Components[comp].yrLE = Components[comp-1].ytLE 
				Components[comp].zrLE = Components[comp-1].ztLE
			#end
			Components[comp]._on_update()
		#end
		
		
		# Geometric Attributes
		self._Geometry()
		
		
		# Parametric Surface
		self._Surface()
		
		
	def _Geometry(self):
		
		'''
		Calculate Geometric Properties
		
		Documentation last updated:  April. 7, 2008 - Ruben E. Perez
		'''
		
		# Inputs
		Symmetry = self.Symmetry
		Components = self._components
		ncomp = len(Components)
		
		# 
		if ncomp == 0:
			print 'Skipped Geometry Calculation - Lifting Surface has no Segments'
			return
		#end
		
		
		# Span
		Span = 0.0
		for comp in xrange(ncomp):
			Span += Components[comp].Span
		#end
		if (Symmetry == True):
			Span = 2.0*Span
		#end
		
		# Area & Exposed_Area
		Area = 0.0
		exposed_Area = 0.0
		for comp in xrange(ncomp):
			Area += Components[comp].Area
			if (Components[comp].Type != "internal"):
				exposed_Area += Components[comp].Area
			#end
		#end
		if (Symmetry == True):
			Area = 2.0*Area
			exposed_Area = 2.0*exposed_Area
		#end
		
		# Aspect Ratio
		AR = Span**2/Area
		
		# Wetted Area
		wetted_Area = 0.0
		for comp in xrange(ncomp):
			if (Components[comp].Type != "internal"):
				wetted_Area += Components[comp].wetted_Area
			#end
		#end
		if (Symmetry == True):
			wetted_Area = 2.0*wetted_Area		
		#end
		
		## Volume
		#Volume = 0.0
		#for comp in xrange(ncomp):
		#	Volume += Components[comp].Volume
		##end
		#if (Symmetry == True):
		#	Volume = 2.0*Volume		
		##end
		
		
		# Equivalent Lifting Surface
		# Calculate entire wing properties from n segments data
		# encapsulate everything using an instance of a lifting segment
		# Equivalent_Surface eqv also has MAC for wing
		#self.equivalent_Surface
		
		
		# Output
		self.Span = Span
		self.Area = Area
		self.AR = AR
		self.exposed_Area = exposed_Area
		self.wetted_Area = wetted_Area
		#self.Volume = Volume
		#self.Equivalent_Surface = Equivalent_Surface
		
		
	def _Surface(self):
		
		'''
		Calculate Parametric Surface
		
		Documentation last updated:  April. 7, 2008 - Ruben E. Perez
		'''
		
		## Inputs
		#Symmetry = self.Symmetry
		#Components = self._components
		#ncomp = len(Components)
		
		## 
		#if ncomp == 0:
		#	print 'Skipped Surface Calculation - Lifting Surface has no Segments'
		#	return
		##end
		
		
		## Parametric Surface Coordinates
		#
		## RHS (Outer to Inner)
		#sseg = 0
		#for seg in xrange(nseg-1,-1,-1):
		#	if (Segment[seg].Type != "internal"):
		#		[ns,nc] = numpy.shape(Segment[seg].Surface_x)
		#		xseg = numpy.zeros([ns,nc],float)
		#		yseg = numpy.zeros([ns,nc],float)
		#		zseg = numpy.zeros([ns,nc],float)			
		#		i = 0
		#		for k in xrange(ns-1,-1,-1):
		#			j = 0
		#			for l in xrange(0,nc):
		#				xseg[i,j] = Segment[seg].Surface_x[k,l]
		#				yseg[i,j] = Segment[seg].Surface_y[k,l]
		#				zseg[i,j] = Segment[seg].Surface_z[k,l]
		#				j += 1
		#			#end
		#			i += 1
		#		#end
		#		if (sseg == 0):
		#			x = xseg
		#			y = yseg
		#			z = zseg
		#		else:
		#			x = numpy.concatenate((x,xseg), axis=0)
		#			y = numpy.concatenate((y,yseg), axis=0)
		#			z = numpy.concatenate((z,zseg), axis=0)
		#		#end
		#		sseg += 1
		#	#end
		##end
		#
		## LHS (Inner to Outer)
		#if (self.Symmetry == True):
		#	for seg in xrange(nseg):
		#		if (Segment[seg].Type != "internal"):
		#			xseg = Segment[seg].Surface_x
		#			yseg = -Segment[seg].Surface_y
		#			zseg = Segment[seg].Surface_z
		#			x = numpy.concatenate((x,xseg), axis=0)
		#			y = numpy.concatenate((y,yseg), axis=0)
		#			z = numpy.concatenate((z,zseg), axis=0)
		#		#end
		#		sseg += 1
		#	#end
		##end
		
		
		## Output
		#self.Surface_x = x
		#self.Surface_y = y
		#self.Surface_z = z
		
		
	def _on_plotPlanform(self,fig,*args,**kwargs):
		
		'''
		Plot Planform (Component Dependent Part)
		
		Keyword Arguments:
		------------------
		fig: OBJECT -> Parent Plot Figure
		
		Documentation last updated:  April. 7, 2008 - Ruben E. Perez
		'''
		
		# Plot Components
		ncomp = len(self._components)
		for comp in xrange(ncomp):
			if (self._components[comp].Type != "internal"):
				fig = self._components[comp].plotPlanform(fig,*args,**kwargs)
			#end
		#end
		
		return fig
		
		
		## Planform Coordinates
		#nseg = len(self.Segment)
		#if (self.Symmetry == True):
		#	ns = 8*(nseg+1)
		#else:
		#	ns = 4*(nseg+1)
		##end
		#x = zeros(ns,float)
		#y = zeros(ns,float)
		#z = zeros(ns,float)
		#
		## RHS (Outer to Inner)
		#n = 0
		#for seg in xrange(nseg):
		#	x[n] = Segment[seg].xrLE
		#	y[n] = Segment[seg].yrLE
		#	z[n] = Segment[seg].zrLE
		#	n += 1
		#	x[n] = Segment[seg].xtLE
		#	y[n] = Segment[seg].ytLE
		#	z[n] = Segment[seg].ztLE
		#	n += 1
		##end
		#for seg in xrange(nseg-1,-1,-1):
		#	x[n] = Segment[seg].xtTE
		#	y[n] = Segment[seg].ytTE
		#	z[n] = Segment[seg].ztTE
		#	n += 1
		#	x[n] = Segment[seg].xrTE
		#	y[n] = Segment[seg].yrTE
		#	z[n] = Segment[seg].zrTE
		#	n += 1
		##end
		#
		## LHS (Inner to Outer)
		#if (self.Symmetry == True):
		#	for seg in xrange(nseg):
		#		x[n] =  Segment[seg].xrTE
		#		y[n] = -Segment[seg].yrTE
		#		z[n] =  Segment[seg].zrTE
		#		n += 1
		#		x[n] =  Segment[seg].xtTE
		#		y[n] = -Segment[seg].ytTE
		#		z[n] =  Segment[seg].ztTE
		#		n += 1
		#	#end
		#	for seg in xrange(nseg-1,-1,-1):
		#		x[n] =  Segment[seg].xtLE
		#		y[n] = -Segment[seg].ytLE
		#		z[n] =  Segment[seg].ztLE
		#		n += 1
		#		x[n] =  Segment[seg].xrLE
		#		y[n] = -Segment[seg].yrLE
		#		z[n] =  Segment[seg].zrLE
		#		n += 1
		#	#end
		##end
		#
		## Plot Planform
		#fig.append(pySciPlot_curve.PlotCurve(x,y,z,wireframe_color=(0,0,1)))
		
		
	def _on_plotSurface(self,fig,*args,**kwargs):
		
		'''
		Plot Surface (Component Dependent Part)
		
		Keyword Arguments:
		------------------
		fig: OBJECT -> Parent Plot Figure
		
		Documentation last updated:  April. 7, 2008 - Ruben E. Perez
		'''
		
		# Plot Components
		ncomp = len(self._components)
		for comp in xrange(ncomp):
			if (self._components[comp].Type != "internal"):
				fig = self._components[comp].plotSurface(fig,*args,**kwargs)
			#end
		#end
		
		return fig
		
		
	def writeOut(self,outfile,format=None):
		
		'''
		Write Formatted Output to File
		
		Keyword Arguments:
		------------------
		outfile: OBJECT -> Output File I/O Object.
		format: STRING -> Output Format [None,TeX].  Default = None
		
		Documentation last updated:  April. 7, 2008 - Ruben E. Perez
		'''
		
		# Inputs
		Components = self._components
		ncomp = len(Components)
		
		
		# Write Title Header
		outfile.write('%s - Lifting Surface\n' %(self.name))
		
		# Write Geometric Properties
		#outfile.write('\n')
		
		
		# Write Lifting Surface Components
		for comp in xrange(ncomp):
			Components[comp].writeOut(outfile,format)
		#end
		
		
	def writeDAT(self,datfile):
		
		'''
		Write to DAT File
		
		Keyword Arguments:
		------------------
		datfile: OBJECT   ->  DAT File I/O Object.
		
		Documentation last updated:  January. 29, 2009 - Ruben E. Perez
		'''
		
		# Inputs
		Components = self._components
		ncomp = len(Components)
		
		
		# Write Surface Segments
		for comp in xrange(ncomp):
			Components[comp].writeDAT(datfile)
		#end
		
		
	def writeGEO(self,geofile,nSWs=None,nCWs=None):
		
		'''
		Write to GEO File
		
		Keyword Arguments:
		------------------
		geofile: OBJECT   ->  GEO File I/O Object.
		nSW: ARRAY:   ->  Segments Spanwise Discretization.  Default = None (use Surface Discretization)
		nCW: ARRAY:   ->  Segments Chordwise Discretization.  Default = None (use Surface Discretization)
		
		Documentation last updated:  October. 20, 2007 - Ruben E. Perez
		'''
		
		# Inputs
		Components = self._components
		ncomp = len(Components)
		
		
		# Discretization Parameters
		if (nSWs == None):
			sum_nSWs = 0
			for comp in xrange(ncomp):
				if Components[comp].Type=='internal':
					sum_nSWs += Components[comp].surface_SW_segments
					sum_nSWs +=1
				else:
					sum_nSWs += Components[comp].surface_SW_segments
				#endif
			#end
		else:
			#if (len(nSWs) != nseg):
			if (len(nSWs) != ncomp):
				raise IOError('Spanwise Discretization Input Does Not Match Number of Segments')
			else:
				sum_nSWs = sum(nSWs)
			#end
		#end
		if (nCWs != None):
			if (len(nCWs) != ncomp):
			#if (len(nCWs) != nseg):
				raise IOError('Chordwise Discretization Input Does Not Match Number of Segments')
			#end
		#end
		
		
		# Undefined Parameters
		IW = 0.0
		#ALPHAW = 0.0
		
		# Write Surface Header
		geofile.write('NSEC  IW\n')
		#geofile.write('       FNC     ALPHAW\n')
		#geofile.write('%6.3f %6.3f\n' %(sum_nSWs+ncomp,IW))
		#self.nsec = sum_nSWs+ncomp
		geofile.write('%6.3f %6.3f\n' %(sum_nSWs,IW))
		self.nsec = sum_nSWs
		#geofile.write('%6.3f %6.3f\n' %(sum_nSWs+nseg,IW))
		#geofile.write('%9.3f %9.3f\n' %(sum_nSWs,ALPHAW))
		
		# Write Surface Segments
		#nseg = len(Segment)
		for comp in xrange(ncomp):
			if (nSWs == None and nCWs == None):
				Components[comp].writeGEO(geofile)
			elif (nSWs != None and nCWs == None):
				Components[comp].writeGEO(geofile,nSWs[comp],None)
			elif (nSWs == None and nCWs != None):
				Components[comp].writeGEO(geofile,None,nCWs[comp])
			elif (nSWs != None and nCWs != None):
				Components[comp].writeGEO(geofile,nSWs[comp],nCWs[comp])
			#end
		#end
	


#==============================================================================
# Class Test
#==============================================================================
if __name__ == '__main__':
	
	# Testing
	print 'Testing ...\n'
	sys.path.append(os.path.abspath('../'))
	sys.path.append(os.path.abspath('../Visualization/pySciPlot'))
	
	# -- Test Lifting Segment --
	
	
#	# Test Lifting Segment 1
#	seg = LiftingSegment()
#	seg.ListAttributes()
#	seg.plotPlanform()
#	seg.plotSurface()
#	#seg.writeGEO(open("mygeo.geo",'w'))
	
	
#	# Test Lifting Segment 2
#	input = {"Name":'Wingtip Segment',"Type":'wingtip',"Area":600.0,"Span":20.0,"Taper":0.5,"SweepLE":25.0,"Dihedral":0.0,"root_Incidence":0.0,"root_Thickness":0.12,"root_Airfoil_type":"naca4","root_Airfoil_ID":"00xx","tip_Incidence":0.0,"tip_Thickness":0.11,"tip_Airfoil_type":"naca4","tip_Airfoil_ID":"64xx"}
#	seg = LiftingSegment(input)
#	seg.ListAttributes()
#	seg.plotPlanform()
#	seg.plotSurface()
	
	
#	# Test Lifting Segment 3
#	seg = LiftingSegment(Name='Wingtip Segment',Type='wingtip',Area=600.0,Span=20.0,Taper=0.5,SweepLE=25.0,Dihedral=0.0,root_Incidence=0.0,root_Thickness=0.12,root_Airfoil_type="naca4",root_Airfoil_ID="00xx",tip_Incidence=0.0,tip_Thickness=0.11,tip_Airfoil_type="naca4",tip_Airfoil_ID="64xx")
#	seg.ListAttributes()
#	seg.plotPlanform()
#	seg.plotSurface()
	
	
#	# Test Lifting Segment 4
#	input = {"Name":'Wingtip Segment',"Type":'wingtip',"Area":600.0,"Span":20.0}
#	seg = LiftingSegment(input,Taper=0.5,SweepLE=25.0,Dihedral=0.0,root_Incidence=0.0,root_Thickness=0.12,root_Airfoil_type="naca4",root_Airfoil_ID="00xx",tip_Incidence=0.0,tip_Thickness=0.11,tip_Airfoil_type="naca4",tip_Airfoil_ID="64xx")
#	seg.ListAttributes()
#	seg.plotPlanform()
#	seg.plotSurface()
	
	
#	# Test Lifting Segment 5
#	input1 = {"Name":'Wingtip Segment',"Type":'wingtip'}
#	input2 = {"Area":600.0,"Span":20.0}
#	seg = LiftingSegment(input1,input2,Taper=0.5,SweepLE=25.0,Dihedral=0.0,root_Incidence=0.0,root_Thickness=0.12,root_Airfoil_type="naca4",root_Airfoil_ID="00xx",tip_Incidence=0.0,tip_Thickness=0.11,tip_Airfoil_type="naca4",tip_Airfoil_ID="64xx")
#	seg.ListAttributes()
#	#seg.plotPlanform()
#	seg.plotSurface()
	
	
	
	
	# -- Test Lifting Surface --
	
	
#	# Test Lifting Surface 1
#	srf = LiftingSurface()
#	srf.ListAttributes()
#	srf.plotPlanform()
#	srf.plotSurface()
#	#srf.writeGEO(open("mygeo.geo",'w'))
	
	
#	# Test Lifting Surface 2
#	input = {
#		"Name":None,
#		"Symmetry":True,
#		"xrLE":0.0,"yrLE":0.0,"zrLE":0.0,
#		"xRot":0.0,"yRot":0.0,"zRot":0.0,
#		'_components':{
#			0:LiftingSegment({"Name":"Wing Segment","Type":"external","Area":600,"Span":40,"Taper":0.6 ,"SweepLE":25.0,"Dihedral":0.0 ,"xc_offset":0.0,"root_Incidence":10.0,"root_Thickness":0.14,"root_Airfoil_type":"naca4","root_Airfoil_ID":"00xx","tip_Incidence":10.0,"tip_Thickness":0.14,"tip_Airfoil_type":"naca4","tip_Airfoil_ID":"00xx"}),
#		},
#	}
#	srf = LiftingSurface(input)
#	srf.ListAttributes()
#	srf.plotPlanform()
#	srf.plotSurface()
#	#srf.writeGEO(open("mygeo.geo",'w'))
	
	
#	# Test Lifting Surface 3
#	srf = LiftingSurface(Name="Lifting Surface",Symmetry=True,_components={0:LiftingSegment({"Name":"Wing Segment","Type":"external","Area":600,"Span":40,"Taper":0.6 ,"SweepLE":25.0,"Dihedral":0.0 ,"xc_offset":0.0,"root_Incidence":10.0,"root_Thickness":0.14,"root_Airfoil_type":"naca4","root_Airfoil_ID":"00xx","tip_Incidence":10.0,"tip_Thickness":0.14,"tip_Airfoil_type":"naca4","tip_Airfoil_ID":"00xx"})})
#	srf.ListAttributes()
#	srf.plotPlanform()
#	srf.plotSurface()
#	#srf.writeGEO(open("mygeo.geo",'w'))
	
	
#	# Test Lifting Surface 4
#	input = {"Name":"Lifting Surface","Symmetry":True}
#	seg = LiftingSegment({"Name":"Wing Segment","Type":"external","Area":600,"Span":40,"Taper":0.6 ,"SweepLE":25.0,"Dihedral":0.0 ,"xc_offset":0.0,"root_Incidence":10.0,"root_Thickness":0.14,"root_Airfoil_type":"naca4","root_Airfoil_ID":"00xx","tip_Incidence":10.0,"tip_Thickness":0.14,"tip_Airfoil_type":"naca4","tip_Airfoil_ID":"00xx"})
#	srf = LiftingSurface(input,seg)
#	srf.ListAttributes()
#	srf.plotPlanform()
#	srf.plotSurface()
#	#srf.writeGEO(open("mygeo.geo",'w'))
	
	
#	# Test Lifting Surface 5
#	input = {"Name":"Lifting Surface","Symmetry":True}
#	seg0 = LiftingSegment({"Name":"Wing Segment","Type":"external","Area":600,"Span":40,"Taper":0.6 ,"SweepLE":25.0,"Dihedral":0.0 ,"xc_offset":0.0,"root_Incidence":10.0,"root_Thickness":0.14,"root_Airfoil_type":"naca4","root_Airfoil_ID":"00xx","tip_Incidence":10.0,"tip_Thickness":0.14,"tip_Airfoil_type":"naca4","tip_Airfoil_ID":"00xx"})
#	seg1 = LiftingSegment({"Name":"Outer Segment"   ,"Type":"external","Area":800.0,"Span":50.0,"Taper":0.40,"SweepLE":30.0,"Dihedral":25.0,"xc_offset":0.0,"root_Incidence": 10.0,"root_Thickness":0.12,"root_Airfoil_type":"naca5","root_Airfoil_ID":"230xx","tip_Incidence":-20.0,"tip_Thickness":0.11,"tip_Airfoil_type":"naca5","tip_Airfoil_ID":"230xx"})	
#	srf = LiftingSurface(input,seg0,seg1,xrLE=0.0,yrLE=0.0,zrLE=0.0,xRot=0.0,yRot=0.0,zRot=0.0)
#	srf.ListAttributes()
#	srf.plotPlanform()
#	srf.plotSurface()
#	#srf.writeGEO(open("mygeo.geo",'w'))
	
	
#	# Test Lifting Surface 6
#	input = {"Name":"Lifting Surface","Symmetry":True}
#	srf = LiftingSurface(input,xrLE=0.0,yrLE=0.0,zrLE=0.0,xRot=0.0,yRot=0.0,zRot=0.0)
#	srf[0] = LiftingSegment({"Name":"Wing Segment","Type":"external","Area":600,"Span":40,"Taper":0.6 ,"SweepLE":25.0,"Dihedral":0.0 ,"xc_offset":0.0,"root_Incidence":10.0,"root_Thickness":0.14,"root_Airfoil_type":"naca4","root_Airfoil_ID":"00xx","tip_Incidence":10.0,"tip_Thickness":0.14,"tip_Airfoil_type":"naca4","tip_Airfoil_ID":"00xx"})
#	srf.ListAttributes()
#	srf.plotPlanform()
#	srf.plotSurface()
#	#srf.writeGEO(open("mygeo.geo",'w'))
	
	
#	# Test Lifting Surface 7
#	srf = LiftingSurface()
#	srf[0] = {"Name":"Wing Segment","Type":"external","Area":600,"Span":40,"Taper":0.6 ,"SweepLE":25.0,"Dihedral":0.0 ,"xc_offset":0.0,"root_Incidence":10.0,"root_Thickness":0.14,"root_Airfoil_type":"naca4","root_Airfoil_ID":"00xx","tip_Incidence":10.0,"tip_Thickness":0.14,"tip_Airfoil_type":"naca4","tip_Airfoil_ID":"00xx"}
#	srf.ListAttributes()
#	srf.plotPlanform()
#	srf.plotSurface()
#	#srf.writeGEO(open("mygeo.geo",'w'))
	
	
#	# Test Lifting Surface 8
#	input = {
#		"Name":None,
#		"Symmetry":True,
#		"xrLE":0.0,"yrLE":0.0,"zrLE":0.0,
#		"xRot":0.0,"yRot":0.0,"zRot":0.0,
#		'_components':{
#			0:{"Name":"Internal Segment","Type":"internal","Area":575.0,"Span":12.3529,"Taper":0.7,"SweepLE":25.0,"Dihedral":0.0,"xc_offset":0.0,"root_Incidence": 1.0,"root_Thickness":0.14,"root_Airfoil_type":"naca5","root_Airfoil_ID":"230xx","tip_Incidence": 1.0,"tip_Thickness":0.14,"tip_Airfoil_type":"naca5","tip_Airfoil_ID":"230xx"},
#			1:{"Name":"Inner Segment","Type":"external","Area":600.0,"Span":19.6875,"Taper":0.6,"SweepLE":25.0,"Dihedral":3.0,"xc_offset":0.0,"root_Incidence": 1.0,"root_Thickness":0.14,"root_Airfoil_type":"naca5","root_Airfoil_ID":"230xx","tip_Incidence":-2.0,"tip_Thickness":0.12,"tip_Airfoil_type":"naca5","tip_Airfoil_ID":"230xx"},
#			2:{"Name":"Outer Segment","Type":"external","Area":800.0,"Span":50.0,"Taper":0.4,"SweepLE":30.0,"Dihedral":6.0,"xc_offset":0.0,"root_Incidence":-2.0,"root_Thickness":0.12,"root_Airfoil_type":"naca5","root_Airfoil_ID":"230xx","tip_Incidence":-3.0,"tip_Thickness":0.11,"tip_Airfoil_type":"naca5","tip_Airfoil_ID":"230xx"},
#			3:{"Name":"Wingtip Segment","Type":"wingtip" ,"Area":12.0 ,"Span":1.5,"Taper":0.75,"SweepLE":60.0,"Dihedral":6.0,"xc_offset":0.0,"root_Incidence":-3.0,"root_Thickness":0.11,"root_Airfoil_type":"naca5","root_Airfoil_ID":"230xx","tip_Incidence": 0.0,"tip_Thickness":0.10,"tip_Airfoil_type":"naca5","tip_Airfoil_ID":"230xx"},
#			4:{"Name":"Winglet Segment","Type":"winglet" ,"Area":50.0 ,"Span":9.0,"Taper":0.5 ,"SweepLE":45.0,"Dihedral":90.0,"xc_offset":0.0,"root_Incidence": 0.0,"root_Thickness":0.10,"root_Airfoil_type":"naca5","root_Airfoil_ID":"230xx","tip_Incidence": 0.0,"tip_Thickness":0.00,"tip_Airfoil_type":"naca4","tip_Airfoil_ID":"22xx"},
#		},
#	}
#	srf = LiftingSurface(input)
#	srf.ListAttributes()
#	srf.plotPlanform()
#	srf.plotSurface()
#	srf.writeGEO(open("mygeo.geo",'w'))
#	pdb.set_trace()
	
	
	
#	# Test Lifting Surface 9
#	input = {
#		"Name":None,
#		"Symmetry":True,
#		"xrLE":0.0,"yrLE":0.0,"zrLE":0.0,
#		"xRot":0.0,"yRot":0.0,"zRot":0.0,
#		'_components':{
#			0:{"Name":"Internal Segment","Type":"internal","Area":575.0,"Span":12.3529,"Taper":0.7,"SweepLE":25.0,"Dihedral":0.0,"xc_offset":0.0,"root_Incidence": 1.0,"root_Thickness":0.14,"root_Airfoil_type":"naca5","root_Airfoil_ID":"230xx","tip_Incidence": 1.0,"tip_Thickness":0.14,"tip_Airfoil_type":"naca5","tip_Airfoil_ID":"230xx"},
#			1:{"Name":"Inner Segment","Type":"external","Area":600.0,"Span":19.6875,"Taper":0.6,"SweepLE":25.0,"Dihedral":3.0,"xc_offset":0.0,"root_Incidence": 1.0,"root_Thickness":0.14,"root_Airfoil_type":"naca5","root_Airfoil_ID":"230xx","tip_Incidence":-2.0,"tip_Thickness":0.12,"tip_Airfoil_type":"naca5","tip_Airfoil_ID":"230xx",'_components':{0:{"Spoiler":{}},1:{"Slat":{}},2:{"Flap":{}}}},
#			2:{"Name":"Outer Segment","Type":"external","Area":800.0,"Span":50.0,"Taper":0.4,"SweepLE":30.0,"Dihedral":6.0,"xc_offset":0.0,"root_Incidence":-2.0,"root_Thickness":0.12,"root_Airfoil_type":"naca5","root_Airfoil_ID":"230xx","tip_Incidence":-3.0,"tip_Thickness":0.11,"tip_Airfoil_type":"naca5","tip_Airfoil_ID":"230xx"},
#			3:{"Name":"Wingtip Segment","Type":"wingtip" ,"Area":12.0 ,"Span":1.5,"Taper":0.75,"SweepLE":60.0,"Dihedral":6.0,"xc_offset":0.0,"root_Incidence":-3.0,"root_Thickness":0.11,"root_Airfoil_type":"naca5","root_Airfoil_ID":"230xx","tip_Incidence": 0.0,"tip_Thickness":0.10,"tip_Airfoil_type":"naca5","tip_Airfoil_ID":"230xx"},
#			4:{"Name":"Winglet Segment","Type":"winglet" ,"Area":50.0 ,"Span":9.0,"Taper":0.5 ,"SweepLE":45.0,"Dihedral":90.0,"xc_offset":0.0,"root_Incidence": 0.0,"root_Thickness":0.10,"root_Airfoil_type":"naca5","root_Airfoil_ID":"230xx","tip_Incidence": 0.0,"tip_Thickness":0.00,"tip_Airfoil_type":"naca4","tip_Airfoil_ID":"22xx"},
#		},
#	}
#	srf = LiftingSurface(input)
#	srf.ListAttributes()
#	srf.plotPlanform()
#	srf.plotSurface()
#	#srf.writeGEO(open("mygeo.geo",'w'))
	
