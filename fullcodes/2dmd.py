#///////////////////////////////////////////////////////////////////////////////////////////////////////////
# This material is released under a Creative Commons License
# (Attribution-Noncommercial-Share Alike License)
# (http://creativecommons.org/licenses/by-nc-sa/3.0/) with the following attribution:
# "Original version of this material was developed by the Center for Computational Biology at UC Merced." 
#///////////////////////////////////////////////////////////////////////////////////////////////////////////
#
# 2D Molecular Dynamics Simulation program
# By Dr. Masa Watanabe, Center for Computational Biology, UC Merced
# Version 0.65
#
# With minor bug fixes and updates by Prof. Mike Colvin UC Merced
#
# This program implements a simple 2D Molecular dynamics simulation of
# atoms with Lennard-Jones and other potentials.
#
# Bugfix
# 05/28/2009: Initialization of kinetic energy for the Berendsen thermostat was fixed
# 05/28/2009: fixed to reflect the correct mass on the slider interface
#

from tkinter import *
from tkinter.filedialog import askopenfilename, asksaveasfilename
import random
import math, sys, string

# The graphical interface
# TKinter 
class GUI:
    global Width, Height, WW, W2, space
    Width, Height = 500, 500    
    WW = Width-25.0
    W2 = 25.0/2.0
    space = 0.25

    def __init__(self):
        self.tk = tk = Tk()
        tk.title('2D Molecular Dynamics Simulator')

        self.selection = 'MD Simulation'

        self.frame3 = frame3 = Frame(tk,width=300,height=20,name='frame3')
        frame3.pack(expand=1,fill=BOTH, side=LEFT)
        self.frame1 = frame1 = Frame(tk,width=Width,height=Height,name='frame1',bg='#FFE7BA')
        frame1.pack(expand=1,fill=BOTH,side=LEFT)
        self.frame2 = frame2 = Frame(tk,width=300,height=Height,name='frame2')
        frame2.pack(expand=1,fill=BOTH, side=LEFT)
        
        self.canvas = canvas = Canvas(frame1,width=Width,height=Height,bg='#CDBA96')
        canvas.pack(side=TOP)

        self.canvas3 = canvas3 = Canvas(frame3,width=275,height=200,bg='#FFE7BA')
        canvas3.pack(side=BOTTOM)

        self.canvas5 = canvas5 = Canvas(frame3,width=275,height=200,bg='#B9D3EE')
        canvas5.pack(side=BOTTOM)
        
        self.canvas2 = canvas2 = Canvas(frame3,width=275,height=150,bg='yellow')
        canvas2.pack(side=TOP)

        canvas2.create_text(125,20,text="Properties\n",font='Helvetica 10 bold')

        self.canvas4 = canvas4 = Canvas(frame1,width=Width,height=50,name='frame4')
        text1 = 'University of California, Merced\nCenter for Computational Biology'
        canvas4.create_text(260,20,text=text1,justify=CENTER,font='Helvetica 10',fill="#A0522D",tag='text1')
        
        canvas4.pack(side=BOTTOM)
        self.print_num = 0
        self.gui_initial = 0
        self.hitting_counter = 0

        self.histogram = histogram(20,0.,2.0)
        self.histogram_veloX = histogram(20,-2.0,2.0)
        self.histogram_veloY = histogram(20,-2.0,2.0)
         
    def Quit(self,md):
        md._running = 0
        self.tk.quit()
        
    def setup(self,md):
        n = md.number_atoms
        canvas = self.canvas
        c3 = self.canvas3
        canvas.delete('lines')
        c3.delete('lines')
        c5 = self.canvas5
        c5.delete('lines2')

        c3.delete('hist1')
        c5.delete('hist')

        BOX1, BOX2, BOX3 = 50.0, 37.5, 12.5
            
        c3.create_line(25,175,250,175,fill='#333366')
        c3.create_line(25,175,25,35,fill='#333366')
        
        c5.create_line(25,175,250,175,fill='#333366')
        c5.create_line(25,175,25,35,fill='#333366')
        
        centerX = 25
        centerY = 175

        bin_width = float(225/self.histogram.nbins)

        c5.create_text(125,15,text='Distribution',tag='lines5',font='Helvetica 10 bold')
         
        for i in range(0,self.histogram.nbins+1,5):
            c3.create_line(centerX+i*bin_width,centerY,centerX+i*bin_width,centerY+10,tag='lines')
            tt0 = i * float(self.histogram.binend/self.histogram.nbins)
            tt1 = str(round(tt0,2))
            c3.create_text(centerX+i*bin_width,centerY+13,text=tt1,tag='lines')
            c5.create_line(centerX+i*bin_width,centerY,centerX+i*bin_width,centerY+10,tag='lines2')
            c5.create_text(centerX+i*bin_width,centerY+13,text=tt1,tag='lines2')
            
        delta_height = float(225.0/md.number_atoms)
        for i in range(5):
            c3.create_line(centerX-5,centerY-i*70.0,centerX,centerY-i*70.0,tag='lines')
            tt1 = round(i*0.25,1)
            if(tt1%0.5==0):
                tt0 = str(tt1)
                c3.create_text(centerX-13,centerY-i*70.0,text=tt0,tag='lines')  # 35.0 as 1 unit
            tt0 = str(int(round(i*md.number_atoms/5,0)))
            c5.create_line(centerX-5,centerY-i*45.0,centerX,centerY-i*45.0,tag='lines2')
            c5.create_text(centerX-10,centerY-i*45.0,text=tt0,tag='lines2')

        W1 = WW/md.BoxL
        T1 = 5.0*10.0/md.BoxL
        box = md.BoxL
        if md.PBC:
            boxinv = 1.0 / box
        else:
            boxinv = 0.0

        for i in range(n):
            local_x = md.atoms[i].atom.x
            local_x -= round(local_x * boxinv) * box
            local_y = md.atoms[i].atom.y
            local_y -= round(local_y * boxinv) * box
            x = W2 + W1 * abs(local_x+md.BoxL_h)
            x1, x2 = x-T1,x+T1
            y = W2 + W1 * abs(local_y-md.BoxL_h)
            y1, y2 = y-T1,y+T1
            w = 1
            if (md.atoms[i].trace):
                w = 2
            canvas.create_oval(x1,y1,x2,y2,fill=md.atoms[i].color,outline='#333333',width=w,tag='points')
        
        self.tk.update()
        
    def md_delete(self,md,n):
        iatom = md.number_atoms
        dn = md.number_atoms - n

        i = iatom - 1
        #for j in range(dn):
        #    ix = int((md.atoms[i].atom.x+md.BoxL_h)/space)
        #    iy = int((md.atoms[i].atom.y+md.BoxL_h)/space)
        #    md.gridX[ix] = md.gridY[iy] = 0
        #    i -= 1

        for j in range(dn):
            del md.atoms[iatom-1]
            iatom -=1
            md.number_atoms -= 1

        md.histogram1.data  = [0 for i in range(md.histogram1.nbins)]
        md.hist1_veloX.data = [0 for i in range(md.hist1_veloX.nbins)]
        md.hist1_veloY.data = [0 for i in range(md.hist1_veloY.nbins)]
         
        self.tk.update()
        
    def md_move(self,md):
        c = self.canvas
        c2 = self.canvas2

        histogram = self.histogram   ## Instantaneous Distribution
        hist_veloX = self.histogram_veloX
        hist_veloY = self.histogram_veloY
        
        c.delete('points')
        c.delete('text')
        c.delete('lines')
        c2.delete('text')
        self.canvas3.delete('hist1')
        self.canvas5.delete('hist')
        self.canvas3.delete('text')
        self.canvas5.delete('lines3')
        self.canvas3.delete('lines3')
        self.canvas5.delete('lines')
        self.canvas3.delete('lines')
        self.canvas5.delete('lines2')
        self.canvas3.delete('lines2')
        
        i = 0
        n = md.number_atoms

        xy1 = -25.0,-25.0,50.0,50.0
        xy2 = -25.0,Height+25.0,50.0,Height-50.0
        xy3 = Width+25.0,-25.0,Width-50.0,50.0
        xy4 = Width+25.0,Height+25.0,Width-50.0,Height-50.0
        xy5 = Width-12.5,12.5+Height/2.0,Width-37.5,-12.5+Height/2.0
        xy6 = 12.5,12.5+Height/2.0,37.5,-12.5+Height/2.0

        W1 = WW/md.BoxL
        T1 = 5.0*10.0/md.BoxL
        box = md.BoxL
        if md.PBC:
            boxinv = 1.0 / box
        else:
            boxinv = 0.0

        histogram.data  = [0 for i in range(histogram.nbins)]
        hist_veloX.data = [0 for i in range(hist_veloX.nbins)]
        hist_veloY.data = [0 for i in range(hist_veloY.nbins)]
        
        for i in range(n):
            local_x = md.atoms[i].atom.x
            local_x -= round(local_x * boxinv) * box
            local_y = md.atoms[i].atom.y
            local_y -= round(local_y * boxinv) * box
            x = W2 + W1 * abs(local_x+md.BoxL_h)
            y = W2 + W1 * abs(local_y-md.BoxL_h)
            x1, x2 = x-T1, x+T1
            y1, y2 = y-T1, y+T1
            w=1
            if (md.atoms[i].trace):
                if (md.trace):
                    last_x = md.atoms[i].atom.last_x
                    last_x -= round(last_x * boxinv) * box
                    last_y = md.atoms[i].atom.last_y
                    last_y -= round(last_y * boxinv) * box
                    if (abs(local_x - last_x) < md.BoxL_h and abs(local_y - last_y) < md.BoxL_h):
                        l_x = W2 + W1 * abs(last_x+md.BoxL_h)
                        l_y = W2 + W1 * abs(last_y-md.BoxL_h)
                        c.create_line(x,y,l_x,l_y,fill=md.atoms[i].color,tag='traces')
                w=2
            md.atoms[i].atom.last_x = md.atoms[i].atom.x
            md.atoms[i].atom.last_y = md.atoms[i].atom.y
            c.create_oval(x1,y1,x2,y2,fill=md.atoms[i].color,outline='#333333',width=w,tag='points')
            vx = md.atoms[i].atom.vx
            vy = md.atoms[i].atom.vy
            spd = math.sqrt(vx*vx+vy*vy)
            histogram.add(spd,0.0)
            hist_veloX.add(vx,2.0)
            hist_veloY.add(vy,2.0)
            
        md.histogram1.statistic()
        md.hist1_veloX.statistic()
            
        if(md.md_anlysis == 0 ) :

            centerX = 25
            centerY = 175

            bin_width = float(225/self.histogram.nbins)
        
            for i in range(0,self.histogram.nbins+1,5):
                self.canvas3.create_line(centerX+i*bin_width,centerY,centerX+i*bin_width,centerY+10,tag='lines')
                tt0 = i * float(self.histogram.binend/self.histogram.nbins)
                tt1 = str(round(tt0,2))
                self.canvas3.create_text(centerX+i*bin_width,centerY+13,text=tt1,tag='lines')
                self.canvas5.create_line(centerX+i*bin_width,centerY,centerX+i*bin_width,centerY+10,tag='lines2')
                self.canvas5.create_text(centerX+i*bin_width,centerY+13,text=tt1,tag='lines2')
                
            for i in range(5):
                self.canvas3.create_line(centerX-5,centerY-i*70.0,centerX,centerY-i*70.0,tag='lines')
                tt1 = round(i*0.25,1)
                if(tt1%0.5==0):
                    tt0 = str(tt1)
                    self.canvas3.create_text(centerX-13,centerY-i*70.0,text=tt0,tag='lines')  # 35.0 as 1 unit
                tt0 = str(int(round(i*(md.number_atoms+1)/5,0)))
                self.canvas5.create_line(centerX-5,centerY-i*35.0,centerX,centerY-i*35.0,tag='lines2')
                self.canvas5.create_text(centerX-10,centerY-i*35.0,text=tt0,tag='lines2')
            
            delta_height = float(150.0/md.number_atoms)

            self.canvas5.create_text(125,30,text='Speed (Instantaneous)',tag='lines3',font='Helvetica 9')
            self.canvas3.create_text(125,30,text='Speed (Accumulated)  ',tag='lines3',font='Helvetica 9')                    

            for i in range(histogram.nbins):
                scaledY = centerY - histogram.data[i]*delta_height
                self.canvas5.create_rectangle(centerX+i*bin_width,centerY,centerX+(i+1)*bin_width,scaledY,fill='yellow', tag='hist')

            sum_dist = 0
            
            for i in range(md.histogram1.nbins):
                sum_dist += md.histogram1.data[i]

            V_Mean = md.histogram1.Mean
            V_STdev = md.histogram1.STdev
            
            delta_height = 2.0*float(150.0/sum_dist)
            for i in range(md.histogram1.nbins): # Normalized histogram
                scaledY = centerY - md.histogram1.data[i]*delta_height
                self.canvas3.create_rectangle(centerX+i*bin_width,centerY,centerX+(i+1)*bin_width,scaledY,fill='#EE7621', tag='hist1')

            text1 = 'Ave.   = ' + str(round(V_Mean,3)) + '\n'
            text2 = 'Stdev. = ' + str(round(V_STdev,3)) + '\n'
            text0 = text1 + text2
            self.canvas3.create_text(190,60,text=text0,justify=LEFT,font='Helvetica 9',tag='text')

# End of md_analysis == 1 (Speed)

        else:
 
            centerX = 25
            centerY = 175

            bin_width = float(225/self.histogram.nbins)
            
            for i in range(0,self.histogram.nbins+1,5):
                self.canvas3.create_line(centerX+i*bin_width,centerY,centerX+i*bin_width,centerY+10,tag='lines')
                tt0 = i * float(self.histogram.binend/self.histogram.nbins)-1.0
                tt1 = str(round(tt0,2))
                self.canvas3.create_text(centerX+i*bin_width,centerY+13,text=tt1,tag='lines')
                self.canvas5.create_line(centerX+i*bin_width,centerY,centerX+i*bin_width,centerY+10,tag='lines2')
                self.canvas5.create_text(centerX+i*bin_width,centerY+13,text=tt1,tag='lines2')
                
            for i in range(5):
                self.canvas3.create_line(centerX-5,centerY-i*70.0,centerX,centerY-i*70.0,tag='lines')
                tt1 = round(i*0.25,1)
                if(tt1%0.5==0):
                    tt0 = str(tt1)
                    self.canvas3.create_text(centerX-13,centerY-i*70.0,text=tt0,tag='lines')  # 35.0 as 1 unit
                tt0 = str(int(round(i*(md.number_atoms+1)/5,0)))
                self.canvas5.create_line(centerX-5,centerY-i*35.0,centerX,centerY-i*35.0,tag='lines2')
                self.canvas5.create_text(centerX-10,centerY-i*35.0,text=tt0,tag='lines2')
                
            self.canvas5.create_text(125,30,text='Velocity (x) (Instantaneous)',tag='lines3',font='Helvetica 9')
            self.canvas3.create_text(125,30,text='Velocity (x) (Accumulated)  ',tag='lines3',font='Helvetica 9')                    

            delta_height = float(150.0/md.number_atoms)
            
            for i in range(hist_veloX.nbins):
                scaledY = centerY - hist_veloX.data[i]*delta_height
                self.canvas5.create_rectangle(centerX+i*bin_width,centerY,centerX+(i+1)*bin_width,scaledY,fill='#DAA520', tag='hist')

            sum_dist = 0
            
            for i in range(md.hist1_veloX.nbins):
                sum_dist += md.hist1_veloX.data[i]

            V_Mean = md.hist1_veloX.Mean
            V_STdev = md.hist1_veloX.STdev
            
            delta_height = 2.0*float(140.0/sum_dist)
            for i in range(md.hist1_veloX.nbins): # Normalized histogram
                scaledY = centerY - md.hist1_veloX.data[i]*delta_height
                self.canvas3.create_rectangle(centerX+i*bin_width,centerY,centerX+(i+1)*bin_width,scaledY,fill='#A2CD5A', tag='hist1')

            text1 = 'Ave.   = %5.3lf \n'%(round(V_Mean,3))
            text2 = 'St. dev. = ' + str(round(V_STdev,3)) + '\n'
            text3 = 'Var. = %5.3lf\n'%(round(V_STdev*V_STdev,3))
            text0 = text1 + text3
            self.canvas3.create_text(190,60,text=text0,justify=LEFT,font='Helvetica 9',tag='text')
                
          
        K_Mean = float(md.K_sum)/float(md.number_of_steps)
        K_STdev = math.sqrt(float(md.number_of_steps)*md.K_sum_2-md.K_sum*md.K_sum)/float(md.number_of_steps)
        
        # Try to print average
        text1  = '\n\nTime = %5.2lf\n'%(round(md.time,2))
        text2  = 'Total Energy = %5.2lf\n'%(round(md.Total_E,1))
        text22 = 'Kinetic Energy(KE) = %5.2lf\n'%(round(md.K,2))
        #text23 = '<KE> = ' + str(round(K_Mean,2))+'   Stdev = ' + str(round(K_STdev,2))+'\n'
        text23 = '<KE> = %5.2lf\n'%(round(K_Mean,2))

        if(md.md_constantT != 0):
            text3 = 'Temperature = %5.2lf\n'%(round(md.T_desire,2))
        else:
            text3 = 'Temperature = %5.2lf\n'%(round(md.K/float(md.number_atoms),2))

        v2=md.histogram1.SRMS*md.histogram1.SRMS
        text4 = '<V*V>   = %5.2lf'%(round(v2,3))
        #text5 = '\nrms Veloc. (x) = ' + str(round(md.hist1_veloX.SRMS,3))
        #text0 = text1 + text2 + text22 + text23 + text3 + text4 + text5
        text0 = text1 + text3 + text2 + text22 + text23 + text4
        self.canvas2.create_text(120,60,text=text0,justify=LEFT,font='Helvetica 10',tag='text')
        
        self.tk.update()

    def md_cleanup(self):
        c = self.canvas
        c.delete('points')
        c.delete('text')
        c.delete('corner')
        c.delete('lines')
        c.delete('traces')
        self.canvas3.delete('phase')
        self.hitting_counter = 0
        
# Utility class for building histograms
class histogram:         
    def __init__(self, nbins, binstart, binend):
        self.maxbars=30     #Maximum number bars in printed histogram
        self.nbins=nbins
        self.binstart=binstart
        self.binend=binend
        self.binwidth=(float(binend)-float(binstart))/float(nbins)
        self.bins=[]
        self.sum_x=0.0
        self.sum_x2=0.0
        self.n=0
        self.Mean = 0.0
        self.STdev = 0.0
        self.SRMS = 0.0
        self.bins.append(binstart)
        for i in range(nbins):
            self.bins.append((i+1)*self.binwidth)
        self.data=[0 for i in range(nbins)]
        
    def add(self,data,p):
##        if (data < self.binstart or data > self.binend):
##            print "data outside of histogram range"
##            #return
            
        i = 1
        while ((data+p) > self.bins[i]):
            i+=1
            if(i > self.nbins):
                i = self.nbins
                break
            
        self.data[i-1]+=1
        self.sum_x+=data
        self.sum_x2+=data*data
        self.n+=1
        
    def reset(self):
        self.data=[0 for i in range(nbins)]
        self.sum_x=0.0
        self.sum_x2=0.0
        self.n=0
        self.Mean = 0.0
        self.STdev = 0.0
        self.SRMS = 0.0
        
    def print_data(self):
        for i in range(self.nbins):
            print('%lf2.2-%lf2.2: %-4d' % (self.bins[i],self.bins[i+1],
                                           self.data[i]))
    def statistic(self):
        self.Mean = float(self.sum_x)/float(self.n)
        self.SRMS = math.sqrt(float(self.sum_x2)/float(self.n))
        self.STdev = math.sqrt(float(self.n)*self.sum_x2-self.sum_x*self.sum_x)/float(self.n)
        
    def print_hist(self):
        maxdata=0;
        mean=float(self.sum_x)/float(self.n)
        stdev=math.sqrt(float(self.n)*self.sum_x2-self.sum_x*self.sum_x)/float(self.n)
        for i in range(self.nbins):
           if self.data[i]>maxdata:
               maxdata=self.data[i]
        print('Mean: %4.2lf' % (mean))
        print('Standard deviation: %4.2lf' % (stdev))
        for i in range(self.nbins):
            bars=int(self.maxbars*float(self.data[i])/float(maxdata))
            print('%2.2lf-%2.2lf: %-4d %s' % (self.bins[i],self.bins[i+1],
                                           self.data[i],bars*"="))
        
# Atom Class
class Atom:
    # Create our objects
    def __init__(self, x, y, mass):
        self.x = x
        self.y = y
        self.last_x = x
        self.last_y = y
        self.mass = mass
        self.mu = 0.0
        self.charge = 0.0
        self.all_stopX = 0
        self.all_stopY = 0
        
    def velocity(self,vx,vy):
        self.vx = vx
        self.vy = vy

    def forces(self,fx,fy):
        self.fx = fx
        self.fy = fy
    
    def moveBall(self,dt): # Velocity Verlet Method - Part I
        dt2 = dt / 2.0
        dtsq2 = dt * dt2
        self.x = self.x + dt * self.vx + dtsq2 * self.fx / self.mass
        self.y = self.y + dt * self.vy + dtsq2 * self.fy / self.mass
        self.vx = self.vx + dt2 * self.fx / self.mass
        self.vy = self.vy + dt2 * self.fy / self.mass

    def moveBall2(self,dt): # Velocity Verlet Method - Part II
        dt2 = dt / 2.0
        self.vx = self.vx + dt2 * self.fx / self.mass
        self.vy = self.vy + dt2 * self.fy / self.mass

class System:
    # Create systems
    def __init__(self,myid,atom_name,atom):
        self.atom = atom
        self.myid = myid
        self.name = atom_name
        self.color = '#FFCC00'
        self.trace = False
        
############### Define Python classes for simulation ########################
# This is a utility class for making random numbers
class ranclass:
    def ranint(self, n):
        return int(n*random.random())
    def ranfloat(self):
        return random.random()
    
###################
# MD
class MD:
    def __init__(self):
        self._running = 0
        self.number_of_steps = 0
        self.time = 0
        self.friction = 0
        self.static = 0
        self.hitting_counter = 0
        self.loop1 = 0
        self.md_potential = 0
        self.md_constantT = 2
        self.interaction_strength = 1.0
        self.timestep = 1.0
        self.T_desire = 1.0
        self.tau_T = 2.0
        self.trace = True

        self.md_anlysis = 1  # Analysis category (0 - Speed, 1 - Velocity)

        self.Total_E = 0.0
        self.K = 0.0
        self.K_sum = 0.0
        self.K_sum_2 = 0.0

        self.histogram1 = histogram(20,0.,2.0)
        self.hist1_veloX = histogram(20,-2.,2.0)
        self.hist1_veloY = histogram(20,-2.,2.0)
        
    def MD_reset(self):
        self.number_of_steps = 0
        self.time = 0
        self.hitting_counter = 0
        self.loop1 = 0
        self.Total_E = 0.0
        self.K = 0.0
        self.K_sum = 0.0
        self.K_sum_2 = 0.0
        self.histogram1 = histogram(20,0.,2.0)
        self.hist1_veloX = histogram(20,-2.,2.0)
        self.hist1_veloY = histogram(20,-2.,2.0)
##        self.md_anlysis = 1
        
    def MD_cleanup(self):
        self.number_of_steps = 0
        self.static = 0
        self.loop1 = 0
        
    def set_friction(self,mu):
        self.friction = 1
        for i in range(self.number_atoms):
            self.atoms[i].atom.mu = mu
            
    def setup(self,n,BoxL,PBC):

        self.number_atoms = n
        self.initial_atoms = n
        self.BoxL = BoxL
        self.BoxL_h = BoxL/2.0
        self.PBC = PBC

        self.atoms=[]
        iatom=0
        
        x, y = 0.0, 0.0
        B1 = BoxL-0.5
        
        self.ngrid = int(BoxL/space)

        jg = []
        self.grid=[]
        for i in range(1,self.ngrid):
            for j in range(1,self.ngrid):
                xy = i * space - self.BoxL_h, j*space - self.BoxL_h
                self.grid.append(xy)
                
        random.seed()
        ngrid2 = (self.ngrid-1) * (self.ngrid-1)
        jg = random.sample(range(ngrid2),ngrid2)
        
        for i in range(n):
            j = jg[i]
            x = self.grid[j][0] #* space - self.BoxL_h
            y = self.grid[j][1] #* space - self.BoxL_h
            self.atoms.append(System(iatom,"LJ",Atom(x,y,20.0))) ## Set up the mass and initial cooridnate
            iatom += 1
            #if(i >= ngrid2): exit

        vx, vy = random.random()-0.5, random.random()-0.5
        for i in range(n):
            self.atoms[i].atom.velocity(vx,vy)
            vx = random.random()-0.5
            vy = random.random()-0.5

        if(self.md_potential == 0):
            self.forces_LJ() # Calculating initial forces
        elif(self.md_potential == 1):
            self.forces_1(10.0)
        elif(self.md_potential == 2):
            self.forces_1(-10.0)

        self.histogram1.data=[0 for i in range(self.histogram1.nbins)]
            
    def add(self,n):

        iatom = n_original = self.number_atoms
        self.initial_atoms = n
        dn = n - iatom
        B1 = self.BoxL - 0.5
                
        for i in range(dn):
            post = 1
            
            while post:
                x = B1 *(random.random()-0.5)
                y = B1 *(random.random()-0.5)
                for j in range(iatom):
                    x1 = x-self.atoms[j].atom.x
                    y1 = y-self.atoms[j].atom.y
                    sq2 = math.sqrt(x1*x1+y1*y1)
                    if(sq2 < space):
                        post = 1
                        break
                    post = 0
            self.atoms.append(System(iatom,"LJ",Atom(x,y,20.0)))
            iatom +=1

        vx, vy = random.random()-0.5, random.random()-0.5
        for i in range(n_original,n):
            self.atoms[i].atom.velocity(vx,vy)
            vx = random.random()-0.5
            vy = random.random()-0.5

        self.number_atoms = n

        
        for i in range(n): # initialize all forces
            self.atoms[i].atom.forces(0.0,0.0)

        if(self.md_potential == 0):
            self.forces_LJ() # Calculating initial forces
        elif(self.md_potential == 1):
            self.forces_1(10.0)
        elif(self.md_potential == 2):
            self.forces_1(-10.0)
        else:
            self.forces_LJ()
        
    def loop(self,step,dt):
    #
    # Main Part of Molecular Dynamic Simulation Loop
    #
        
        n = self.number_atoms
        self.number_of_steps += 1
        self.loop1 = 1
        box = self.BoxL
        boxinv = 1.0 / box
        Kinetic_E = 0.0

        self.time += dt

        if not self.PBC:
            for i in range(n):
                local_x = self.atoms[i].atom.x
                local_y = self.atoms[i].atom.y
                if (local_x > self.BoxL_h) or (local_x < -self.BoxL_h):
                    self.atoms[i].atom.vx = -1.0 * self.atoms[i].atom.vx
                if (local_y > self.BoxL_h) or (local_y < -self.BoxL_h):
                    self.atoms[i].atom.vy = -1.0 * self.atoms[i].atom.vy

        for i in range(n):
            Kinetic_E += (self.atoms[i].atom.vx * self.atoms[i].atom.vx + self.atoms[i].atom.vy * self.atoms[i].atom.vy)
            self.atoms[i].atom.moveBall(dt)

        self.V = 0.0
        if(self.md_potential == 0):
            self.V = self.forces_LJ() # Potential Energy and force calculation
        elif(self.md_potential == 1):
            self.V = self.forces_1(10.0)
        elif(self.md_potential == 2):
            self.V = self.forces_1(-10.0)
        else:
            self.V = self.forces_LJ()
                    
        self.K = 0.0 # Kinetic Energy
        self.static = 0
        Kinetic_E1 = 0.0
        for i in range(n):
            self.atoms[i].atom.moveBall2(dt)
            self.static = self.static + self.atoms[i].atom.all_stopX + self.atoms[i].atom.all_stopY
            Kinetic_E1 += self.atoms[i].atom.vx * self.atoms[i].atom.vx + self.atoms[i].atom.vy * self.atoms[i].atom.vy
            self.K = self.K + 0.5 * self.atoms[i].atom.mass*((self.atoms[i].atom.vx)**2 + (self.atoms[i].atom.vy)**2)

        #self.K = 0.0
        #if((self.md_potential == 3) or (self.md_constantT == 1)):
        if(self.md_constantT == 1):
            Kinetic_E = n * self.T_desire
            if(Kinetic_E != 0.0):
                velocity_scale = math.sqrt(Kinetic_E1/Kinetic_E)
                #velocity_scale = 1.0
            else:
                velocity_scale = 1.0

            self.K = 0.0 # fixed 05/28/2009
            for i in range(n):
                self.atoms[i].atom.vx = velocity_scale * self.atoms[i].atom.vx
                self.atoms[i].atom.vy = velocity_scale * self.atoms[i].atom.vy
                self.K = self.K + 0.5 * self.atoms[i].atom.mass*((self.atoms[i].atom.vx)**2 + (self.atoms[i].atom.vy)**2)
        elif(self.md_constantT == 2):
            # Berendesen constant T routine
            Temp = self.K/n
            velocity_scale = 1 + dt * (self.T_desire/Temp - 1.0) / (2.0 * self.tau_T)
            self.K = 0.0 # fixed 05/28/2009
            for i in range(n):
                self.atoms[i].atom.vx = velocity_scale * self.atoms[i].atom.vx
                self.atoms[i].atom.vy = velocity_scale * self.atoms[i].atom.vy
                self.K = self.K + 0.5 * self.atoms[i].atom.mass*((self.atoms[i].atom.vx)**2 + (self.atoms[i].atom.vy)**2)

        for i in range(n):
            vx = self.atoms[i].atom.vx
            vy = self.atoms[i].atom.vy
            spd = math.sqrt(vx*vx+vy*vy)
            self.histogram1.add(spd,0)
            self.hist1_veloX.add(vx,2.0)
            self.hist1_veloY.add(vy,2.0)
            
        self.K_sum += self.K
        self.K_sum_2 += self.K * self.K
        self.Total_E = self.K + self.V
        
    def forces_LJ(self):
    #
    # Simple Lennard-Jones Force Routine
    #

        n = self.number_atoms

        sigma = 0.21
        epslon = 1.0 * self.interaction_strength

        sigsq = sigma**2
        eps4 = epslon * 4.0
        eps24 = epslon * 24.0

        box = self.BoxL
        if self.PBC:
            boxinv = 1.0 / box
        else:
            boxinv = 0.0

        V, W = 0.0, 0.0

        for i in range(n):
            self.atoms[i].atom.forces(0.0,0.0)
            self.atoms[i].all_stop = 0

        for i in range(n-1):
            RXI = self.atoms[i].atom.x
            RYI = self.atoms[i].atom.y
            FXI = self.atoms[i].atom.fx
            FYI = self.atoms[i].atom.fy

            for j in range(i+1,n):
                
                rxij = RXI - self.atoms[j].atom.x
                ryij = RYI - self.atoms[j].atom.y
                ## JLP - needed for PBC
                rxij -= round(rxij * boxinv) * box
                ryij -= round(ryij * boxinv) * box
                rijsq = rxij**2 + ryij**2

                if ((rijsq < 10.0) and (rijsq > 0.0)):
                    SR2 = sigsq/rijsq
                    SR6 = SR2 * SR2 * SR2
                    SR12 = SR6 * SR6
                    VIJ  = SR12 - SR6
                    V = V + VIJ         
                    WIJ = VIJ + SR12
                    W = W + WIJ
                    fij = WIJ/rijsq
                    fxij = eps24*fij * rxij
                    fyij = eps24*fij * ryij
                    FXI  = FXI + fxij
                    FYI  = FYI + fyij
                    self.atoms[j].atom.fx -= fxij
                    self.atoms[j].atom.fy -= fyij
            self.atoms[i].atom.fx = FXI
            self.atoms[i].atom.fy = FYI

        V = V * eps4
        return V

    def forces_1(self,i):
    #
    # Simple 1/r potentials
    # i > 0 - Repulsive forces  V = + 1/r
    # i < 0 - Attractive forces V = - 1/r
        n = self.number_atoms
 

        n = self.number_atoms
        #self.md_constantT = 2
        if( i > 0):
            sig1 = 2.5 * self.interaction_strength
            sig2 = 0.0
        else:
            #self.md_constantT = 1
            sig1 = 1.0
            sig2 = 2.1 + 0.1 *(self.interaction_strength-1.0)

        sigma = 0.21
        epslon = 1.0

        sigsq = sigma**2
        eps4 = epslon * 4.0
        eps24 = epslon * 24.0

        box = self.BoxL
        if self.PBC:
            boxinv = 1.0 / box
        else:
            boxinv = 0.0

        V, W = 0.0, 0.0

        for i in range(n):
            self.atoms[i].atom.forces(0.0,0.0)
            self.atoms[i].all_stop = 0

        for i in range(n-1):
            RXI = self.atoms[i].atom.x
            RYI = self.atoms[i].atom.y
            FXI = self.atoms[i].atom.fx
            FYI = self.atoms[i].atom.fy

            for j in range(i+1,n):
                rxij = RXI - self.atoms[j].atom.x
                ryij = RYI - self.atoms[j].atom.y
                ## JLP - needed for PBC
                rxij -= round(rxij * boxinv) * box
                ryij -= round(ryij * boxinv) * box
                rijsq = rxij**2 + ryij**2

                if ((rijsq < 9.0) and (rijsq > 0.0)):
                    SR2 = sigsq/rijsq
                    SR6 = SR2 * SR2 * SR2
                    SR12 = SR6 * SR6
                    VIJ  = sig1*SR12 - sig2*SR6
                    V = V + VIJ         
                    WIJ = VIJ + sig1*SR12
                    W = W + WIJ
                    fij = WIJ/rijsq
                    fxij = eps24*fij * rxij
                    fyij = eps24*fij * ryij
                    FXI  = FXI + fxij
                    FYI  = FYI + fyij
                    self.atoms[j].atom.fx -= fxij
                    self.atoms[j].atom.fy -= fyij

            self.atoms[i].atom.fx = FXI
            self.atoms[i].atom.fy = FYI

        V = V * eps4
        return V
    
class simulation(Frame):
    def __init__(self,master=None):
        Frame.__init__(self,master)
        
    def run_simulation(self,md,gui):
        i = 0
        if not md._running:  md._running = 1
        
        if(md.md_potential == 0):
            text0 = 'Lennard-Jones'
        elif(md.md_potential == 1):
            text0 = 'Repulsive'
        elif(md.md_potential == 2):
            text0 = 'Attractive'
        
        gui.canvas4.delete('text1')
        
        text2 = "Potential: " + text0 +"; "
        text1 = text2 + gui.selection
        gui.canvas4.create_text(205,15,text=text1,justify=LEFT,font='Helvetica 12',fill="#993366",tag='text1')
        
        while md._running == 1:
            
            md.loop(i,0.005*md.timestep)
            
            if(md.static == 2*md.number_atoms):
                for i in range(md.number_atoms):
                    md.atoms[i].atom.velocity(0.0,0.0)
                self.static = 0
                md._running = 0
                
            i += 1
            if((i%20)==0):
                gui.md_move(md)

    def stop_simulation(self,md):
        if md._running:
            md._running = 0

    def quit_application(self,md,g):
        if md._running:
            md._running = 0
        g.tk.quit()

def setup_timestep(step,md):
    md.timestep = step
    
def CueBall_go(sim,g,md):
    mu = md.atoms[0].atom.mu
    if(mu > 0.0):
        md.set_friction(mu)
    else:
        md.friction = 0
    sim.run_simulation(md,g)

def setup(n,BoxL,PBC,h,g):
    mu = h.friction
    if(g.gui_initial == 1):
        if(h.BoxL != BoxL):
            h.MD_cleanup()
            h.setup(n,BoxL,PBC)
        elif(n > h.number_atoms):
            if(h.loop1 == 1):
                #Cation - Warning
                h.add(n)
            else:
                h.MD_cleanup()
                h.setup(n,BoxL,PBC)
        elif (h.PBC != PBC):
            h.MD_cleanup()
            h.setup(n,BoxL,PBC)
        else:
            g.md_delete(h,n)
    else:
        h.MD_cleanup()
        h.setup(n,BoxL,PBC)

    h.MD_reset()
    g.md_cleanup()
    g.setup(h)

    if(g.gui_initial == 0):
        g.gui_initial = 1
##    else:
##        g.md_move(h)
##        g.tk.update()
        
def restart(n,box,PBC,h,g):
    if(h.md_potential == 3):
        del h.atoms[0]
        h.number_atoms -= 1
#    h.md_potential = 0
    h.MD_reset()
    setup(n,box,PBC,h,g)
    h.md_anlysis = 1

def read_gro_file_data(filename):
    f = open(filename,'r')
    line = f.readline()
    try:
        time = float(line[line.rindex('t=')+2:].split()[0])
    except:
        time = 0.0
    try:
        natoms = int(f.readline())
    except:
        natoms = 0

    names = []
    pos_x = []
    pos_y = []
    vel_x = []
    vel_y = []
    color = []
    trace = []

    for x in range(natoms):
        line = f.readline()
        names.append(line[10:15].lstrip().rstrip())
        pos_x.append(float(line[20:28]))
        pos_y.append(float(line[28:36]))
        vel_x.append(float(line[44:52]))
        vel_y.append(float(line[52:60]))
        try:
            color.append(line[69:76])
        except:
            color.append('#FFCC00')

        try:
            trace.append(bool(int(line[76:78])))
        except:
            trace.append(False)

    line = f.readline()

    try:
        box = float(line.split()[0])
    except:
        box = max(pos_x + pos_y) - min(pos_x + pos_y)

    f.close()

    return (time,natoms,box,names,pos_x,pos_y,vel_x,vel_y,color,trace)

def write_gro_file_data(filename,md):
    f = open(filename,'w')
    print("2DMD Simulation, t=%f"%(md.time), file=f)
    print("%5d"%(md.number_atoms), file=f)
    for x in range(md.number_atoms):
        print("%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f%8s%2d"%(0,"None",md.atoms[x].name,x,
                                                                      md.atoms[x].atom.x,
                                                                      md.atoms[x].atom.y,
                                                                      0.0,
                                                                      md.atoms[x].atom.vx,
                                                                      md.atoms[x].atom.vy,
                                                                      0.0,
                                                                      md.atoms[x].color,
                                                                      md.atoms[x].trace), file=f)
    print("%10.5f%10.5f%10.5f"%(md.BoxL,0.0,0.0), file=f)
    f.close()

def load(PBC,h,g):
    try:
        filename = askopenfilename(filetypes=[('GRO Files','*.gro')],multiple=False,initialfile='2D_md.gro')
    except:
        return()

    if filename == '':
        return()

    time, natoms, box, names, pos_x, pos_y, vel_x, vel_y, color, trace = read_gro_file_data(filename)

    h.MD_cleanup()
    h.setup(natoms,box,PBC)
    h.MD_reset()

    ## Opposite here - if we are PBC, then don't reimage
    ## Note that loading a PBC run into a non-PBC condition can
    ## be catastrophic!!!
    if not PBC:
        boxinv = 1.0 / box
    else:
        boxinv = 0.0

    for i in range(h.number_atoms):
        local_x = pos_x[i]
        local_x -= round(local_x * boxinv) * box
        local_y = pos_y[i]
        local_y -= round(local_y * boxinv) * box
        h.atoms[i].atom.x = local_x
        h.atoms[i].atom.y = local_y
        h.atoms[i].atom.last_x = local_x
        h.atoms[i].atom.last_y = local_y
        h.atoms[i].atom.vx = vel_x[i]
        h.atoms[i].atom.vy = vel_y[i]
        h.atoms[i].name = names[i]
        h.atoms[i].color = color[i]
        h.atoms[i].trace = trace[i]

    h.time = time
    g.md_cleanup()
    g.setup(h)
    scale.set(h.number_atoms)
    scale1.set(h.BoxL)

    if(g.gui_initial == 0):
        g.gui_initial = 1

def save(h):
    try:
        filename = asksaveasfilename(filetypes=[('GRO Files','*.gro')],initialfile='2D_md.gro')
    except:
        return()

    if filename == '':
        return()

    write_gro_file_data(filename,h)

def seedone(i,h,g):
    g.selection = i
    print(i)
    h._running = 0
        
    g.md_cleanup()
    g.setup(h)
    
    scale7.config(fg='#EEEEEE',state='disabled')
    
    h.md_potential = 0
    scale.config(fg='black',state='normal')
    scale1.config(fg='black',state='normal')
    
    if(h.md_constantT != 0):
        scale7.config(fg='black',state='normal')
        scale7.set(1.0)

    if(h.md_potential == 0):
        text0 = 'Lennard-Jones'
    elif(h.md_potential == 1):
        text0 = 'Repulsive'
    elif(h.md_potential == 2):
        text0 = 'Attractive'

    g.canvas4.delete('text1')

    text2 = "Potential: " + text0 +"; "
    text1 = text2 + g.selection
    g.canvas4.create_text(205,15,text=text1,justify=LEFT,font='Helvetica 12',fill="#993366",tag='text1') 

def set_dyna(i,h,g):
    h.md_constantT = 0
    if(i == 'Microcanonical'):
        h.md_constantT = 0
    elif(i == 'Scaling'):
        h.md_constantT = 1
    elif(i == 'Berendsen'):
        h.md_constantT = 2

    print(i)
    
    if(h.md_constantT != 0):
        scale7.config(fg='black',state='normal')
        scale7.set(1.0)
    else:
        scale7.config(fg='#EEEEEE',state='disabled')


def set_analysis(i,h,g):
    h.md_anlysis = 0
    if(i == 'Speed'):
        h.md_anlysis = 0
    elif(i == 'Velo'):
        h.md_anlysis = 1
    else:
        h.md_anlysis = 0

    print(i)
        
def seedone1(i,h,g):

    if(i == 'Lennard-Jones'):
        h.md_potential = 0
        h.interaction_strength = 1.0
        scale5.config(fg='black', state='normal')
        scale5.set(1.0)
        text0 = 'Lennard-Jones'
    elif(i== 'Repulsive'):
        if(h.md_potential == 2):
            for i in range(h.number_atoms):
                h.atoms[i].atom.velocity(0.0,0.0)
        h.md_potential = 1
        if(g.selection == 'MD Simulation'):
            h.md_constantT = 2
        scale5.config(fg='black',state='normal')
        scale5.set(1.0)
        text0 = 'Repulsive'
    elif(i== 'Attractive'):
        h.md_potential = 2
        if(g.selection == 'MD Simulation'):
            h.md_constantT = 2
        scale5.config(fg='black',state='normal')
        scale5.set(1.0)
        text0 = 'Attractive'
    else:
        h.md_potential = 0
        h.interaction_strength = 1.0
        scale5.config(fg='#EEEEEE',state='disabled')
        text0 = 'Lennard-Jonse'
        
    g.canvas4.delete('text1')
    text2 = "Potential: " + text0 +"; "
    text1 = text2 + g.selection
    g.canvas4.create_text(205,15,text=text1,justify=LEFT,font='Helvetica 12',fill="#993366",tag='text1')
    g.tk.update()

def setup_mass(Temp,h):
    if(g.gui_initial == 1):
        n = h.number_atoms
        for i in range(n):
            h.atoms[i].atom.mass = Temp * 20.0

    h.MD_reset()

def setup_mass2(Temp,h):
    if(g.gui_initial == 1):
        n = h.number_atoms
        for i in range(n):
            h.atoms[i].atom.mass = Temp * 1.0

    h.MD_reset()
    
def setup_temp(Temp,h):
    h.T_desire = Temp
    h.MD_reset()

def setup_trace(trace,h):
    h.trace = trace
    h.MD_reset()
    
def Interaction_strength(h,force):
    h.interaction_strength = force
    h.MD_reset()
    
def Print(g,i):
    postscript = g.canvas.postscript(fontmap='fontmap')
    filename = 'histogram_'+str(i)+'.ps'
    fd = open(filename,'w')
    fd.write(postscript)
    fd.close()
    i += 1


# Call main when run as script
if __name__ == '__main__':

    global scale, scale1, scale5, scale7
    global print_num, interaction_strength
    print_num = 0
    
    h = MD()
    g = GUI()
    
    natoms = IntVar()
    ndegree = IntVar()
    BoxL = DoubleVar()
    PBC = IntVar()
    Trace = IntVar()
    friction = DoubleVar()
    speed = DoubleVar()

    scale = Scale(g.frame2,label='Number of atoms',orient='horizontal',font='Helvetica 9',
                  variable=natoms, from_=1, to=100,
                  command=(lambda i=natoms.get(): setup(natoms.get(),BoxL.get(),PBC.get(),h,g)))
    scale.set(40)
    scale.pack(side=TOP)
    
    scale1 = Scale(g.frame2,label='Box Size',orient='horizontal',font='Helvetica 9',
                   variable=BoxL, from_=3, to=15, resolution=0.5,
                   command=(lambda i=BoxL: setup(natoms.get(),BoxL.get(),PBC.get(),h,g)))
    scale1.set(4.0)
    scale1.pack(side=TOP)

    ndegree = IntVar()
    friction = DoubleVar()
    speed = DoubleVar()
    interaction_strength=DoubleVar()
    time_step = DoubleVar()
    mass_p = DoubleVar()
     
    scale5 = Scale(g.frame2,label='Interaction strength',orient='horizontal',font='Helvetica 8',
                   variable=interaction_strength, from_=0.2, to=10.0, resolution = 0.2,
                   command=(lambda i=interaction_strength:  Interaction_strength(h,interaction_strength.get())))

#    scale8 = Scale(g.frame2,label='Mass (x20)',orient='horizontal',
#                   variable=mass_p, from_=0.2, to=7,resolution=0.1,font='Helvetica 9',
#                   command=(lambda i=mass_p: setup_mass(mass_p.get(),h)))

#
# fix to reflect the correct mass on the bar
    scale8 = Scale(g.frame2,label='Mass',orient='horizontal',
                   variable=mass_p, from_=10, to=60,resolution=10,font='Helvetica 9',
                   command=(lambda i=mass_p: setup_mass2(mass_p.get(),h)))
    
    scale6 = Scale(g.frame2,label='Time Step',orient='horizontal',
                   variable=time_step, from_=0.2, to=3,resolution=0.2,font='Helvetica 9',
                   command=(lambda i=time_step: setup_timestep(time_step.get(),h)))

    Temp_desired = DoubleVar()
    scale7 = Scale(g.frame2,label='Temperature',orient='horizontal',
                   variable=Temp_desired, from_=0.1, to= 5,resolution=0.1,font='Helvetica 9',
                   command=(lambda i=Temp_desired: setup_temp(Temp_desired.get(),h)))

    scale5.set(1.0)
    scale5.pack(side=TOP)
    scale8.set(20.0)
    scale8.pack(side=TOP)    
    scale7.set(1.0)
    scale7.pack(side=TOP)
    scale6.set(2.0)
    scale6.pack(side=TOP)

    check1 = Checkbutton(g.frame2,text="Periodic Boundary",font='Helvetica 9',variable=PBC,
                         command=(lambda : setup(natoms.get(),BoxL.get(),PBC.get(),h,g)))
    check1.deselect()
    check1.pack(side=TOP)

    check1 = Checkbutton(g.frame2,text="Trace Selected",font='Helvetica 9',variable=Trace,
                         command=(lambda : setup_trace(Trace.get(),h)))
    check1.select()
    check1.pack(side=TOP)

    sim = simulation(g.tk)
    
    run = Button(g.frame2,text='STOP Simulation',command=(lambda : sim.stop_simulation(h)),
                 font="Helvetica 10 bold italic", fg='gold', bg='darkgrey', width=15)
    run.pack(side=BOTTOM)
    
    run1 = Button(g.frame2,text='RUN Simulation',command=(lambda : sim.run_simulation(h,g)),
                  font="Helvetica 10 bold italic", fg='gold', bg='darkgrey', width=15)
    run1.pack(side=BOTTOM)
    
    top = Menu(g.tk)
    g.tk.config(menu=top)
    
##    file4 = Menu(top, font='Helvetica 9 italic')
    file5 = Menu(top, font='Helvetica 9')
##    file2 = Menu(top, font='Helvetica 9 italic')
##    file3 = Menu(top, font='Helvetica 9 italic')
    
    cascad1 = Menu(top, tearoff=0, font='Helvetica 9')
    top.add_cascade(label='Run',menu=cascad1)
    titre1 = Menu(cascad1, tearoff=0, font='Helvetica 9')
    cascad1.add_command(label='Restart', command=lambda : restart(natoms.get(),BoxL.get(),PBC.get(),h,g))
    cascad1.add_command(label='Load', command=lambda : load(PBC.get(),h,g))
    cascad1.add_command(label='Save', command=lambda : save(h))
    cascad1.add_command(label='Quit', command=lambda : g.Quit(h))    

#
# Potential Selections
#
    menu5 = ['Lennard-Jones','Repulsive','Attractive']
    for j in menu5:
        file5.add_command(label=j, command=(lambda i=j: seedone1(i,h,g)))
    top.add_cascade(label='Potential Form', menu=file5)

#
# Simulation System Setup
#
    cascad2 = Menu(top, tearoff=0, font='Helvetica 9')    
    top.add_cascade(label="Dynamics Setup",menu=cascad2)
    tire2 = Menu(cascad2, tearoff=0, font='Helvetica 9 italic')
    cascad2.add_command(label='Constant Energy', command=lambda : set_dyna('Microcanonical',h,g))
    color2 = Menu(cascad2, font='Helvetica 9 italic')
    cascad2.add_cascade(label='Constant Temp', menu=color2)
    color2.add_command(label="Simple Rescaling", state='disabled')
    color2.add_command(label="Berendsen Thermostat", command=lambda : set_dyna('Berendsen',h,g))

#
# Analysis Menu Construction
#
    cascad3 = Menu(top, tearoff=0, font='Helvetica 9')    
    top.add_cascade(label="Analysis",menu=cascad3)
    color3 = Menu(cascad3, font='Helvetica 9 italic')
    cascad3.add_cascade(label='Kinetics', menu=color3)
    color3.add_command(label="Speed Distribution", command=lambda : set_analysis('Speed',h,g))
    color3.add_command(label="Velocity Distribution (x-component)", command=lambda : set_analysis('Velo',h,g))

    
    mainloop()

