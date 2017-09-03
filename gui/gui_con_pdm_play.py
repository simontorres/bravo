import sys
import matplotlib
matplotlib.use('QT4Agg')
from matplotlib.widgets import Button
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
from matplotlib.widgets import MultiCursor
#import asciitable # hay que instalar
#from astroML.time_series import lomb_scargle #hay que instalar
#from gatspy import datasets, periodic #hay que instalar
#from astroML.plotting import setup_text_plots
#setup_text_plots(fontsize=8, usetex=True)
import os


import argparse
from astropy.stats import LombScargle



def get_args(arguments=None):
    parser = argparse.ArgumentParser(
        description="Interactive tool to search for variable objects")

    parser.add_argument('table',
                        action='store',
                        help="Table of photometry or radial velocity data.")
    args = parser.parse_args(args=arguments)

    return args


def cargar_datos(tabla):
    data = np.genfromtxt(tabla)
    jda = data[:,0]
    maga = data[:,1]
    erra = data[:,2]
    return jda, maga, erra

def calculo_fase(mag, date, er, per, T0):
    fase2, tt2, err, t = [], [], [], []
    for i in range(len(mag)):
        fa2 = ((float(date[i]) * 1.0 - T0) / per) - int((float(date[i]) - T0) / per)
        if fa2 > 0:
            fase2.append(fa2)
            tt2.append(fa2 + 1.0)
            t.append(date[i])
            err.append(er[i])
        else:
            fase2.append(fa2 + 1)
            tt2.append(fa2 + 2)
            t.append(date[i])
            err.append(er[i])
    fase2 = np.array(fase2)
    tt2 = np.array(tt2)
    err = np.array(err)
    re2 = np.concatenate((fase2, tt2), axis=0)
    mag3 = np.concatenate((mag, mag), axis=0)
    t2 = np.concatenate((t, t), axis=0)
    err2 = np.concatenate((err, err), axis=0)

    return re2, mag3, t2, err2

def spline(jda,maga,orden,splineYes=True):
    spl = UnivariateSpline(jda, maga, k=orden)
    return spl,splineYes

class GuiExample(object):

    def __init__(self):
        self.tabla = None
        self.fig = None
        self.ax1 = None
        self.ax2 = None
        self.ax3 = None
        self.ax1_bb = None
        self.line_plot = None
        self.jda = None
        self.maga = None
        self.erra = None
        self.per = None
        self.t0 = None
        self.min_per = 0.1
        self.max_per = 10.0
  #      self.step_per = 80000.0
	self.step_pdm = 50000.0
        self.periodos = None
        self.omega = None
        self.PS = None
        self.model = None
        self.power = None
        self.freqs = None
        self.spl = None
        self.splineYes = None
        self.multi = None
        self.name_star = "V1216 Sco"
        self.auto_per=None

    def __call__(self, tabla, *args, **kwargs):
        self.tabla = tabla
        self.jda, self.maga, self.erra = cargar_datos(self.tabla)
        self.fig, (self.ax1, self.ax2 , self.ax3) = plt.subplots(nrows=3)
        self.fig.subplots_adjust(hspace=0.27, bottom=0.07, top=0.99, left=0.08, right=0.98)
        self.multi = MultiCursor(self.fig.canvas, (self.ax1,self.ax2), color='r', \
                                 lw=.5, horizOn=None, vertOn=True)

        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()

        #jd/mag
        if self.splineYes == True:
            self.spl, self.splineYes = spline(self.jda, self.maga, 5)
            print "Aplicando spline"
            self.maga = self.maga - self.spl(self.jda)
            #self.ax1.plot(self.jda, self.spl(self.jda), 'r--', lw=3)
            self.spl, self.splineYes = spline(self.jda, self.maga, 5)
            self.ax1.plot(self.jda, self.spl(self.jda), 'r--', lw=3)
            self.maga = self.maga - self.spl(self.jda)
        else:
            print "NO aplicando spline"
            self.spl,self.splineYes = spline(self.jda, self.maga,5)
            self.ax1.plot(self.jda, self.spl(self.jda), 'g--', lw=3)



        self.ax1.plot(self.jda, self.maga, 'o', c='black')
        self.ax1.set_ylim(max(self.maga+0.01), min(self.maga)-0.01)
        self.ax1.set_xlim(min(self.jda)-10, max(self.jda)+0.01)
        self.ax1.set_xlabel(r'Time', fontsize=20)
        self.ax1.set_ylabel(r"Magnitud", fontsize=20)

        if self.name_star!=True:
            self.ax1.text(0.03, 0.145, "%s"%self.name_star, ha='left', va='top', \
                          transform=self.ax1.transAxes, fontsize=25,color="red")




	
        #LS astropy
        frequency2, power3 = LombScargle(self.jda, self.maga, self.erra).autopower\
            (minimum_frequency=1 / self.max_per,maximum_frequency=1 / self.min_per)
        self.ax2.plot(1.0 / frequency2, power3, '-', c='cyan', lw=1, zorder=1, \
                      label="PyLS")

        if self.auto_per == True:
            self.step_per = len(frequency2)
        """
        #GLS
        self.periodos=np.linspace(self.min_per, self.max_per, self.step_per)
        self.omega = 2 * np.pi / self.periodos
        self.PS = lomb_scargle(self.jda, self.maga, self.erra, self.omega, generalized=True)
        self.ax2.plot(self.periodos, self.PS, '-', c='black', lw=1, zorder=1,label="GLS")

        #LS
        model = periodic.LombScargle().fit(self.jda, self.maga, self.erra)
        fmin = 1.0 / self.max_per
        fmax = 1.0 / self.min_per
        df = (fmax - fmin) / self.step_per
        self.power = model.score_frequency_grid(fmin, df, self.step_per)
        self.freqs = fmin + df * np.arange(self.step_per)
        self.ax2.plot(1.0/self.freqs, self.power, '-', c='red', lw=1, zorder=1,label="LS")
        self.ax2.legend(fontsize = 'x-large')
        """


        #PDM

        os.system("awk '{print $1,$2,$3,$1,$2}' %s > borrar.dat"%tabla)
        self.tabla="borrar.dat"
	longi=open(self.tabla, 'r').read().count("\n")
        pdm1=float(os.popen("./pdmmm %s %0.1f %0.3f %0.3f 10 5 %0.3f"%(self.tabla,longi,self.min_per,self.max_per,self.step_pdm)).readlines()[0])
	f_11=np.genfromtxt(self.tabla+".pdm")
        Pdm1t,Ppdm1=f_11[:,0],1./f_11[:,1]
        print len(f_11[0]),1./min(f_11[1])

        self.ax2.plot(Ppdm1, Pdm1t, '-', c='green', lw=1, zorder=1,label="PDM",alpha=0.8)
        self.ax2.legend(fontsize = 'x-large')


        self.ax2.set_xlabel(r'Periodo', fontsize=20)
        self.ax2.set_ylabel(r"Power", fontsize=20)


        self.ax2.set_xlim(self.min_per, self.max_per)

        self.ax1_bb = self.ax2.get_position()

        self.fig.canvas.mpl_connect('motion_notify_event', self.on_mouse_over)
        plt.tight_layout()
#        plt.show()


    def on_mouse_over(self, event):
        ax1_x, ax1_y = \
            self.fig.transFigure.inverted().transform((event.x, event.y))

        if self.ax1_bb.contains(ax1_x, ax1_y):
            if self.line_plot is not None:
                try:
                    self.line_plot.remove()
                    self.ax3.relim()
                    
                except:
                    pass
            if event.ydata is not None:
                self.per=event.xdata
                print event.xdata
                self.t0=0
                fasA,magniA,t_A,er_A=calculo_fase(self.maga, self.jda, self.erra, self.per, self.t0)
                self.line_plot, = self.ax3.plot(fasA,magniA,"o",color='k')
                self.ax3.set_xlim(0,2)
                self.ax3.set_ylim(max(magniA+0.01), min(magniA)-0.01)
                self.fig.canvas.draw()
                self.ax3.set_xlabel(r'Fase', fontsize=20)
                self.ax3.set_ylabel(r"Magnitud", fontsize=20)



if __name__ == '__main__':

    args =  get_args()
    # print(args.table)
    gui = GuiExample()
    gui(tabla=args.table)
