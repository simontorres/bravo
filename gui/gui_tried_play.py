import matplotlib
matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
import numpy as np
import asciitable
from astroML.time_series import lomb_scargle
from gatspy import datasets, periodic

def cargar_datos(tabla):
    data = asciitable.read(tabla)
    jda = data['col1']
    maga = data['col2']
    erra = data['col3']
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

class GuiExample(object):

    def __init__(self):
        self.tabla = None
        self.fig = None
        self.ax1 = None
        self.ax2 = None
        self.ax1_bb = None
        self.line_plot = None
        self.jda = None
        self.maga = None
        self.erra = None
        self.per = None
        self.t0 = None
        self.min_per = 0.1
        self.max_per = 10.0
        self.step_per = 10000.0
        self.periodos = None
        self.omega = None
        self.PS = None
        self.model = None
        self.power = None
        self.freqs = None

    def __call__(self, *args, **kwargs):
        self.tabla = "V1216Sco-165458-4356.5.asas.pdm0"
        self.jda, self.maga, self.erra = cargar_datos(self.tabla)
        self.fig, (self.ax1, self.ax2) = plt.subplots(nrows=2)
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
        #GLS
        self.periodos=np.linspace(self.min_per, self.max_per, self.step_per)
        self.omega = 2 * np.pi / self.periodos
        self.PS = lomb_scargle(self.jda, self.maga, self.erra, self.omega, generalized=True)
        self.ax1.plot(self.periodos, self.PS, '-', c='black', lw=1, zorder=1)

        #LS
        model = periodic.LombScargle().fit(self.jda, self.maga, self.erra)
        fmin = 1.0 / self.max_per
        fmax = 1.0 / self.min_per
        df = (fmax - fmin) / self.step_per
        self.power = model.score_frequency_grid(fmin, df, self.step_per)
        self.freqs = fmin + df * np.arange(self.step_per)
        self.ax1.plot(1.0/self.freqs, self.power, '-', c='red', lw=1, zorder=1)

        #PDM


        self.ax1_bb = self.ax1.get_position()

        self.fig.canvas.mpl_connect('motion_notify_event', self.on_mouse_over)
        plt.tight_layout()
        plt.show()

    def on_mouse_over(self, event):
        ax1_x, ax1_y = \
            self.fig.transFigure.inverted().transform((event.x, event.y))

        if self.ax1_bb.contains(ax1_x, ax1_y):
            #print('Se movio')
            # A, B = event.xdata, event.ydata
            if self.line_plot is not None:
                try:
                    self.line_plot.remove()
                    self.ax2.relim()
                except:
                    pass
            if event.ydata is not None:
                self.per=event.xdata
                print event.xdata
                self.t0=0
                fasA,magniA,t_A,er_A=calculo_fase(self.maga, self.jda, self.erra, self.per, self.t0)
                self.line_plot, = self.ax2.plot(fasA,magniA,"o",color='k')
                self.ax2.set_xlim(0,2)
                self.ax2.set_ylim(max(magniA+0.01), min(magniA)-0.01)
                self.fig.canvas.draw()


if __name__ == '__main__':
    gui = GuiExample()
    gui()