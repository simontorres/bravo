import matplotlib
matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
from ccdproc import CCDData


class GuiExample(object):

    def __init__(self):
        self.ccd = None
        self.fig = None
        self.ax1 = None
        self.ax2 = None
        self.ax1_bb = None
        self.line_plot = None

    def __call__(self, *args, **kwargs):
        self.ccd = CCDData.read('/home/simon/data/soar/work/aller/2017-06-09/RED/cfzsto_0111_CuHeAr_G1200M2_slit103.fits')
        self.fig, (self.ax1, self.ax2) = plt.subplots(nrows=2, sharex=True)
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
        self.ax1.imshow(self.ccd.data)

        self.ax1_bb = self.ax1.get_position()

        self.fig.canvas.mpl_connect('motion_notify_event', self.on_mouse_over)
        plt.tight_layout()
        plt.show()

    def on_mouse_over(self, event):
        ax1_x, ax1_y = \
            self.fig.transFigure.inverted().transform((event.x, event.y))

        if self.ax1_bb.contains(ax1_x, ax1_y):
            # print('Se movio')
            # A, B = event.xdata, event.ydata
            if self.line_plot is not None:
                try:
                    self.line_plot.remove()
                    self.ax2.relim()
                except:
                    pass
            if event.ydata is not None:
                data = self.ccd.data[int(event.ydata), :]
                self.line_plot, = self.ax2.plot(data,
                                                color='k')
                self.fig.canvas.draw()













if __name__ == '__main__':
    gui = GuiExample()

    gui()