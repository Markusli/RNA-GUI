import sys
from matplot import *
import random
import math
from PyQt4 import QtCore, QtGui
from mpltools import style

style.use('ggplot')

sys.path.append('/home/markus/git/plotter')
import hdf5_gen

class GUIForm(QtGui.QDialog):

    def __init__(self, parent=None):

        QtGui.QWidget.__init__(self,parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        sys.stdout = OutLog(self.ui.text, sys.stdout)
        sys.stderr = OutLog(self.ui.text, sys.stderr, QtGui.QColor(255,0,0) )

        QtCore.QObject.connect(self.ui.plotButton, QtCore.SIGNAL('clicked()'), self.PlotFunc)
        QtCore.QObject.connect(self.ui.clearButton, QtCore.SIGNAL('clicked()'), self.ui.text.clear)


    def PlotFunc(self):
        index1 = self.ui.combo1.currentText()
        index2 = self.ui.combo2.currentText()
        index3 = self.ui.combo3.currentText()
        index4 = self.ui.combo4.currentIndex()
        index1 = str(index1).split(' ')
        index2 = str(index2).split(' ')
        index3 = str(index3)

        data_dic = hdf5_gen.get_hdf_data([index3, index3], [index1[0], index2[0]],[index1[1], index2[1]],['_input_', 'some'])


        subunit = ['y_pos_16S', 'y_pos_23S']
        subunitc = ['colour_16S', 'colour_23S']
        nucl_data = ['nucl_data_16S', 'nucl_data_23S']
        subunit_neg = ['y_neg_16S', 'y_neg_23S']
        plotname = ['16S Scatterplot', '23S Scatterplot']

        MA_X_16S = [(math.log(float(y1), 2) + math.log(float(y2), 2))/2 for y1, y2 in zip(data_dic['data1'][subunit[index4]], data_dic['data2'][subunit[index4]])]
        MA_Y_16S = [math.log((float(y1)/y2), 2) for y1, y2 in zip(data_dic['data1'][subunit[index4]], data_dic['data2'][subunit[index4]])]
        self.ui.widget1.canvas.ax.clear()
        self.ui.widget1.canvas.ax.set_ylabel('M', fontsize=20)
        self.ui.widget1.canvas.ax.set_xlabel('A', fontsize=20)
        self.ui.widget1.canvas.ax.scatter(MA_X_16S, MA_Y_16S, alpha=0.5, c=data_dic["data1"][subunitc[index4]], linewidths=( 0, 0, 0), picker=True)
        self.ui.widget1.canvas.fig.tight_layout()
        self.ui.widget1.canvas.draw()
        self.ui.widget1.canvas.mpl_connect('pick_event', lambda event: onpick(event, data_dic, index4, index1, index2))

        self.ui.widget2.canvas.ax.clear()
        self.ui.widget2.canvas.ax.set_ylabel('Read Counts', fontsize=14)
        self.ui.widget2.canvas.ax.set_xlabel('Nucleotide Position', fontsize=14)
        self.ui.widget2.canvas.ax.set_title(plotname[index4], fontsize=14)
        for data_key, data_nt in data_dic.items():
            self.ui.widget2.canvas.ax.scatter(data_nt[nucl_data[index4]], data_nt[subunit[index4]], alpha=0.5, c=data_nt[subunitc[index4]], linewidths=( 0, 0, 0), picker=True, marker = data_nt['symbol'])
            self.ui.widget2.canvas.ax.scatter(data_nt[nucl_data[index4]], [-1 * data for data in data_nt[subunit_neg[index4 ]]], alpha=0.5, c=data_nt[subunitc[index4]], linewidths=( 0, 0, 0), picker=True, marker = data_nt['symbol'])
        self.ui.widget2.canvas.fig.tight_layout()
        self.ui.widget2.canvas.draw()
        self.ui.widget2.canvas.mpl_connect('pick_event', lambda event: onpick(event, data_dic, index4, index1, index2))


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = GUIForm()
    sys.exit(app.exec_())