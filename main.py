import sys
from matplot import *
import random
import math
from PyQt4 import QtCore, QtGui

sys.path.append('/home/markus/git/plotter')
import hdf5_gen

class GUIForm(QtGui.QDialog):

    def __init__(self, parent=None):

        QtGui.QWidget.__init__(self,parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)



        QtCore.QObject.connect(self.ui.plotButton, QtCore.SIGNAL('clicked()'), self.PlotFunc)


    def PlotFunc(self):
        index1 = self.ui.combo1.currentText()
        index2 = self.ui.combo2.currentText()
        index1 = str(index1).split(' ')
        index2 = str(index2).split(' ')
        data_dic = hdf5_gen.get_hdf_data(['5prim', '5prim'], [index1[0], index2[0]],[index1[1], index2[1]],['_input_', 'some'])
        MA_X_16S = [(math.log(float(y1), 2) + math.log(float(y2), 2))/2 for y1, y2 in zip(data_dic['data1']['y_pos_16S'], data_dic['data2']['y_pos_16S'])]
        MA_Y_16S = [math.log((float(y1)/y2), 2) for y1, y2 in zip(data_dic['data1']['y_pos_16S'], data_dic['data2']['y_pos_16S'])]
        self.ui.widget1.canvas.ax.clear()
        self.ui.widget1.canvas.ax.set_ylabel('M', fontsize=20)
        self.ui.widget1.canvas.ax.set_xlabel('A', fontsize=20)
        self.ui.widget1.canvas.ax.scatter(MA_X_16S, MA_Y_16S, alpha=0.5, c=data_dic["data1"]['colour_16S'], linewidths=( 0, 0, 0), picker=True)
        self.ui.widget1.canvas.draw()

        self.ui.widget2.canvas.ax.clear()
        self.ui.widget2.canvas.ax.set_ylabel('Read Counts', fontsize=14)
        self.ui.widget2.canvas.ax.set_xlabel('Nucleotide Position', fontsize=14)
        self.ui.widget2.canvas.ax.set_title('16S Scatterplot', fontsize=14)
        for data_key, data_nt in data_dic.items():
            self.ui.widget2.canvas.ax.scatter(data_nt['nucl_data_16S'], data_nt['y_pos_16S'], alpha=0.5, c=data_nt['colour_16S'], linewidths=( 0, 0, 0), picker=True, marker = data_nt['symbol'])
            self.ui.widget2.canvas.ax.scatter(data_nt['nucl_data_16S'], [-1 * data for data in data_nt['y_neg_16S']], alpha=0.5, c=data_nt['colour_16S'], linewidths=( 0, 0, 0), picker=True, marker = data_nt['symbol'])
        self.ui.widget2.canvas.draw()


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = GUIForm()
    myapp.show()
    sys.exit(app.exec_())