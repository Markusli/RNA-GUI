from PyQt4 import QtCore, QtGui
import numpy as np
from matplotlibwidget import maplot, scatter
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import hdf5_gen


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(QtGui.QWidget):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))

        self.plotButton = QtGui.QPushButton(self.centralwidget)
        self.plotButton.setObjectName(_fromUtf8("plotButton"))

        self.widget1 = maplot(self.centralwidget)
        self.widget1.setObjectName(_fromUtf8("MA plot"))
        self.toolbar1 = NavigationToolbar(self.widget1.canvas, self)
        self.widget1.canvas.mpl_connect('pick_event', lambda event: self.onpick(event))


        self.widget2 = scatter(self.centralwidget)
        self.widget2.setObjectName(_fromUtf8("Scatter"))
        self.toolbar2 = NavigationToolbar(self.widget2.canvas, self)
        self.widget2.canvas.mpl_connect('pick_event', lambda event: self.onpick(event))


        QtCore.QMetaObject.connectSlotsByName(MainWindow)

        self.combo1 = QtGui.QComboBox(self.centralwidget)
        self.combo1.addItems('MazF2h none,MazF2h PNK,MazF2h TAP,MqsR2h none,MqsR2h PNK,MqsR2h TAP,MG1655log none,MG1655log PNK,MG1655log TAP,delta3 none,delta3 PNK,delta3 TAP'.split(','))

        self.combo2 = QtGui.QComboBox(self.centralwidget)
        self.combo2.addItems('MazF2h none,MazF2h PNK,MazF2h TAP,MqsR2h none,MqsR2h PNK,MqsR2h TAP,MG1655log none,MG1655log PNK,MG1655log TAP,delta3 none,delta3 PNK,delta3 TAP'.split(','))

        self.combo3 = QtGui.QComboBox(self.centralwidget)
        self.combo3.addItems('5prime,3primeok,3primeshady'.split(','))
        self.combo3.currentIndexChanged[str].connect(self.on_combo_prime_change)

        self.combo4 = QtGui.QComboBox(self.centralwidget)
        self.combo4.addItems('16S,23S'.split(','))

        self.text = QtGui.QPlainTextEdit(self.centralwidget)
        self.text.setMaximumSize(400,400)
        self.clearButton = QtGui.QPushButton(self.centralwidget)
        self.clearButton.setObjectName(_fromUtf8("clearButton"))
        self.retranslateUi(MainWindow)

        pbox = QtGui.QVBoxLayout()
        pbox.setSpacing(1)
        pbox.addWidget(self.toolbar1)
        pbox.addWidget(self.widget1)
        pbox.addWidget(self.toolbar2)
        pbox.addWidget(self.widget2)

        tbbox = QtGui.QHBoxLayout()
        tbbox.addWidget(self.plotButton)
        tbbox.addWidget(self.clearButton)

        cbbox = QtGui.QVBoxLayout()
        cbbox.addWidget(self.combo1)
        cbbox.addWidget(self.combo2)
        cbbox.addWidget(self.combo3)
        cbbox.addWidget(self.combo4)

        bbox = QtGui.QVBoxLayout()
        bbox.addLayout(cbbox)
        bbox.addWidget(self.text)
        bbox.addLayout(tbbox)

        abox = QtGui.QHBoxLayout()
        abox.addLayout(bbox)
        abox.addLayout(pbox)

        self.setLayout(abox)

        self.show()



    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "RNA-Tools", None))
        self.plotButton.setText(_translate("MainWindow", "Plot", None))
        self.clearButton.setText(_translate("MainWindow", "Clear Text", None))

    def on_combo_prime_change(self, index):
        items_5 = 'MazF2h none,MazF2h PNK,MazF2h TAP,MqsR2h none,MqsR2h PNK,MqsR2h TAP,MG1655log none,MG1655log PNK,MG1655log TAP,delta3 none,delta3 PNK,delta3 TAP'
        items_3 = 'MazF2h PNK-,MazF2h PNK+,MqsR2h PNK-,MqsR2h PNK+,MG1655log PNK-,MG1655log PNK+'
        text = self.combo3.currentText()
        if text == '5prime':
            self.combo2.clear()
            self.combo2.addItems(items_5.split(','))
            self.combo1.clear()
            self.combo1.addItems(items_5.split(','))
        else:
            self.combo2.clear()
            self.combo2.addItems(items_3.split(','))
            self.combo1.clear()
            self.combo1.addItems(items_3.split(','))

    def onpick(self, event):

        index1 = self.combo1.currentText()
        index2 = self.combo2.currentText()
        index3 = self.combo3.currentText()
        index4 = self.combo4.currentIndex()
        index1 = str(index1).split(' ')
        index2 = str(index2).split(' ')
        index3 = str(index3)
        data_dic = hdf5_gen.get_hdf_data([index3, index3], [index1[0], index2[0]],[index1[1], index2[1]],['_input_', 'some'])

        #Gets the index of datapoint (index of X and Y)
        ind = event.ind

        subunit = ['y_pos_16S', 'y_pos_23S']
        nucl_data = ['nucl_data_16S', 'nucl_data_23S']
        subunit_neg = ['y_neg_16S', 'y_neg_23S']
        #Retrieves information from lists based on index of the datapoint.
        #The information to be displayed for datapoint has the same index as datapint.

        nucleotide_pos = np.take(data_dic['data1'][nucl_data[index4]], ind)
        sample_1_pos_read_count = np.take(data_dic['data1'][subunit[index4]], ind)
        sample_2_pos_read_count = np.take(data_dic['data2'][subunit[index4]], ind)
        sample_1_neg_read_count = np.take(data_dic['data1'][subunit_neg[index4]], ind)
        sample_2_neg_read_count = np.take(data_dic['data2'][subunit_neg[index4]], ind)

        for array_ind in range(len(nucleotide_pos)):
            print ("Nucleotide position: {0} \n{5} pos: {1}\n{6} pos {2}\
                            \n{5} neg: {3}\n{6} neg {4}".format(nucleotide_pos[array_ind], sample_1_pos_read_count[array_ind],
                                                             sample_2_pos_read_count[array_ind], sample_1_neg_read_count[array_ind],
                                                             sample_2_neg_read_count[array_ind], ' '.join(index1), ' '.join(index2)))

class OutLog:
    def __init__(self, edit, out=None, color=None):

        self.edit = edit
        self.out = None
        self.color = color

    def write(self, m):
        if self.color:
            tc = self.edit.textColor()
            self.edit.setTextColor(self.color)

        self.edit.moveCursor(QtGui.QTextCursor.Start)
        self.edit.insertPlainText( m )

        if self.color:
            self.edit.setTextColor(tc)

        if self.out:
            self.out.write(m)







def testing(x):
    print x
