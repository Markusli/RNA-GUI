from PyQt4 import QtCore, QtGui
import numpy as np
from matplotlibwidget import maplot, scatter
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from Bio import SeqIO
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
        self.combo1.addItems('MazF none,MazF PNK,MazF TAP,MqsR none,MqsR PNK,MqsR TAP,Log none,Log PNK,Log TAP,Stats none,Stats PNK,Stats TAP,Delta3 none,Delta3 PNK,Delta3 TAP'.split(','))

        self.combo2 = QtGui.QComboBox(self.centralwidget)
        self.combo2.addItems('MazF none,MazF PNK,MazF TAP,MqsR none,MqsR PNK,MqsR TAP,Log none,Log PNK,Log TAP,Stats none,Stats PNK,Stats TAP,Delta3 none,Delta3 PNK,Delta3 TAP'.split(','))

        self.combo3 = QtGui.QComboBox(self.centralwidget)
        self.combo3.addItems('5prime,3prime'.split(','))
        self.combo3.currentIndexChanged[str].connect(self.on_combo_prime_change)

        self.combo4 = QtGui.QComboBox(self.centralwidget)
        self.combo4.addItems('16S,23S'.split(','))

        self.combo5 = QtGui.QComboBox(self.centralwidget)
        self.combo5.addItems('Percentages,Absolute values'.split(','))

        self.text = QtGui.QTextEdit(self.centralwidget)
        self.text.setMaximumWidth(400)
        self.cursor = self.text.textCursor()
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
        cbbox.addWidget(self.combo5)

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
        self.combo5.setEnabled(False)
        self.combo5.setCurrentIndex(1)
        items_5 = 'MazF none,MazF PNK,MazF TAP,MqsR none,MqsR PNK,MqsR TAP,Log none,Log PNK,Log TAP,Stats none,Stats PNK,Stats TAP,Delta3 none,Delta3 PNK,Delta3 TAP'
        items_3 = 'MazF PNK-,MazF PNK+,MqsR PNK-,MqsR PNK+,Log PNK-,Log PNK+,Delta3 PNK-,Delta3 PNK+'
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

        name_conv = {'MazF': 'MazF2h',
            'MqsR': 'MqsR2h',
            'Log': 'MG1655log',
            'Stat': 'MG1655stats',
            'Delta3': 'delta3'}

        sample1 = self.combo1.currentText()
        sample2 = self.combo2.currentText()
        prim_select = self.combo3.currentText()
        sub_select = self.combo4.currentIndex()
        value_select = self.combo5.currentIndex()
        sample1 = str(sample1).split(' ')
        sample2 = str(sample2).split(' ')
        prim_select = str(prim_select)
        data_dic = hdf5_gen.get_hdf_data([prim_select, prim_select], [name_conv[sample1[0]], name_conv[sample2[0]]],[sample1[1], sample2[1]],['_input_', 'some'])

        subunit = ['y_pos_16S', 'y_pos_23S']
        nucl_data = ['nucl_data_16S', 'nucl_data_23S']
        subunit_neg = ['y_neg_16S', 'y_neg_23S']
        fastas = ['data/16S_new.fasta', 'data/23S_new.fasta']

        #Gets the index of datapoint (index of X and Y)
        ind = event.ind

        #Retrieves information from lists based on index of the datapoint.
        #The information to be displayed for datapoint has the same index as datapint.
        for fasta in SeqIO.parse(fastas[sub_select], "fasta"):
            fasta = fasta.seq
            pass

        if value_select == 1:
            nucleotide_pos = np.take(data_dic['data1'][nucl_data[sub_select]], ind)
            sample_1_pos_read_count = np.take(data_dic['data1'][subunit[sub_select]], ind)
            sample_2_pos_read_count = np.take(data_dic['data2'][subunit[sub_select]], ind)
            sample_1_neg_read_count = np.take(data_dic['data1'][subunit_neg[sub_select]], ind)
            sample_2_neg_read_count = np.take(data_dic['data2'][subunit_neg[sub_select]], ind)

            for array_ind in range(len(nucleotide_pos)):
                position = int(nucleotide_pos[array_ind])
                sequence = str(fasta[position-3:position] + '<b>' + fasta[position] + '</b>' + fasta[position+1:position+4])
                message = "Primary sequence:<br>{7}<br>Nucleotide position: {0:,} <br>{5} '+' strand: {1:,}<br>{6} '+' strand {2:,}\
                                <br>{5} '-' strand: {3:,}<br>{6} '-' strand {4:,}<br>======================="\
                                .format(position,
                                        sample_1_pos_read_count[array_ind],
                                        sample_2_pos_read_count[array_ind],
                                        sample_1_neg_read_count[array_ind],
                                        sample_2_neg_read_count[array_ind],
                                        ' '.join(sample1),
                                        ' '.join(sample2),
                                        sequence)

                self.text.setHtml(message)

        elif value_select == 0:

            hundred_1 = 0
            hundred_2 = 0
            for i in range(-19,22):
                index = data_dic['data1'][nucl_data[sub_select]].index(float(i))
                hundred_1 += data_dic['data1'][subunit[sub_select]][index]
                hundred_2 += data_dic['data2'][subunit[sub_select]][index]
            heights_1 = [(x / hundred_1 * 100) for x in data_dic['data1'][subunit[sub_select]]]
            heights_2 = [(x / hundred_2 * 100) for x in data_dic['data2'][subunit[sub_select]]]
            #Gets the index of datapoint (index of X and Y)
            ind = event.ind


            #Retrieves information from lists based on index of the datapoint.
            #The information to be displayed for datapoint has the same index as datapint.

            nucleotide_pos = np.take(data_dic['data1'][nucl_data[sub_select]], ind)
            sample_1_pos_read_count = np.take(heights_1, ind)
            sample_2_pos_read_count = np.take(heights_2, ind)
            sample_1_pos_abs_count = np.take(data_dic['data1'][subunit[sub_select]], ind)
            sample_2_pos_abs_count = np.take(data_dic['data2'][subunit[sub_select]], ind)
            #sample_1_neg_read_count = np.take(data_dic['data1'][subunit_neg[sub_select]], ind)
            #sample_2_neg_read_count = np.take(data_dic['data2'][subunit_neg[sub_select]], ind)

            for array_ind in range(len(nucleotide_pos)):
                position = int(nucleotide_pos[array_ind])
                sequence = str(fasta[position-3:position] + '<b>' + fasta[position] + '</b>' + fasta[position+1:position+4])
                message = "Primary sequence:<br>{7}<br>Nucleotide position: {0:,}<br>Percentages:<br>{3}: {1}<br>{4}:\
                            {2}<br>Absolute values:<br>{3}: {5:,}<br>{4}: {6:,}<br>======================="\
                            .format(int(nucleotide_pos[array_ind]),
                                    round(sample_1_pos_read_count[array_ind],4),
                                    round(sample_2_pos_read_count[array_ind], 4),
                                    ' '.join(sample1), ' '.join(sample2),
                                    int(sample_1_pos_abs_count[array_ind]),
                                    int(sample_2_pos_abs_count[array_ind]),
                                    sequence)

                self.text.setHtml(message)








