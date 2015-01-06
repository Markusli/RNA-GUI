from PySide import QtCore, QtGui
import numpy as np
import math
import os
import csv
from matplotlibwidget import maplot, scatter
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.patches as mpatches
import matplotlib.collections as collections
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



items3prime = 'MazF PNK-,MazF PNK+,MqsR PNK-,MqsR PNK+,Log PNK-,Log PNK+,Stat PNK-,Stat PNK+,Delta3 PNK-,Delta3 PNK+'
items5prime = 'MazF none,MazF PNK,MazF TAP,MqsR none,MqsR PNK,MqsR TAP,Log none,Log PNK,Log TAP,Stat none,Stat PNK,Stat TAP,Delta3 none,Delta3 PNK,Delta3 TAP'

class Plotter(QtGui.QWidget):
    def setupUi(self):

        self.plotButton = QtGui.QPushButton('Plot', self)
        QtCore.QObject.connect(self.plotButton, QtCore.SIGNAL('clicked()'), self.PlotFunc)

        self.widget1 = maplot(self)
        self.toolbar1 = NavigationToolbar(self.widget1.canvas, self)
        self.widget1.canvas.mpl_connect('pick_event', lambda event: self.onpick(event))

        self.widget2 = scatter(self)
        self.widget2.setObjectName(_fromUtf8("Scatter"))
        self.toolbar2 = NavigationToolbar(self.widget2.canvas, self)
        self.widget2.canvas.mpl_connect('pick_event', lambda event: self.onpick(event))

        self.proc1combo = QtGui.QComboBox(self)
        self.proc1combo.addItems(items5prime.split(','))

        self.proc2combo = QtGui.QComboBox(self)
        self.proc2combo.addItems(items5prime.split(','))

        self.primcombo = QtGui.QComboBox(self)
        self.primcombo.addItems('5prime,3prime'.split(','))
        self.primcombo.currentIndexChanged[str].connect(self.on_combo_prime_change)

        self.subcombo = QtGui.QComboBox(self)
        self.subcombo.addItems('16S,23S'.split(','))

        self.valcombo = QtGui.QComboBox(self)
        self.valcombo.addItems('Percentages,Absolute values'.split(','))

        self.text = QtGui.QTextEdit(self)
        self.text.setMaximumWidth(400)
        self.clearButton = QtGui.QPushButton('clear', self)
        QtCore.QObject.connect(self.clearButton, QtCore.SIGNAL('clicked()'), self.text.clear)


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
        cbbox.addWidget(self.proc1combo)
        cbbox.addWidget(self.proc2combo)
        cbbox.addWidget(self.primcombo)
        cbbox.addWidget(self.subcombo)
        cbbox.addWidget(self.valcombo)

        bbox = QtGui.QVBoxLayout()
        bbox.addLayout(cbbox)
        bbox.addWidget(self.text)
        bbox.addLayout(tbbox)

        abox = QtGui.QHBoxLayout()
        abox.addLayout(bbox)
        abox.addLayout(pbox)

        self.setLayout(abox)

    data_dic = {}

    def PlotFunc(self):

        global data_dic

        subunit = ['y_pos_16S', 'y_pos_23S']
        subunitc = ['colour_16S', 'colour_23S']
        nucl_data = ['nucl_data_16S', 'nucl_data_23S']
        subunit_neg = ['y_neg_16S', 'y_neg_23S']
        plotname = ['16S Scatterplot', '23S Scatterplot']

        name_conv = {'MazF': 'MazF2h',
                    'MqsR': 'MqsR2h',
                    'Log': 'MG1655log',
                    'Stat': 'MG1655stats',
                    'Delta3': 'delta3'}

        proc1 = self.proc1combo.currentText()
        proc2 = self.proc2combo.currentText()
        prim_select = self.primcombo.currentText()
        sub_select = self.subcombo.currentIndex()
        value_select = self.valcombo.currentIndex()
        proc1 = str(proc1).split(' ')
        proc2 = str(proc2).split(' ')
        prim_select = str(prim_select)

        data_dic = hdf5_gen.get_hdf_data([prim_select, prim_select], [name_conv[proc1[0]], name_conv[proc2[0]]],[proc1[1], proc2[1]],['_input_', 'some'])
        colors = []
        for c in data_dic["data1"][subunitc[sub_select]]:
            if c == 'b':
                c = 'black'
            colors.append(c)

        MA_X_16S = [(math.log(float(y1), 2) + math.log(float(y2), 2))/2 for y1, y2 in zip(data_dic['data1'][subunit[sub_select]], data_dic['data2'][subunit[sub_select]])]
        MA_Y_16S = [math.log((float(y1)/y2), 2) for y1, y2 in zip(data_dic['data1'][subunit[sub_select]], data_dic['data2'][subunit[sub_select]])]
        self.widget1.canvas.ax.clear()
        self.widget1.canvas.ax.set_ylabel('M', fontsize=20)
        self.widget1.canvas.ax.set_xlabel('A', fontsize=20)
        series1 = self.widget1.canvas.ax.scatter(MA_X_16S, MA_Y_16S, alpha=0.5, c=colors, linewidths=( 0, 0, 0), picker=True, label='Datapoints')
        mazf = self.widget1.canvas.ax.scatter(0,0, alpha=0.5, c='red', marker = 'o', label = 'MazF')
        mqsr = self.widget1.canvas.ax.scatter(0,0, alpha=0.5, c='magenta', marker = 'o', label = 'MqsR')
        self.widget1.canvas.ax.legend(handles=[series1,mazf,mqsr],loc='best', scatterpoints = 1)
        self.widget1.canvas.draw()

        if value_select == 0:
            self.widget2.canvas.ax.clear()
            hundred_1 = 0
            hundred_2 = 0
            for i in range(-19,22):
                index = data_dic['data1'][nucl_data[sub_select]].index(float(i))
                hundred_1 += data_dic['data1'][subunit[sub_select]][index]
                hundred_2 += data_dic['data2'][subunit[sub_select]][index]
            heights_1 = [(x / hundred_1 * 100) for x in data_dic['data1'][subunit[sub_select]]]
            heights_2 = [(x / hundred_2 * 100) for x in data_dic['data2'][subunit[sub_select]]]
            self.widget2.canvas.ax.set_ylabel('Relative percentage of reads', fontsize=14)
            self.widget2.canvas.ax.set_xlabel('Nucleotide Position', fontsize=14)
            self.widget2.canvas.ax.set_title(plotname[sub_select], fontsize=14)
            series1 = self.widget2.canvas.ax.scatter(data_dic['data1'][nucl_data[sub_select]], heights_1, alpha=0.5, c=data_dic['data1'][subunitc[sub_select]], linewidths=( 0, 0, 0), picker=True, marker = data_dic['data1']['symbol'], label=' '.join(proc1))
            series2 = self.widget2.canvas.ax.scatter(data_dic['data1'][nucl_data[sub_select]], heights_2, alpha=0.5, c=data_dic['data2'][subunitc[sub_select]], linewidths=( 0, 0, 0), picker=True, marker = data_dic['data2']['symbol'], label=' '.join(proc2))
            self.widget2.canvas.ax.set_ylim(-20,max(heights_1 + heights_2) + 10)
            mazf = self.widget2.canvas.ax.scatter(0,0, alpha=0.5, c='red', marker = 'o', label = 'MazF')
            mqsr = self.widget2.canvas.ax.scatter(0,0, alpha=0.5, c='magenta', marker = 'o', label = 'MqsR')
            self.widget2.canvas.ax.legend(handles=[series1, series2, mazf,mqsr],loc='best', scatterpoints = 1)
            self.widget2.canvas.draw()

        elif value_select == 1:
            self.widget2.canvas.ax.clear()
            self.widget2.canvas.ax.set_ylabel('Read Counts', fontsize=14)
            self.widget2.canvas.ax.set_xlabel('Nucleotide Position', fontsize=14)
            self.widget2.canvas.ax.set_title(plotname[sub_select], fontsize=14)
            series1 = self.widget2.canvas.ax.scatter(data_dic['data1'][nucl_data[sub_select]], data_dic['data1'][subunit[sub_select]], alpha=0.5, c=data_dic['data1'][subunitc[sub_select]], linewidths=( 0, 0, 0), picker=True, marker = data_dic['data1']['symbol'], label=' '.join(proc1))
            self.widget2.canvas.ax.scatter(data_dic['data1'][nucl_data[sub_select]], [-1 * data for data in data_dic['data1'][subunit_neg[sub_select ]]], alpha=0.5, c=data_dic['data1'][subunitc[sub_select]], linewidths=( 0, 0, 0), picker=True, marker = data_dic['data1']['symbol'])
            series2 = self.widget2.canvas.ax.scatter(data_dic['data2'][nucl_data[sub_select]], data_dic['data2'][subunit[sub_select]], alpha=0.5, c=data_dic['data2'][subunitc[sub_select]], linewidths=( 0, 0, 0), picker=True, marker = data_dic['data2']['symbol'], label=' '.join(proc2))
            self.widget2.canvas.ax.scatter(data_dic['data2'][nucl_data[sub_select]], [-1 * data for data in data_dic['data2'][subunit_neg[sub_select ]]], alpha=0.5, c=data_dic['data2'][subunitc[sub_select]], linewidths=( 0, 0, 0), picker=True, marker = data_dic['data2']['symbol'])
            self.widget2.canvas.fig.tight_layout()
            mazf = self.widget2.canvas.ax.scatter(0,0, alpha=0.5, c='red', marker = 'o', label = 'MazF')
            mqsr = self.widget2.canvas.ax.scatter(0,0, alpha=0.5, c='magenta', marker = 'o', label = 'MqsR')
            self.widget2.canvas.ax.legend(handles=[series1, series2, mazf, mqsr],loc='best', scatterpoints = 1)
            self.widget2.canvas.draw()

    def on_combo_prime_change(self, index):
        self.valcombo.setEnabled(False)
        self.valcombo.setCurrentIndex(1)
        text = self.primcombo.currentText()
        if text == '5prime':
            self.proc2combo.clear()
            self.proc2combo.addItems(items5prime.split(','))
            self.proc1combo.clear()
            self.proc1combo.addItems(items5prime.split(','))
        else:
            self.proc2combo.clear()
            self.proc2combo.addItems(items3prime.split(','))
            self.proc1combo.clear()
            self.proc1combo.addItems(items3prime.split(','))



    def onpick(self, event):

        global data_dic

        name_conv = {'MazF': 'MazF2h',
            'MqsR': 'MqsR2h',
            'Log': 'MG1655log',
            'Stat': 'MG1655stats',
            'Delta3': 'delta3'}

        proc1 = self.proc1combo.currentText()
        proc2 = self.proc2combo.currentText()
        prim_select = self.primcombo.currentText()
        sub_select = self.subcombo.currentIndex()
        value_select = self.valcombo.currentIndex()
        proc1 = str(proc1).split(' ')
        proc2 = str(proc2).split(' ')
        prim_select = str(prim_select)
        #data_dic = hdf5_gen.get_hdf_data([prim_select, prim_select], [name_conv[proc1[0]], name_conv[proc2[0]]],[proc1[1], proc2[1]],['_input_', 'some'])

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
                position = int(nucleotide_pos[array_ind]) + 20
                sequence = str(fasta[position-3:position]+ '_' + '<b>' + fasta[position] + '</b>' + fasta[position+1:position+4])
                message = "Primary sequence:<br>{7}<br>Nucleotide position: {0:,} <br>{5} '+' strand: {1:,}<br>{6} '+' strand {2:,}\
                                <br>{5} '-' strand: {3:,}<br>{6} '-' strand {4:,}<br>=======================<br>"\
                                .format(int(nucleotide_pos[array_ind]),
                                        sample_1_pos_read_count[array_ind],
                                        sample_2_pos_read_count[array_ind],
                                        sample_1_neg_read_count[array_ind],
                                        sample_2_neg_read_count[array_ind],
                                        ' '.join(proc1),
                                        ' '.join(proc2),
                                        sequence)

                self.text.insertHtml(message)
                self.text.moveCursor(QtGui.QTextCursor.Start)

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
                position = int(nucleotide_pos[array_ind]) + 20
                sequence = str(fasta[position-3:position]+ '_' + '<b>' + fasta[position] + '</b>' + fasta[position+1:position+4])
                message = "Primary sequence:<br>{7}<br>Nucleotide position: {0:,}<br>Percentages:<br>{3}: {1}<br>{4}:\
                            {2}<br>Absolute values:<br>{3}: {5:,}<br>{4}: {6:,}<br>=======================<br>"\
                            .format(int(nucleotide_pos[array_ind]),
                                    round(sample_1_pos_read_count[array_ind],4),
                                    round(sample_2_pos_read_count[array_ind], 4),
                                    ' '.join(proc1), ' '.join(proc2),
                                    int(sample_1_pos_abs_count[array_ind]),
                                    int(sample_2_pos_abs_count[array_ind]),
                                    sequence)

                self.text.insertHtml(message)
                self.text.moveCursor(QtGui.QTextCursor.Start)

    def save(self):

        proc1 = self.proc1combo.currentText()
        proc2 = self.proc2combo.currentText()
        prim_select = self.primcombo.currentText()
        sub_select = self.subcombo.currentIndex()
        value_select = self.valcombo.currentIndex()
        proc1 = str(proc1).split(' ')
        proc2 = str(proc2).split(' ')
        prim_select = str(prim_select)

        subunit = ['y_pos_16S', 'y_pos_23S']
        nucl_data = ['nucl_data_16S', 'nucl_data_23S']
        subunit_neg = ['y_neg_16S', 'y_neg_23S']

        global data_dic

        path, _ = QtGui.QFileDialog.getSaveFileName(self, 'Save data', os.getcwd(), selectedFilter='*.csv')

        index = ['data1', 'data2']
        with open(path, 'w') as outfile:
            write = csv.writer(outfile, dialect='excel')
            col1 = data_dic['data1'][nucl_data[sub_select]]
            col2 = data_dic['data1'][subunit[sub_select]]
            col3 = data_dic['data1'][subunit_neg[sub_select]]
            col4 = data_dic['data2'][nucl_data[sub_select]]
            col5 = data_dic['data2'][subunit[sub_select]]
            col6 = data_dic['data2'][subunit_neg[sub_select]]
            data = zip(col1, col2, col3, col4, col5, col6)
            for line in data:
                write.writerow(line)




class Compare(QtGui.QWidget):
    def setupUi(self):

        self.plotButton = QtGui.QPushButton('Plot', self)
        QtCore.QObject.connect(self.plotButton, QtCore.SIGNAL('clicked()'), self.plotfunc)

        self.widget = scatter(self)
        self.toolbar = NavigationToolbar(self.widget.canvas, self)
        self.widget.canvas.mpl_connect('pick_event', lambda event: self.onpick(event))

        self.f_prime_combo = QtGui.QComboBox(self)
        self.f_prime_combo.addItems(items5prime.split(','))

        self.t_prime_combo = QtGui.QComboBox(self)
        self.t_prime_combo.addItems(items3prime.split(','))

        self.subcombo = QtGui.QComboBox(self)
        self.subcombo.addItems('16S,23S'.split(','))

        self.valcombo = QtGui.QComboBox(self)
        self.valcombo.addItems('Percentages,Absolute values'.split(','))

        self.text = QtGui.QTextEdit(self)
        self.text.setMaximumWidth(400)
        self.cursor = self.text.textCursor()
        self.clearButton = QtGui.QPushButton('clear', self)
        QtCore.QObject.connect(self.clearButton, QtCore.SIGNAL('clicked()'), self.text.clear)

        self.valcombo.setEnabled(False)
        self.valcombo.setCurrentIndex(0)

        pbox = QtGui.QVBoxLayout()
        pbox.setSpacing(1)
        pbox.addWidget(self.widget)
        pbox.addWidget(self.toolbar)

        tbbox = QtGui.QHBoxLayout()
        tbbox.addWidget(self.plotButton)
        tbbox.addWidget(self.clearButton)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.f_prime_combo)
        vbox.addWidget(self.t_prime_combo)
        vbox.addWidget(self.subcombo)
        vbox.addWidget(self.valcombo)
        vbox.addWidget(self.text)
        vbox.addLayout(tbbox)

        abox = QtGui.QHBoxLayout()
        abox.addLayout(vbox)
        abox.addLayout(pbox)

        self.setLayout(abox)


    data_dic = {}

    def plotfunc(self):

        global data_dic

        name_conv = {'MazF': 'MazF2h',
            'MqsR': 'MqsR2h',
            'Log': 'MG1655log',
            'Stat': 'MG1655stats',
            'Delta3': 'delta3'}


        proc1 = self.f_prime_combo.currentText()
        proc2 = self.t_prime_combo.currentText()
        sub_select = self.subcombo.currentIndex()
        value_select = self.valcombo.currentIndex()
        proc1 = str(proc1).split(' ')
        proc2 = str(proc2).split(' ')

        data_dic = hdf5_gen.get_hdf_data(['5prime', '3prime'], [name_conv[proc1[0]], name_conv[proc2[0]]],[proc1[1], proc2[1]],['_input_', 'some'])

        subunit = ['y_pos_16S', 'y_pos_23S']
        subunitc = ['colour_16S', 'colour_23S']
        nucl_data = ['nucl_data_16S', 'nucl_data_23S']
        subunit_neg = ['y_neg_16S', 'y_neg_23S']
        plotname = ['16S Scatterplot', '23S Scatterplot']


        if value_select == 0:
            self.widget.canvas.ax.clear()
            hundred_1 = sum(data_dic['data1'][subunit[sub_select]])
            hundred_2 = sum(data_dic['data2'][subunit[sub_select]])
            heights_1 = [(x / hundred_1 * 100) for x in data_dic['data1'][subunit[sub_select]]]
            heights_2 = [(x / hundred_2 * 100) for x in data_dic['data2'][subunit[sub_select]]]
            self.widget.canvas.ax.set_ylabel('Relative percentage of reads', fontsize=14)
            self.widget.canvas.ax.set_xlabel('Nucleotide Position', fontsize=14)
            self.widget.canvas.ax.set_title(plotname[sub_select], fontsize=14)
            series1 = self.widget.canvas.ax.scatter(data_dic['data1'][nucl_data[sub_select]], heights_1, alpha=0.5, c=data_dic['data1'][subunitc[sub_select]], linewidths=( 0, 0, 0), picker=True, marker = data_dic['data1']['symbol'], label='5prime ' + ' '.join(proc1))
            series2 = self.widget.canvas.ax.scatter(data_dic['data1'][nucl_data[sub_select]], heights_2, alpha=0.5, c=data_dic['data2'][subunitc[sub_select]], linewidths=( 0, 0, 0), picker=True, marker = data_dic['data2']['symbol'], label='3prime ' + ' '.join(proc2))
            self.widget.canvas.ax.set_ylim(-20,max(heights_1 + heights_2) + 10)
            mazf = self.widget.canvas.ax.scatter(0,0, alpha=0.5, c='red', marker = 'o', label = 'MazF')
            mqsr = self.widget.canvas.ax.scatter(0,0, alpha=0.5, c='magenta', marker = 'o', label = 'MqsR')
            self.widget.canvas.ax.legend(handles=[series1, series2, mazf,mqsr],loc='best', scatterpoints = 1)
            self.widget.canvas.draw()

        elif value_select == 1:
            self.widget.canvas.ax.clear()
            self.widget.canvas.ax.set_ylabel('Read Counts', fontsize=14)
            self.widget.canvas.ax.set_xlabel('Nucleotide Position', fontsize=14)
            self.widget.canvas.ax.set_title(plotname[sub_select], fontsize=14)
            for data_key, data_nt in data_dic.items():
                self.widget.canvas.ax.scatter(data_nt[nucl_data[sub_select]], data_nt[subunit[sub_select]], alpha=0.5, c=data_nt[subunitc[sub_select]], linewidths=( 0, 0, 0), picker=True, marker = data_nt['symbol'])
                self.widget.canvas.ax.scatter(data_nt[nucl_data[sub_select]], [-1 * data for data in data_nt[subunit_neg[sub_select ]]], alpha=0.5, c=data_nt[subunitc[sub_select]], linewidths=( 0, 0, 0), picker=True, marker = data_nt['symbol'])
            self.widget.canvas.fig.tight_layout()
            self.widget.canvas.draw()


    def onpick(self, event):

        global data_dic

        name_conv = {'MazF': 'MazF2h',
            'MqsR': 'MqsR2h',
            'Log': 'MG1655log',
            'Stat': 'MG1655stats',
            'Delta3': 'delta3'}

        proc1 = self.f_prime_combo.currentText()
        proc2 = self.t_prime_combo.currentText()
        sub_select = self.subcombo.currentIndex()
        value_select = self.valcombo.currentIndex()
        proc1 = str(proc1).split(' ')
        proc2 = str(proc2).split(' ')

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
                position = int(nucleotide_pos[array_ind]) + 20
                sequence = str(fasta[position-3:position]+ '_' + '<b>' + fasta[position] + '</b>' + fasta[position+1:position+4])
                message = "Primary sequence:<br>{7}<br>Nucleotide position: {0:,} <br>{5} '+' strand: {1:,}<br>{6} '+' strand {2:,}\
                                <br>{5} '-' strand: {3:,}<br>{6} '-' strand {4:,}<br>=======================<br>"\
                                .format(position,
                                        sample_1_pos_read_count[array_ind],
                                        sample_2_pos_read_count[array_ind],
                                        sample_1_neg_read_count[array_ind],
                                        sample_2_neg_read_count[array_ind],
                                        ' '.join(proc1),
                                        ' '.join(proc2),
                                        sequence)

                self.text.insertHtml(message)
                self.text.moveCursor(QtGui.QTextCursor.Start)

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
                position = int(nucleotide_pos[array_ind]) + 20
                sequence = str(fasta[position-3:position]+ '_' + '<b>' + fasta[position] + '</b>' + fasta[position+1:position+4])
                message = "Primary sequence:<br>{7}<br>Nucleotide position: {0:,}<br>Percentages:<br>{3}: {1}<br>{4}:\
                            {2}<br>Absolute values:<br>{3}: {5:,}<br>{4}: {6:,}<br>=======================<br>"\
                            .format(int(nucleotide_pos[array_ind]),
                                    round(sample_1_pos_read_count[array_ind],4),
                                    round(sample_2_pos_read_count[array_ind],4),
                                    ' '.join(proc1), ' '.join(proc2),
                                    int(sample_1_pos_abs_count[array_ind]),
                                    int(sample_2_pos_abs_count[array_ind]),
                                    sequence)

                self.text.insertHtml(message)
                self.text.moveCursor(QtGui.QTextCursor.Start)













