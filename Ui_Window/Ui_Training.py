# -*- coding: utf-8 -*-


'''
author: bin

训练子窗体类

'''

import sys
import psutil

from PyQt5 import QtCore, QtGui, QtWidgets

from Ui_Window.Ui_Graph import Ui_Graph

class Ui_Training(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()

        self.label, self.graphWidget, self.progressBar = {}, {}, {}


        gridLayout = QtWidgets.QGridLayout(self)

        self.progressBar["batch"] = QtWidgets.QProgressBar(self)
        self.progressBar["batch"].setFormat("当前批次进度为: %v / %m")
        self.progressBar["batch"].setMinimum(0)
        gridLayout.addWidget(self.progressBar["batch"], 0, 0, 1, 2)


        self.graphWidget["gen"] = Ui_Graph("Rewards",{   
                "gen_mean":[],
                "gen_std":[],
                "gen_min":[],
                "gen_max":[]}) 
        gridLayout.addWidget(self.graphWidget["gen"], 1, 0, 2, 1)

        self.progressBar["gen"] = QtWidgets.QProgressBar(self)
        self.progressBar["gen"].setMinimum(0)
        self.progressBar["gen"].setFormat("%v / %m")
        gridLayout.addWidget(self.progressBar["gen"], 1, 1)
 
        self.label["gen"] = QtWidgets.QLabel(self)
        self.graphWidget["gen"].signal_data[dict].connect(
            lambda data : self.label["gen"].setText(
                "生成器训练\n" + "\n".join([
                    "{}: {:10f}".format(k,v)  for k,v in data.items()])))
        gridLayout.addWidget(self.label["gen"], 2, 1)


        self.graphWidget["dis"] = Ui_Graph(
            "判别器训练", {
                "dis_loss":[],
                "dis_accuracy":[]}) 
        gridLayout.addWidget(self.graphWidget["dis"], 3, 0, 2, 1)

        self.progressBar["dis"] = QtWidgets.QProgressBar(self)
        self.progressBar["dis"].setMinimum(0)
        self.progressBar["dis"].setFormat("%v / %m")
        gridLayout.addWidget(self.progressBar["dis"], 3, 1)

        self.label["dis"] = QtWidgets.QLabel(self)
        self.graphWidget["dis"].signal_data[dict].connect(
            lambda data : self.label["dis"].setText(
                "生成器训练\nloss: {:10f}\naccuracy: {:10f}".format(
                    data["dis_loss"],data["dis_accuracy"])))
        gridLayout.addWidget(self.label["dis"], 4, 1)



        self.label["train"] = QtWidgets.QLabel(self)
        gridLayout.addWidget(self.label["train"], 5, 0, 1, 2)