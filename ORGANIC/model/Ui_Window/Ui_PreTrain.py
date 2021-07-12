# -*- coding: utf-8 -*-


'''
author: bin

预训练子窗体类

'''

import sys

from PyQt5 import QtCore, QtGui, QtWidgets

from Ui_Window.Ui_Graph import Ui_Graph

class Ui_PreTrain(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()

        self.label, self.graphWidget, self.progressBar = {}, {}, {}


        verticalLayout = QtWidgets.QVBoxLayout(self)

        # label = QtWidgets.QLabel("生成器预训练", self)
        # label.setAlignment(QtCore.Qt.AlignCenter)
        # verticalLayout.addWidget(label)

        self.graphWidget["pre_gen"] = Ui_Graph("生成器预训练", {
            "pre_gen_loss":[]}) 
        verticalLayout.addWidget(self.graphWidget["pre_gen"])

        self.progressBar["pre_gen"] = QtWidgets.QProgressBar(self)
        self.progressBar["pre_gen"].setMinimum(0)
        self.progressBar["pre_gen"].setFormat("%v / %m")
        verticalLayout.addWidget(self.progressBar["pre_gen"])
 
        self.label["pre_gen"] = QtWidgets.QLabel(self)
        self.graphWidget["pre_gen"].signal_data[dict].connect(
            lambda data : self.label["pre_gen"].setText(
                "生成器预训练 loss: {:10f} ".format(data["pre_gen_loss"])))
        verticalLayout.addWidget(self.label["pre_gen"])


        verticalLayout.addStretch(1)


        # label = QtWidgets.QLabel("判别器预训练", self)
        # label.setAlignment(QtCore.Qt.AlignCenter)
        # verticalLayout.addWidget(label)

        self.graphWidget["pre_dis"] = Ui_Graph("判别器预训练", {
                "pre_dis_loss":[],
                "pre_dis_accuracy":[]}) 
        verticalLayout.addWidget(self.graphWidget["pre_dis"])

        self.progressBar["pre_dis"] = QtWidgets.QProgressBar(self)
        self.progressBar["pre_dis"].setMinimum(0)
        self.progressBar["pre_dis"].setFormat("%v / %m")
        verticalLayout.addWidget(self.progressBar["pre_dis"])

        self.label["pre_dis"] = QtWidgets.QLabel(self)
        self.graphWidget["pre_dis"].signal_data[dict].connect(
            lambda data : self.label["pre_dis"].setText(
                "生成器预训练 loss: {:10f} accuracy: {:10f}".format(
                    data["pre_dis_loss"],data["pre_dis_accuracy"])))
        verticalLayout.addWidget(self.label["pre_dis"])


        verticalLayout.addStretch(1)


        self.label["pre"] = QtWidgets.QLabel(self)
        verticalLayout.addWidget(self.label["pre"])


        verticalLayout.addStretch(1)



if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    w = Ui_PreTrain()
    w.show()
    sys.exit(app.exec_())
