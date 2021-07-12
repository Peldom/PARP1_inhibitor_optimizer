# -*- coding: utf-8 -*-


'''
author: bin

训练过程父窗体Ui设计类

'''

import sys
from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_TrainBase(QtWidgets.QDockWidget):
    def __init__(self, Form, name):
        super().__init__(Form)
        self.setFeatures(QtWidgets.QDockWidget.AllDockWidgetFeatures)
        self.setWindowTitle(name)


        verticalLayoutWidget = QtWidgets.QWidget()
        verticalLayout = QtWidgets.QVBoxLayout(verticalLayoutWidget)
        verticalLayout.setContentsMargins(10, 10, 10, 0)


        self.stackedWidget = QtWidgets.QStackedWidget(self)
        verticalLayout.addWidget(self.stackedWidget)


        # verticalLayout.addStretch(1)


        horizontalLayoutWidget = QtWidgets.QWidget()
        horizontalLayout = QtWidgets.QHBoxLayout(horizontalLayoutWidget)

        horizontalLayout.addStretch(1)

        horizontalLayout.addWidget(
            QtWidgets.QLabel("当前显示:", horizontalLayoutWidget))

        self.spinBox = QtWidgets.QSpinBox(horizontalLayoutWidget)
        self.spinBox.setMinimum(1)
        self.spinBox.valueChanged[int].connect(lambda i:
            self.stackedWidget.setCurrentIndex(i - 1))
        horizontalLayout.addWidget(self.spinBox)

        horizontalLayout.addStretch(1)

        nextBtn = QtWidgets.QPushButton("正在运行页面", horizontalLayoutWidget)
        nextBtn.clicked.connect(
            lambda : self.spinBox.setValue(self.stackedWidget.count()))
        horizontalLayout.addWidget(nextBtn)

        horizontalLayout.addStretch(1)

        verticalLayout.addWidget(horizontalLayoutWidget)


        self.setWidget(verticalLayoutWidget)

        self.close()

        QtCore.QMetaObject.connectSlotsByName(self)



    @QtCore.pyqtSlot(object)
    def addWidget(self, widget):
        self.show()
        self.stackedWidget.addWidget(widget)
        self.stackedWidget.setCurrentWidget(widget)
        self.spinBox.setMaximum(self.stackedWidget.count())
        self.spinBox.setValue(self.stackedWidget.currentIndex() + 1)

