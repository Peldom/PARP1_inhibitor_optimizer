# -*- coding: utf-8 -*-


'''
author: bin

主页面Ui设计类

'''

import json
import time
import sys
import os
import re

from PyQt5 import QtCore, QtGui, QtWidgets

from Ui_Window.Ui_ToolBar import Ui_ToolBar
from Ui_Window.Ui_Logo import Ui_Logo

from Ui_Window.Ui_Setting import Ui_Setting
from Ui_Window.Ui_Schedule import Ui_Schedule


class Ui_Main(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()

        self.setGeometry(100,50,1500,900)
        self.setWindowTitle("ORGANIC")
        self.setDockNestingEnabled(True)
        self.statusBar().showMessage('准备就绪')


        self.toolbar = Ui_ToolBar(self)
        self.logo = Ui_Logo(self)

        
        self.trainSetting = Ui_Setting(self)
        self.trainParams = self.trainSetting.trainParams
        self.trainParticular = self.trainSetting.trainParticular


        self.trainSchedule = Ui_Schedule(self)
        self.trainInf = self.trainSchedule.trainInf
        self.perTrain = self.trainSchedule.perTrain
        self.training = self.trainSchedule.training


        # 此处需严格遵守 从左到右，从上倒下

        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, self.trainParams)
        self.splitDockWidget(self.trainParams, self.trainInf, QtCore.Qt.Horizontal)
        self.splitDockWidget(self.trainInf, self.toolbar, QtCore.Qt.Horizontal)
        self.splitDockWidget(self.toolbar, self.trainSchedule, QtCore.Qt.Horizontal)
        self.splitDockWidget(self.trainSchedule, self.trainSetting, QtCore.Qt.Horizontal)

        self.splitDockWidget(self.toolbar, self.logo, QtCore.Qt.Vertical)
        self.splitDockWidget(self.trainSchedule, self.perTrain, QtCore.Qt.Vertical)
        self.splitDockWidget(self.trainSetting, self.trainParticular, QtCore.Qt.Vertical)
        
        self.splitDockWidget(self.logo, self.training, QtCore.Qt.Vertical)

        QtCore.QMetaObject.connectSlotsByName(self)

        self.show()


        





