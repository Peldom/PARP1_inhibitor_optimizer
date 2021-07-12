# -*- coding: utf-8 -*-


'''
author: bin

总体进度显示子窗体类

'''

import sys
from PyQt5 import QtCore, QtGui, QtWidgets

from Ui_Window.Ui_TrainBase import Ui_TrainBase
from Ui_Window.Ui_TrainInformation import Ui_TrainInformation
from Ui_Window.Ui_PreTrain import Ui_PreTrain
from Ui_Window.Ui_Training import Ui_Training



class Ui_Schedule(QtWidgets.QDockWidget):
    def __init__(self, Form):
        super().__init__(Form)
        self.setFeatures(
            QtWidgets.QDockWidget.DockWidgetMovable|
            QtWidgets.QDockWidget.DockWidgetFloatable)
        self.setWindowTitle("总体进度")

        self.trainInf = Ui_TrainBase(self, "模型训练信息")
        self.perTrain = Ui_TrainBase(self, "预训练过程")
        self.training = Ui_TrainBase(self, "训练过程")

        self.progressBar  = {}


        horizontalLayoutWidget = QtWidgets.QWidget()
        horizontalLayout = QtWidgets.QHBoxLayout(horizontalLayoutWidget)

        self.progressBar["global"] = QtWidgets.QProgressBar(horizontalLayoutWidget)
        self.progressBar["global"].setMinimum(0)
        self.progressBar["global"].setFormat("总循环进度为: %v / %m")
        horizontalLayout.addWidget(self.progressBar["global"])

        nextBtn = QtWidgets.QPushButton("正在运行页面", horizontalLayoutWidget)
        nextBtn.clicked.connect(self.currentIndex)
        horizontalLayout.addWidget(nextBtn)

        self.setWidget(horizontalLayoutWidget)
        
        self.close()

        QtCore.QMetaObject.connectSlotsByName(self)



    def show(self):
        super().show()
        self.trainInf.show()
        self.perTrain.show()
        self.training.show()


    def close(self):
        super().close()
        self.trainInf.close()
        self.perTrain.close()
        self.training.close()


    def currentIndex(self):
        value = self.progressBar["global"].maximum()
        self.trainInf.stackedWidget.setCurrentIndex(value)
        self.perTrain.stackedWidget.setCurrentIndex(value)
        self.training.stackedWidget.setCurrentIndex(value)


    @QtCore.pyqtSlot(list)
    def loadInf(self, inf):
        """ 训练相关信息，信息显示页面更新

            inf:信息 list(tuple) [(,)]
        """
        self.trainInf.addWidget(Ui_TrainInformation(inf))



    @QtCore.pyqtSlot(str, int, int)
    def progress(self, name, step, loop):
        """ 进度条分发，页面生成，引用更新

            name:进度条名称(需保证各页面各进度条无重名) str

            step:当前进度 int

            loop:总进度 int
        """
        self.progressBar[name].setValue(step)
        self.progressBar[name].setMaximum(loop)

        if step == loop: 
            if name in self.graphWidget:
                self.graphWidget[name].lastData()
            return
        elif name != "global":
            return 

        # 新建页面
        preTrain = Ui_PreTrain()
        training = Ui_Training()

        # 更新引用至当前页面
        self.progressBar = {
            "global": self.progressBar["global"],
            **preTrain.progressBar,
            **training.progressBar}
        self.graphWidget = {
            **preTrain.graphWidget,
            **training.graphWidget}
        self.lable = {
            "pre": preTrain.label["pre"],
            "train": training.label["train"]}

        # 添加页面
        self.perTrain.addWidget(preTrain)
        self.training.addWidget(training)



    @QtCore.pyqtSlot(str, dict)
    def graphUpdate(self, name, tmp):
        """ 实时绘图分发

            name:绘图窗口名 str

            tmp:实时数据 dict
        """
        self.graphWidget[name].update(tmp)



    @QtCore.pyqtSlot(str, str)
    def message(self, name, message):
        """ 消息显示

            name:显示标签名 str

            message:消息 str
        """
        self.lable["pre"].setText(message)
        
        









