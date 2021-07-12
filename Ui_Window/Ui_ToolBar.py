'''
author: bin

训练参数子窗体Ui设计类

'''

import sys
from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_ToolBar(QtWidgets.QDockWidget):

    def __init__(self, Form):
        super().__init__(Form)
        self.setFeatures(QtWidgets.QDockWidget.DockWidgetMovable|QtWidgets.QDockWidget.DockWidgetFloatable)
        self.setWindowTitle("工具栏")

        horizontalLayoutWidget = QtWidgets.QWidget()
        horizontalLayoutWidget.setMaximumHeight(85)
        horizontalLayout = QtWidgets.QHBoxLayout(horizontalLayoutWidget)

        def toolbarAddAction(name, icon, tip):
            button = QtWidgets.QPushButton(QtGui.QIcon(icon), "", horizontalLayoutWidget)
            button.setObjectName(name)
            button.setIconSize(QtCore.QSize(60, 60))
            button.setStatusTip(tip)
            horizontalLayout.addWidget(button)
        #     horizontalLayout.addStretch(0)

        # horizontalLayout.addStretch(0)
        toolbarAddAction("settingTool", "./img/setting.png", "设置模型训练参数")
        toolbarAddAction("startTool", "./img/start.png", "开始训练")


        self.setWidget(horizontalLayoutWidget)


        

