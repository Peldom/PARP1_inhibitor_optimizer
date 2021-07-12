# -*- coding: utf-8 -*-


'''
author: bin

页面基础逻辑类

'''

import json
import time
import sys
import os
import re

from PyQt5 import QtCore, QtGui, QtWidgets

from Ui_Window.Ui_Main import Ui_Main


from Ui_Window.Td_Organic import Td_Organic


class mainWindow(Ui_Main):
    def __init__(self):
        super().__init__()

        self.threads = {}



    @QtCore.pyqtSlot()
    def on_settingTool_clicked(self):
        self.trainSetting.show()



    @QtCore.pyqtSlot()
    def on_startTool_clicked(self):
        self.trainSetting.close()
        self.logo.close()
        self.trainSchedule.show()

        setting = self.trainSetting.getSetting()
        self.threads["worker"] = Td_Organic("worker", setting)
        self.threads["worker"].signal_progress[str, int, int].connect(self.trainSchedule.progress)
        self.threads["worker"].signal_inf[list].connect(self.trainSchedule.loadInf)
        self.threads["worker"].signal_draw[str, dict].connect(self.trainSchedule.graphUpdate)
        self.threads["worker"].signal_message[str, str].connect(self.trainSchedule.message)
        self.threads["worker"].signal_finished[str, str].connect(self.threadFinished)
        self.threads["worker"].start()



    @QtCore.pyqtSlot(str, str)
    def threadFinished(self, name, str):
        if not name:
            self.threads[name].quit()
        _ = QtWidgets.QMessageBox.information(self, '提示', str)




if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    w = mainWindow()
    sys.exit(app.exec_())