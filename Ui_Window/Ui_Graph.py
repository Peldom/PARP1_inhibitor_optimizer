# -*- coding: utf-8 -*-


'''
author: bin

实时数据可视化框类

'''

import sys
import psutil

from PyQt5 import QtCore, QtGui, QtWidgets
import pyqtgraph as pg


class Ui_Graph(pg.PlotWidget, QtWidgets.QWidget):
    signal_data = QtCore.pyqtSignal(dict)

    def __init__(self, title, data, pad=10, background='w'):
        """ 实时数据可视化框

            tltle:图表名

            data:数据型 dict{str:[]}

            pad:信号发出间隔

            background:背景色
        """

        super().__init__(background=background, title=title)
        QtWidgets.QWidget().__init__()

        self.data = data
        self.pad = pad
        self.setMinimumSize(400, 260)
        self.addLegend(offset=(-10, 1))

        i, num = 0, len(self.data)
        for key in self.data:
            self.plot(pen=(i, num), name=key)
            i += 1


    def update(self, tmp):
        """ 实时绘图

            tmp:字典，键需与构造时相同,值为数值型

            signal_data[dict]:每pad次数据更新，发信号返回数据
        """
        i, num = 0, len(self.data)
        for key in self.data:
            if type(tmp[key]) == list:
                self.data[key] += tmp[key]
            else:
                if i == 1 and len(self.data[key]) % self.pad == 0:
                    self.signal_data.emit(tmp)
                self.data[key].append(tmp[key])

            self.plot(self.data[key], pen=(i, num))
            i += 1


    @QtCore.pyqtSlot()
    def lastData(self):
        """ 信号发出并返回最后更新的数据

            return: dict
        """
        t = {}
        for key, value in self.data.items():
            t[key] = value[-1]
        self.signal_data.emit(t)
        return t





    

   
