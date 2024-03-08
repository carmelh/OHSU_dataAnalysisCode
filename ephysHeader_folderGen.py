# -*- coding: utf-8 -*-
"""
Created on Wed Feb 8 2023

Script to run through all the coverslip/cell folders to make the standard folders
for currentClamp and voltageClamp

@author: Carmel Howe
"""

import sys
sys.path.insert(1, r'Z:\Labs\Frank Lab\SHARED\000Frank lab shared\Data Analysis Scripts\Python\ephys_functions_dont_change')
import ephys_analysis_funcs_dontChange as ef
from PyQt5 import QtCore, QtGui, uic, QtWidgets


qtCreatorFile = r'Z:\Labs\Frank Lab\SHARED\000Frank lab shared\Data Analysis Scripts\Python\ephys_functions_dont_change\folderGenerator.ui'

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
 
class MyApp(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.browseButton.clicked.connect(self.getFolderPath)
        self.makeFoldersButton.clicked.connect(self.generateFolders)
        self.closeButton.clicked.connect(self.closeEvent)

#  user selects parent folder, i.e. the experiment day
    def getFolderPath(self):
        folder_path = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select Folder')
        self.pathTextBox.setText(folder_path)
        return folder_path

# generates the folders within the parent folder selected by user
    def generateFolders(self): 
        folder_path = self.pathTextBox.toPlainText()
        if self.expComboBox.currentText() == "Current Clamp":
            self.printNoteBox.setText("Current clamp")
            expParam = "cClamp"
            ef.makeFoldersAnalysis(folder_path,expParam)
        elif self.expComboBox.currentText() == "Voltage Clamp":
            self.printNoteBox.setText("Voltage clamp")
            expParam = "vClamp"
            ef.makeFoldersAnalysis(folder_path,expParam)
        else: 
            self.printNoteBox.setText("Select an experiment")   
        
        self.printNoteBox_2.setText("Go organise trials into folders then open analysis GUI")

    def closeEvent(self):
        self.close()
        
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MyApp()
    window.show()
    sys.exit(app.exec_())