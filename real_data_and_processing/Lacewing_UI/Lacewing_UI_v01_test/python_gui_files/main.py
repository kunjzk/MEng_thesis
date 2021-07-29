# File:		main.py
# Author: 	Lei Kuang
# Date:		2nd May 2020
# @ Imperial College London

from PyQt5.QtWidgets import QApplication
from Lacewing import Lacewing_Dialog

if __name__ == "__main__":
    import sys

    app = QApplication(sys.argv)
    ui = Lacewing_Dialog()
    ui.show()
    
    sys.exit(app.exec_())
