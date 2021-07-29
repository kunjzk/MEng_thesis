# Lacewing_UI
Python GUI for Lacewing Platform

# Version
### \# v00
* Initial release
### \# v01
* Add UDP mode for high frame rate (When UDP mode is active, no debug functions will be supported)
* Fix bluetooth bug (if the UI fails to set up the link, following steps in "Quick Start for Bluetooth"

### Quick Start for Bluetooth
For devices you have paired, ignore steps [01]~[04]

* [01] !!! DO NOT OPEN THE UI (KEEP IT CLOSED)
* [02] Power on the Lacewing board
* [03] Pair the board through "WIN -> Settings -> Devices -> Add bluetooths & other devices -> Bluetooth"
* [04] Choose a Lacewing device, e.g., "LACEWING_XXXX", where XXXX is the bluetooth device (PIN: 1234 if applicable)
* [05] !!! POWER OFF THE LACEWING BOARD
* [06] Launch Lacewing_UI
* [07] Wait for scanning, once it finished, the "Serial Status" will become GREEN (Idle)
* [08] Power on the LacewingPXE board you have paired
* [09] Choose a Lacewing deice from the "Serial Port"
* [10] Wait for the confirmation of connection (If timeout, power cycle LacewingPXE and try again)
* [!!!] Golden Rule: Do not click on any button when the "Serial Status" is RED (Busy), if you are not sure



### IDE
> Eric6  
> https://eric-ide.python-projects.org/eric-download.html  
> Python 3.7.7  
> https://www.python.org/downloads/release/python-377/  

### Required Python Module (After you install Eric6)
* pyserial (pip install pyserial==3.3)  
* pybluez
* numpy
* matplotlib
* opencv-python
* pyqt5-tools (You need to copy the bin folder from this package and replace the same folder in PyQt5)  
[1] if you install the Python in folder "C:\Python37"  
[2] copy the folder "C:\Python37\Lib\site-packages\pyqt5_tools\Qt\bin"  
[3] then go to "C:\Python37\Lib\site-packages\PyQt5\Qt" and paste (√ Replace the files in the destination)  
[!!] this step is necessary if and only if you want to open the ui design file (Lacewing.ui) for qt  
