# primer-design

First, make sure you have the most recent version of Python installed in your computer.

Windows:
	Visit https://www.python.org/downloads/windows/ and download the installer.
	Follow instructions from the GUI

Mac OS:
	Before installing Python, you’ll need to install GCC. GCC can be obtained by downloading XCode or the smaller Command Line Tools (must have an Apple account).
	Xcode 10 can be found at https://developer.apple.com/xcode/
	Command Line Tools can be found at https://developer.apple.com/downloads/
	Open Terminal and type the following command:
	
		$ ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
		$ brew install python
	
Linux:
	To verify what version of Python 3 is installed, open the Terminal and type:
		
		$ python3 --version
	If you are using Ubuntu 16.10 or newer, then you can easily install Python 3.7 with the following commands:
		
		$ sudo apt-get update
		$ sudo apt-get install python3.7

	
Now, to run the script, open the terminal (mac / linux, command prompt on windows) and enter the following:
	
		$ python3 prmrdsgn.py [enter the arguments you need, from the list below]

			LIST OF ARGUMENTS
			-i [required]		'path to input file in FASTA format'
			-n or --number		'number of candidate primer pairs to pick, default = 5'
			-e or --extension	'number of bases from the start and end of the sequence to look for primers in, default = 100'
			-s or --short		'shortest acceptable primer, default = 20'
			-l or --long		'longest acceptable primer, default = 30'
			-m or --mintemp		'min Tm in celsius, default = 55'
			-x or --maxtemp		'max Tm in celsius, default = 75'
			-M or --mingc		'min GC percentage, default = 40'
			-X or --maxgc		'max GC percentage, default = 65'
			-D or --tmdiff		'accepted TM difference to form primer pair, default = 0.5'
			-o or --output		'path to output file in FASTA format, default = output'
			-v or --verbose		'prints each step of each iteration (for debugging)'







