/** @file */
#include "displayFunctions.h"

/// Prints progress bar for given percent completion. To be called in a loop for percent 0 to 100.
void printProgressBar(double percent){
	float progress = 0.0;
	while (progress < 1.0) {
		int barWidth = 70;

	        cout << "[";
	        int pos = percent/100.*barWidth;
	        for (Index i = 0; i < barWidth; ++i) {
		        if (i < pos) cout << "=";
		        else if (i == pos) cout << ">";
		        else cout << " ";
	        };
	        cout << "] " << int(percent) << " %\r";
	        cout.flush();
		progress += 0.16;
	}
};

/// Prints progress bar for given percent completion. To be called in a loop for percent 0 to 100.
/// Prints given 'message' at the end of the progress bar (such as the Courant number or max RES etc.)
void printProgressBar(double percent, string message){
	float progress = 0.0;
	while (progress < 1.0) {
		int barWidth = 70;

	        cout << "[";
	        int pos = percent/100.*barWidth;
	        for (Index i = 0; i < barWidth; ++i) {
		        if (i < pos) cout << "=";
		        else if (i == pos) cout << ">";
		        else cout << " ";
	        };
	        cout << "] " << int(percent) << " %" << "   " << message << "\r";
	        cout.flush();
		progress += 0.16;
	}
};

/// Prints a horizontal line on the screen consisting of the arg string. 
/// example: printLine("_") will have building blocks of the line as "_" characters.
void printLine(string a){
	for (Index i=0; i<120; i++){
		cout << a;
	};
	cout << endl;

};
