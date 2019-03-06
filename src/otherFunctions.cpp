#include "otherFunctions.h"

using namespace std;


/// Converts string to int
int stringToInt(string input)
{
	stringstream convert(input);
	int output;
	convert >> output;
	return output;
};

/// Converts string to float
float stringToFloat(string input)
{
	stringstream convert(input);
	float output;
	convert >> output;
	return output;
};

/// Converts string to double
double stringToDouble(string input)
{
	stringstream convert(input);
	double output;
	convert >> output;
	return output;
};

/// Returns factorial of first k integers
int fact(int k){
	int p = 1;
	for (Index i=0; i<k; i++){
		p = p * (k - i);
	};
	return p;
};


/// Checks whether the given file exists
bool fileExists(const std::string& filename){
	std::ifstream ifile(filename.c_str());
	return (bool)ifile;
}


void generatePlottingScript(string name, int order){
	// Generates a plotting script for plotting the polynomials.
	// the script is generated inside dir tests/
	// first change directory to tests/ Then just run the script with $python plot_name_polynomials.py

	ofstream file;
	file.open(("tests/plot_"+name+"_polynomials.py").c_str());
		file <<  "import numpy as np" << endl;
		file << "import pylab\n";
		file << "import sys\n";

		//file << "f = sys.argv[1]" << endl;
		file << "f = " << "\"tests/" << name << "_polynomials.dat" << "\"" << endl;

		file << "data = np.loadtxt(f)"<< endl;
		file << "n = np.shape(data)[0]"<< endl;
		file << "x = np.linspace(-1,1,n)"<< endl;
		file << "pylab.figure "<< endl;
		file << "for i in range(np.shape(data)[1]): "<< endl;
		file << "	pylab.plot(x,data[:,i],label=r" << "\"$\\phi_{%d}$\"%i)" << endl;
		file << "pylab.legend(loc=1) "<< endl;
		file << "pylab.title (" << "\"" << name   << " " << "polynomials" << "\")"<< endl;
		file << "pylab.xlabel('x') "<< endl;
		file << "pylab.xlim(-1,1.2) "<< endl;
		file << "pylab.ylabel(r'$\\phi$') "<< endl;
		file << "pylab.grid()"<< endl;
		file << "pylab.savefig('tests/" << name << "_" << varToString(order) <<".eps')"<< endl;
	file.close();
};

void generatePlottingScriptAssignment1(string name){
	// Generates a plotting script for plotting the polynomials.
	// the script is generated inside dir tests/
	// first change directory to tests/ Then just run the script with $python plot_name_polynomials.py

	ofstream file;
	file.open(("tests/plot_error_"+name+".py").c_str());
		file <<  "import numpy as np" << endl;
		file << "import pylab\n";
		file << "import sys\n";

		//file << "f = sys.argv[1]" << endl;
		file << "f = " << "\"tests/Assignment1_interpolation.dat" << "\"" << endl;

		file << "data = np.loadtxt(f)"<< endl;
		//file << "n = np.shape(data)[0]"<< endl;
		file << "x = data[:,0]"<< endl;
		file << "fig1 = pylab.figure() "<< endl;
		file << "ax1 = fig1.add_subplot(111)" << endl;
		file << "for i in range(1,np.shape(data)[1]): "<< endl;
		file << "	if (i != 3):\n";
		file << "		ax1.semilogy(x,data[:,i],label=r" << "\"$L_{%d}$\"%i"<< "+ \" error\")" << endl;
		file << "	else:\n";
		file << "		ax1.semilogy(x,data[:,i],label=r" << "\"$L_{\\infty}$\""<< "+\" error\")" << endl;
		file << "pylab.legend(loc=1) "<< endl;
		file << "pylab.title (" << "\" Interpolation error" << "\")"<< endl;
		file << "pylab.xlabel('N (order of accuracy)') "<< endl;
		//file << "pylab.xlim(-1,1.2) "<< endl;
		file << "pylab.ylabel(r'$||error||$') "<< endl;
		file << "pylab.grid()"<< endl;
		file << "pylab.savefig('tests/interpolation_error.eps')"<< endl;

		file << "fig2 = pylab.figure() "<< endl;
		file << "ax2 = fig2.add_subplot(111)" << endl;
		file << "f = " << "\"tests/Assignment1_differentiation.dat" << "\"" << endl;
		file << "data = np.loadtxt(f)"<< endl;
		//file << "n = np.shape(data)[0]"<< endl;
		file << "x = data[:,0]"<< endl;
		file << "for i in range(1,np.shape(data)[1]): "<< endl;
		file << "	if (i != 3):\n";
		file << "		ax2.semilogy(x,data[:,i],label=r" << "\"$L_{%d}$\"%i"<< "+ \" error\")" << endl;
		file << "	else:\n";
		file << "		ax2.semilogy(x,data[:,i],label=r" << "\"$L_{\\infty}$\""<< "+\" error\")" << endl;
		file << "pylab.legend(loc=1) "<< endl;
		file << "pylab.title (" << "\" Differentiation error" << "\")"<< endl;
		file << "pylab.xlabel('N (order of accuracy)') "<< endl;
		//file << "pylab.xlim(-1,1.2) "<< endl;
		file << "pylab.ylabel(r'$||error||$') "<< endl;
		file << "pylab.grid()"<< endl;
		file << "pylab.savefig('tests/differentiation_error.eps')"<< endl;
	file.close();
};

void generatePlottingScriptAssignment2(string name){
	// Generates a plotting script for plotting the polynomials.
	// the script is generated inside dir tests/
	// first change directory to tests/ Then just run the script with $python plot_name_polynomials.py

	ofstream file;
	file.open(("tests/plot_error_"+name+".py").c_str());
		file <<  "import numpy as np" << endl;
		file << "import pylab\n";
		file << "import sys\n";

		//file << "f = sys.argv[1]" << endl;
		file << "f = " << "\"tests/" << name <<"_interpolation.dat" << "\"" << endl;

		file << "data = np.loadtxt(f)"<< endl;
		//file << "n = np.shape(data)[0]"<< endl;
		file << "x = data[:,0]"<< endl;
		file << "fig1 = pylab.figure() "<< endl;
		file << "ax1 = fig1.add_subplot(111)" << endl;
		file << "for i in range(1,np.shape(data)[1]): "<< endl;
		file << "	if (i != 3):\n";
		file << "		ax1.semilogy(x,data[:,i],label=r" << "\"$L_{%d}$\"%i"<< "+ \" error\")" << endl;
		file << "	else:\n";
		file << "		ax1.semilogy(x,data[:,i],label=r" << "\"$L_{\\infty}$\""<< "+\" error\")" << endl;
		file << "pylab.legend(loc=1) "<< endl;
		file << "pylab.title (" << "\" Interpolation error: 2D domain" << "\")"<< endl;
		file << "pylab.xlabel('N (order of accuracy)') "<< endl;
		//file << "pylab.xlim(-1,1.2) "<< endl;
		file << "pylab.ylabel(r'$||error||$') "<< endl;
		file << "pylab.grid()"<< endl;
		file << "pylab.savefig('tests/interpolation_error_" << name << ".eps')"<< endl;

	file.close();
};



void generatePlottingScript2D(string name, int order){
	// Generates a plotting script for plotting the polynomials.
	// the script is generated inside dir tests/
	// first change directory to tests/ Then just run the script with $python plot_name_polynomials.py
	int Np = order+1;

	ofstream file;
	file.open(("tests/plot_"+name+"_polynomials.py").c_str());
		file <<  "from mpl_toolkits.mplot3d import Axes3D" << endl;
		file << "import matplotlib.pyplot as plt"<<endl;
		file << "import numpy as np"<<endl;
		file << "import sys"<<endl;

		file << "f = " << "\"tests/" << name << "_polynomials.dat" << "\"" << endl;

		file << "data = np.loadtxt(f)"<< endl;
		file << "n = int(np.sqrt(np.shape(data)[0]))"<< endl;
		file << "m = int(np.sqrt(np.shape(data)[1]))"<< endl;
		file << "X = np.linspace(-1,1,n)"<< endl;
		file << "Y = np.linspace(-1,1,n)"<< endl;
		file << "X, Y = np.meshgrid(X, Y)"<< endl;
		file << "fig = plt.figure(figsize=("<<4*Np<<","<<4*Np<<"))"<< endl;
		file << "fig.suptitle('2D "<<name<<" polynomials')"<< endl;
		file << "counter = 0"<< endl;
		file << "for i in xrange(m):"<< endl;
		file << "\tfor j in xrange(m):"<< endl;
		file << "\t\tax = fig.add_subplot(m,m,counter+1, projection='3d')"<< endl;
		file << "\t\tphi = data[:,counter]"<< endl;
		file << "\t\tphi = phi.reshape((n,n))"<< endl;
		file << "\t\tsurf = ax.plot_surface(X, Y, phi, rstride=1, cstride=1,linewidth=0, antialiased=False)"<< endl;
		file << "\t\tax.set_zlim3d(-1.1, 1.1)"<< endl;
		file << "\t\tcounter+=1"<< endl;
		file << "plt.savefig('tests/"<<name<<"_"<<order<<".eps',bbox_inches='tight')"<< endl;
		file << "plt.close()"<< endl;
	file.close();
};

