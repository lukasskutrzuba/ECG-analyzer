#include "ecganalyzer.h"
#include <QtGui/QApplication>

int main(int argc, char *argv[])
{
	//ECGController kontroler;
	//kontroler.runECGBaseline();
	//kontroler.runRPeaks();
    //kontroler.runHRV2();
	QApplication a(argc, argv);
	ECGanalyzer w;
	w.show();
	return a.exec();
	//return 0;
}
