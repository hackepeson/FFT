#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "kiss_fftr.h"


float dt = 0.001;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_iSpectrumSelected(0)
{
    ui->setupUi(this);

    connect(ui->actionOpen, SIGNAL(triggered()), SLOT(openFile()));
    connect(ui->actionOpen_MS_Acc, SIGNAL(triggered()), SLOT(openMSAcc()));

    // QCustomPlot
    ui->widgetDataPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes |
                                        QCP::iSelectLegend | QCP::iSelectPlottables);

    ui->widgetFFTPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes |
                                       QCP::iSelectLegend | QCP::iSelectPlottables);
    // connect slot that ties some axis selections together (especially opposite axes):
    connect(ui->widgetDataPlot, SIGNAL(selectionChangedByUser()), this, SLOT(selectionChanged()));
    connect(ui->widgetFFTPlot, SIGNAL(selectionChangedByUser()), this, SLOT(selectionChanged()));

    connect(ui->widgetDataPlot, SIGNAL(mouseWheel(QWheelEvent*)),  SLOT(mouseWheel(QWheelEvent*)));
    connect(ui->widgetFFTPlot, SIGNAL(mouseWheel(QWheelEvent*)),  SLOT(mouseWheel(QWheelEvent*)));

    // connect slots that takes care that when an axis is selected, only that direction can be dragged and zoomed:
    connect(ui->widgetDataPlot, SIGNAL(mousePress(QMouseEvent*)), this, SLOT(mousePress(QMouseEvent*)));
    connect(ui->widgetDataPlot, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(mouseMove(QMouseEvent*)));
    connect(ui->widgetDataPlot, SIGNAL(mouseRelease(QMouseEvent*)), this, SLOT(mouseRelease(QMouseEvent*)));

    connect(ui->widgetFFTPlot, SIGNAL(mousePress(QMouseEvent*)), this, SLOT(mousePress(QMouseEvent*)));
    connect(ui->widgetFFTPlot, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(mouseMove(QMouseEvent*)));
    connect(ui->widgetFFTPlot, SIGNAL(mouseRelease(QMouseEvent*)), this, SLOT(mouseRelease(QMouseEvent*)));

    connect(ui->radioButtonRMS, SIGNAL(clicked()), SLOT(spectrumTypeChanged()));
    connect(ui->radioButtonPSD, SIGNAL(clicked()), SLOT(spectrumTypeChanged()));
    connect(ui->radioButtonESD, SIGNAL(clicked()), SLOT(spectrumTypeChanged()));

    connect(ui->radioButtonXLin, SIGNAL(clicked()), SLOT(XAxisScaleType()));
    connect(ui->radioButtonXLog, SIGNAL(clicked()), SLOT(XAxisScaleType()));

    connect(ui->radioButtonYLin, SIGNAL(clicked()), SLOT(YAxisScaleType()));
    connect(ui->radioButtonYLog, SIGNAL(clicked()), SLOT(YAxisScaleType()));


    ui->widgetDataPlot->addGraph();
    ui->widgetDataPlot->graph(0);

    ui->widgetFFTPlot->addGraph();
    ui->widgetFFTPlot->graph(0);

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::openFile()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Load data"), "", tr("Files (*.dat *.txt)"));

    timeVector.clear();
    dataVector.clear();

    double time = 0.0;
    if (!fileName.isEmpty())
    {
        QFile inFile(fileName);

        if (inFile.open(QIODevice::ReadOnly))
        {
            while (!inFile.atEnd())
            {
                QString line = inFile.readLine();
                timeVector.append(time);
                dataVector.append(line.toDouble());
                time += dt;
            }
        }
        inFile.close();

        ui->widgetDataPlot->graph(0)->setData(timeVector, dataVector);
        ui->widgetDataPlot->rescaleAxes();
        ui->widgetDataPlot->replot();
    }
    calculateFFT(0.0, time);
}


void MainWindow::openMSAcc()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Load data"), "", tr("Files (*.txt)"));

    timeVector.clear();
    dataVector.clear();

    double time = 0.0;

    if (!fileName.isEmpty())
    {
        QFile inFile(fileName);

        if (inFile.open(QIODevice::ReadOnly))
        {
            while (!inFile.atEnd())
            {
                QString line = inFile.readLine();
                line = line.simplified();
                QStringList pieces = line.split(" ");
                timeVector.append(pieces.at(0).toDouble());
                dataVector.append(pieces.at(1).toDouble());
                time += dt;
            }
        }
        inFile.close();

        ui->widgetDataPlot->graph(0)->setData(timeVector, dataVector);
        ui->widgetDataPlot->rescaleAxes();
        ui->widgetDataPlot->replot();
    }
    calculateFFT(0.0, time);
}


/*

void MainWindow::calculateFFT()
{
    fftVector.clear();
    fftVector.resize(0);
    freqVector.clear();
    freqVector.resize(0);

    dt = (timeVector.last() - timeVector.first())/timeVector.length();


    //kiss_fft_cpx outData[dataVector.length()/2+1];
    //kiss_fft_scalar inData[dataVector.length()];
    kiss_fft_scalar* inData;
    inData = (kiss_fft_scalar*)malloc(dataVector.length()*sizeof(kiss_fft_scalar));

    kiss_fft_cpx* outData;

    outData = (kiss_fft_cpx*)malloc((dataVector.length()/2 + 1) * sizeof(kiss_fft_cpx));

    for (int i = 0; i < dataVector.length(); i++)
    {
        inData[i] = dataVector[i];
    }

    kiss_fftr_cfg cfg = kiss_fftr_alloc(dataVector.length(),0, NULL, NULL);

    kiss_fftr(cfg, inData, outData);

    double T = dt*dataVector.length();
    double df = 1/T;
    double f = 0;

    for (int i = 0; i < dataVector.length()/2+1; i++)
    {
         fftVector.append(sqrtf(outData[i].r*outData[i].r + outData[i].i*outData[i].i)*sqrtf(2)/dataVector.length());

         freqVector.append(f);
         f += df;

    }


    ui->widgetFFTPlot->graph(0)->setData(freqVector, fftVector);


    char xLabel[1024];
    sprintf(xLabel, "Frequency [Hz]   df = %1.3f", df);
    ui->widgetFFTPlot->xAxis->setLabel(xLabel);

    ui->widgetFFTPlot->rescaleAxes();
    ui->widgetFFTPlot->replot();

    free(inData);
    free(outData);

}
*/

void MainWindow::selectionChanged()
{

    // make top and bottom axes be selected synchronously, and handle axis and tick labels as one selectable object:
    if (ui->widgetDataPlot->xAxis->selectedParts().testFlag(QCPAxis::spAxis) ||
            ui->widgetDataPlot->xAxis->selectedParts().testFlag(QCPAxis::spTickLabels) ||
            ui->widgetDataPlot->xAxis2->selectedParts().testFlag(QCPAxis::spAxis) ||
            ui->widgetDataPlot->xAxis2->selectedParts().testFlag(QCPAxis::spTickLabels))
    {
        ui->widgetDataPlot->xAxis2->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
        ui->widgetDataPlot->xAxis->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
    }

    if (ui->widgetFFTPlot->xAxis->selectedParts().testFlag(QCPAxis::spAxis) ||
            ui->widgetFFTPlot->xAxis->selectedParts().testFlag(QCPAxis::spTickLabels) ||
            ui->widgetFFTPlot->xAxis2->selectedParts().testFlag(QCPAxis::spAxis) ||
            ui->widgetFFTPlot->xAxis2->selectedParts().testFlag(QCPAxis::spTickLabels))
    {
        ui->widgetFFTPlot->xAxis2->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
        ui->widgetFFTPlot->xAxis->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
    }

    // make left and right axes be selected synchronously, and handle axis and tick labels as one selectable object:
    if (ui->widgetDataPlot->yAxis->selectedParts().testFlag(QCPAxis::spAxis) ||
            ui->widgetDataPlot->yAxis->selectedParts().testFlag(QCPAxis::spTickLabels) ||
            ui->widgetDataPlot->yAxis2->selectedParts().testFlag(QCPAxis::spAxis) ||
            ui->widgetDataPlot->yAxis2->selectedParts().testFlag(QCPAxis::spTickLabels))
    {
        ui->widgetDataPlot->yAxis2->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
        ui->widgetDataPlot->yAxis->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
    }
    if (ui->widgetFFTPlot->yAxis->selectedParts().testFlag(QCPAxis::spAxis) ||
            ui->widgetFFTPlot->yAxis->selectedParts().testFlag(QCPAxis::spTickLabels) ||
            ui->widgetFFTPlot->yAxis2->selectedParts().testFlag(QCPAxis::spAxis) ||
            ui->widgetFFTPlot->yAxis2->selectedParts().testFlag(QCPAxis::spTickLabels))
    {
        ui->widgetFFTPlot->yAxis2->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
        ui->widgetFFTPlot->yAxis->setSelectedParts(QCPAxis::spAxis|QCPAxis::spTickLabels);
    }

    // synchronize selection of graphs with selection of corresponding legend items:
    for (int i=0; i<ui->widgetDataPlot->graphCount(); ++i)
    {
        QCPGraph *graph = ui->widgetDataPlot->graph(i);
        QCPPlottableLegendItem *item = ui->widgetDataPlot->legend->itemWithPlottable(graph);
        if (item->selected() || graph->selected())
        {
            item->setSelected(true);
            graph->setSelected(true);
        }
    }

    // synchronize selection of graphs with selection of corresponding legend items:
    for (int i=0; i<ui->widgetFFTPlot->graphCount(); ++i)
    {
        QCPGraph *graph = ui->widgetFFTPlot->graph(i);
        QCPPlottableLegendItem *item = ui->widgetFFTPlot->legend->itemWithPlottable(graph);
        if (item->selected() || graph->selected())
        {
            item->setSelected(true);
            graph->setSelected(true);
        }
    }
}

void MainWindow::mousePress(QMouseEvent* mevent)
{
    // if an axis is selected, only allow the direction of that axis to be dragged
    // if no axis is selected, both directions may be dragged

    if (ui->widgetDataPlot->xAxis->selectedParts().testFlag(QCPAxis::spAxis))
        ui->widgetDataPlot->axisRect()->setRangeDrag(ui->widgetDataPlot->xAxis->orientation());
    else if (ui->widgetDataPlot->yAxis->selectedParts().testFlag(QCPAxis::spAxis))
        ui->widgetDataPlot->axisRect()->setRangeDrag(ui->widgetDataPlot->yAxis->orientation());
    else
        ui->widgetDataPlot->axisRect()->setRangeDrag(Qt::Horizontal|Qt::Vertical);


    if (ui->widgetFFTPlot->xAxis->selectedParts().testFlag(QCPAxis::spAxis))
        ui->widgetFFTPlot->axisRect()->setRangeDrag(ui->widgetFFTPlot->xAxis->orientation());
    else if (ui->widgetFFTPlot->yAxis->selectedParts().testFlag(QCPAxis::spAxis))
        ui->widgetFFTPlot->axisRect()->setRangeDrag(ui->widgetFFTPlot->yAxis->orientation());
    else
        ui->widgetFFTPlot->axisRect()->setRangeDrag(Qt::Horizontal|Qt::Vertical);
}

void MainWindow::mouseWheel(QWheelEvent* mevent)
{
    if (ui->widgetDataPlot->xAxis->selectedParts().testFlag(QCPAxis::spAxis))
        ui->widgetDataPlot->axisRect()->setRangeZoom(ui->widgetDataPlot->xAxis->orientation());
    else if (ui->widgetDataPlot->yAxis->selectedParts().testFlag(QCPAxis::spAxis))
        ui->widgetDataPlot->axisRect()->setRangeZoom(ui->widgetDataPlot->yAxis->orientation());
    else
        ui->widgetDataPlot->axisRect()->setRangeZoom(Qt::Horizontal|Qt::Vertical);


    if (ui->widgetFFTPlot->xAxis->selectedParts().testFlag(QCPAxis::spAxis))
        ui->widgetFFTPlot->axisRect()->setRangeZoom(ui->widgetDataPlot->xAxis->orientation());
    else if (ui->widgetFFTPlot->yAxis->selectedParts().testFlag(QCPAxis::spAxis))
        ui->widgetFFTPlot->axisRect()->setRangeZoom(ui->widgetDataPlot->yAxis->orientation());
    else
        ui->widgetFFTPlot->axisRect()->setRangeZoom(Qt::Horizontal|Qt::Vertical);
}

void MainWindow::mouseRelease(QMouseEvent *mevent)
{
    if (ui->widgetDataPlot->xAxis->selectedParts().testFlag(QCPAxis::spAxis))
        ui->widgetDataPlot->axisRect()->setRangeDrag(ui->widgetDataPlot->xAxis->orientation());
    else
        ui->widgetDataPlot->axisRect()->setRangeDrag(Qt::Horizontal|Qt::Vertical);

    double xMin = ui->widgetDataPlot->xAxis->range().center() - ui->widgetDataPlot->xAxis->range().size()/2.0;
    double xMax = ui->widgetDataPlot->xAxis->range().center() + ui->widgetDataPlot->xAxis->range().size()/2.0;

    qDebug() << xMin << xMax;

    calculateFFT(xMin, xMax);
}

void MainWindow::mouseMove(QMouseEvent *mevent)
{

}


void MainWindow::calculateFFT(double xMin, double xMax)
{
    fftVector.clear();
    fftVector.resize(0);
    freqVector.clear();
    freqVector.resize(0);

    int indexStart = xMin/dt;
    if (indexStart < 0)
    {
        indexStart = 0;
    }

    int indexStop = xMax/dt;
    if (indexStop > dataVector.length())
    {
        indexStop = dataVector.length();
    }

    if ((indexStop-indexStart) % 2 == 1)
    {
        indexStop -=1;
    }

    int noOfSample = indexStop-indexStart;

    qDebug() << indexStart << indexStop << noOfSample;
    if (noOfSample<2)
    {
        return;
    }


    kiss_fft_scalar* inData;
    inData = (kiss_fft_scalar*)malloc(dataVector.length()*sizeof(kiss_fft_scalar));

    kiss_fft_cpx* outData;

    outData = (kiss_fft_cpx*)malloc((dataVector.length()/2 + 1) * sizeof(kiss_fft_cpx));



    for (int i = indexStart; i < indexStop; i++)
    {
        inData[i-indexStart] = dataVector[i];
    }

    kiss_fftr_cfg cfg = kiss_fftr_alloc(noOfSample,0, NULL, NULL);

    kiss_fftr(cfg, inData, outData);

    double T = dt*noOfSample;
    double df = 1/T;
    double f = 0;
    char yLabel[1024];

    // RMS
    if (m_iSpectrumSelected == 0)
    {

        for (int i = 0; i < noOfSample/2+1; i++)
        {
            fftVector.append(sqrtf(outData[i].r*outData[i].r + outData[i].i*outData[i].i)*sqrtf(2)/noOfSample);
            freqVector.append(f);
            f += df;
        }
        sprintf(yLabel, "Unit");
        ui->widgetFFTPlot->yAxis->setLabel(yLabel);

    }
    // PSD
    else if (m_iSpectrumSelected == 1)
    {
        for (int i = 0; i < noOfSample/2+1; i++)
        {
            fftVector.append(sqrtf(outData[i].r*outData[i].r + outData[i].i*outData[i].i)*sqrtf(2)/noOfSample);
            fftVector.replace(i, (fftVector.at(i)*fftVector.at(i)*dt*noOfSample));
            /*
             if (i > 0)
             {

                 fftVector.replace(i, (fftVector.at(i-1) + fftVector.at(i)));
             }
*/
            freqVector.append(f);
            f += df;
        }

        sprintf(yLabel, "Unit^2/Hz");
        ui->widgetFFTPlot->yAxis->setLabel(yLabel);
    }
    // ESD
    else if (m_iSpectrumSelected == 2)
    {
        for (int i = 0; i < noOfSample/2+1; i++)
        {
            fftVector.append(sqrtf(outData[i].r*outData[i].r + outData[i].i*outData[i].i)*sqrtf(2)/noOfSample);



            fftVector.replace(i, (fftVector.at(i)*fftVector.at(i)*dt*noOfSample));

            // Cumsum
            if (i > 0)
            {

                fftVector.replace(i, (fftVector.at(i-1)) + (fftVector.at(i)));
            }

            freqVector.append(f);
            f += df;
        }
        for (int i = 0; i < noOfSample/2+1; i++)
        {
            fftVector.replace(i, sqrt(fftVector.at(i)));
        }
    }



    ui->widgetFFTPlot->graph(0)->setData(freqVector, fftVector);
    char xLabel[1024];
    sprintf(xLabel, "Frequency [Hz]   df = %1.3f", df);
    ui->widgetFFTPlot->xAxis->setLabel(xLabel);

    ui->widgetFFTPlot->rescaleAxes();
    ui->widgetFFTPlot->replot();

    free(inData);
    free(outData);


}

void MainWindow::spectrumTypeChanged()
{
    if (ui->radioButtonRMS->isChecked())
    {
        m_iSpectrumSelected = 0;
        qDebug() << "0";
    }
    else if (ui->radioButtonPSD->isChecked())
    {
        m_iSpectrumSelected = 1;
        qDebug() << "1";
    }
    else if (ui->radioButtonESD->isChecked())
    {
        m_iSpectrumSelected = 2;
        qDebug() << "2";
    }
    else
    {
        m_iSpectrumSelected = 0;
        qDebug() << "Nop";
    }
}

void MainWindow::XAxisScaleType()
{
    if (ui->radioButtonXLin->isChecked())
    {
        ui->widgetFFTPlot->xAxis->setScaleType(QCPAxis::stLinear);
        qDebug() << "X lin";
    }
    else if (ui->radioButtonXLog->isChecked())
    {
        ui->widgetFFTPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);
        qDebug() << "X log";
    }
    ui->widgetFFTPlot->rescaleAxes();
    ui->widgetFFTPlot->replot();
}

void MainWindow::YAxisScaleType()
{
    if (ui->radioButtonYLin->isChecked())
    {
        ui->widgetFFTPlot->yAxis->setScaleType(QCPAxis::stLinear);
        qDebug() << "Y lin";
    }
    else if (ui->radioButtonYLog->isChecked())
    {
        ui->widgetFFTPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
        qDebug() << "Y log";
    }
    ui->widgetFFTPlot->rescaleAxes();
    ui->widgetFFTPlot->replot();

}
