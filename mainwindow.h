#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    QVector<double> dataVector;
    QVector<double> timeVector;

    QVector<double> fftVector;
    QVector<double> freqVector;

    int m_iSpectrumSelected;



private slots:
    void openFile();
    void openMSAcc();
    //void calculateFFT();
    void calculateFFT(double xMin, double xMax);
    void XAxisScaleType();
    void YAxisScaleType();


    void mousePress(QMouseEvent* mevent);
    void mouseMove(QMouseEvent *mevent);
    void mouseRelease(QMouseEvent *mevent);
    void mouseWheel(QWheelEvent* mevent);
    void selectionChanged();
    void spectrumTypeChanged();
};

#endif // MAINWINDOW_H
