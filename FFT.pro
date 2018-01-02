#-------------------------------------------------
#
# Project created by QtCreator 2015-09-15T08:00:36
#
#-------------------------------------------------

QT       += core gui printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = FFT
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp\
        qcustomplot.cpp\
        kiss_fft.c\
        kiss_fftr.c

HEADERS  += mainwindow.h\
            qcustomplot.h\
            kiss_fft.h\
            kiss_fftr.h\
            _kiss_fft_guts.h


FORMS    += mainwindow.ui
