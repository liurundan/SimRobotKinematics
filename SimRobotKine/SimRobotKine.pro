TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        DHParam.cpp \
        Kinematics.cpp \
        Utils.cpp \
        main.cpp

HEADERS += \
    DHParam.h \
    Kinematics.h \
    Utils.h
