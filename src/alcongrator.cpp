#include <QApplication>
#include <QFile>
#include <QStyle>
#include "MainWindow.h"

int main(int argc, char* argv[]) {
    const QString title = "Alcongrator";

    QApplication app(argc, argv);
    QApplication::setApplicationName(title);
    QApplication::setApplicationDisplayName(title);
    QApplication::setOrganizationName("Humboldt-Universität zu Berlin");
    QApplication::setWindowIcon(QIcon(":/images/aci.ico"));

    MainWindow win;

    win.setWindowTitle(title);
    win.showMaximized();

    return QApplication::exec();
}
