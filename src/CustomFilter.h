#ifndef CUSTOMFILTER_H
#define CUSTOMFILTER_H

#include <QDialog>
#include <armadillo>

namespace Ui {
    class CustomFilter;
}

class CustomFilter : public QDialog {
    Q_OBJECT

public:
    explicit CustomFilter(QWidget* = nullptr);
    ~CustomFilter();

    void set_coefficient(const arma::vec& h);
    arma::vec get_coefficient();

private:
    Ui::CustomFilter* ui;
};

#endif // CUSTOMFILTER_H
