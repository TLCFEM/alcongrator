#include "CustomFilter.h"
#include <QIcon>
#include <regex>
#include "ui_CustomFilter.h"

CustomFilter::CustomFilter(QWidget* parent)
    : QDialog(parent), ui(new Ui::CustomFilter) {
    ui->setupUi(this);

    setWindowTitle("Custom Filter");
    setWindowIcon(QIcon(":/images/aci.ico"));
}

CustomFilter::~CustomFilter() { delete ui; }

void CustomFilter::set_coefficient(const arma::vec& h) {
    std::ostringstream output;
    if(!h.empty()) h.print(output);
    ui->coef_text->setPlainText(output.str().c_str());
}

arma::vec CustomFilter::get_coefficient() {
    const auto coef_string = ui->coef_text->toPlainText().toStdString();

    const std::regex scientific{"((\\+|-)?[[:digit:]]+)(\\.(([[:digit:]]+)?))?((e|E)((\\+|-)?)[[:digit:]]+)?"};

    std::vector<double> number;

    for(auto I = std::sregex_iterator(coef_string.begin(), coef_string.end(), scientific); I != std::sregex_iterator(); ++I)
        number.emplace_back(std::stod({I->str()}));

    return number;
}
