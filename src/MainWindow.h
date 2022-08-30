#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <armadillo>
#include <qcustomplot.h>

QT_BEGIN_NAMESPACE
namespace Ui {
    class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget* = nullptr);
    ~MainWindow() override;

private:
    Ui::MainWindow* ui;

    QColor background_color = Qt::white;
    QColor source_color = QColor(202, 0, 32);
    QColor target_color = QColor(5, 113, 176);

    arma::mat source_data;
    arma::vec time, acceleration, velocity, displacement;

    static void initialise_canvas(QCustomPlot*, const char*, const char*);
    void plot_time_curve(QCustomPlot*, const arma::mat&, const char*);
    void plot_frequency_curve(QCustomPlot*, const arma::mat&, const char*);
    void plot_curve(QCustomPlot*, const arma::vec&, const arma::vec&);
    arma::mat perform_transform(const arma::vec&, const bool full = false);
    void update_data();
    void set_label();
    void replot();
    void load_data(const QString&);
    void save_data(const QString&);
    void interpolate(arma::vec&);

private slots:
    void on_load_data_clicked();
    void on_scheme_currentTextChanged(const QString&);
    void on_acceleration_clicked(bool);
    void on_velocity_clicked(bool);
    void on_displacement_clicked(bool);
    void on_apply_filter_clicked();
    void on_filter_currentTextChanged(const QString&);
    void on_window_currentTextChanged(const QString&);
    void on_window_length_valueChanged(int);
    void on_low_bound_valueChanged(double);
    void on_high_bound_valueChanged(double);
    void on_spline_clicked();
    void on_natural_clicked();
    void on_refresh_clicked();
    void on_step_size_valueChanged(double);
    void on_pa_valueChanged(double);
    void on_pb_valueChanged(double);
    void on_pc_valueChanged(double);
    void on_v_offset_valueChanged(double);
    void on_u_offset_valueChanged(double);
    void on_freq_a_clicked();
    void on_freq_v_clicked();
    void on_freq_u_clicked();
    void on_frequency_clicked(bool);
    void on_light_clicked(bool);
    void on_save_clicked();
    void on_listen_clicked();
};

#endif // MAINWINDOW_H
