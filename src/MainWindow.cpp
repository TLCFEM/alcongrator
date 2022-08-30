#include "MainWindow.h"
#include "resampling.h"
#include "ui_MainWindow.h"

#include <CubicSpline.h>
#include <IntegrationScheme.h>
#include <QAudioDeviceInfo>
#include <QAudioFormat>
#include <QAudioOutput>

QVector<double> to_vector(const arma::vec& in) {
    QVector<double> out(int(in.n_elem), 0.);
    for(auto I = 0; I < int(in.n_elem); ++I) out[I] = in(I);
    return out;
}

arma::uword nextpow2(const arma::uword in) { return arma::uword(std::ceil(log2(double(in)))); }

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
    set_label();
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::on_load_data_clicked() {
    const QString name = QFileDialog::getOpenFileName(this, tr("Load Data"), "", tr("All Files (*)"));

    if(name.isEmpty())
        return;

    load_data(name);

    if(source_data.empty()) return;

    plot_time_curve(ui->source_canvas, source_data, "Source Amplitude");

    ui->step_size->setRange(.0001, source_data.col(0).max() / std::min(source_data.n_elem, 2llu));
    ui->step_size->setValue(arma::diff(source_data.col(0)).min());

    on_refresh_clicked();
}

void MainWindow::on_save_clicked() {
    if(time.empty()) return;

    const QString name = QFileDialog::getSaveFileName(this, tr("Save to"), "", tr("All Files (*)"));

    if(!name.isEmpty())
        save_data(name);
}

void MainWindow::on_scheme_currentTextChanged(const QString&) {
    set_label();
    on_refresh_clicked();
}

void MainWindow::on_acceleration_clicked(const bool checked) {
    ui->v_offset->setEnabled(checked);
    ui->u_offset->setEnabled(checked);
    on_refresh_clicked();
}

void MainWindow::on_velocity_clicked(const bool checked) {
    ui->v_offset->setDisabled(checked);
    ui->u_offset->setEnabled(checked);
    on_refresh_clicked();
}

void MainWindow::on_displacement_clicked(const bool checked) {
    ui->v_offset->setDisabled(checked);
    ui->u_offset->setDisabled(checked);
    on_refresh_clicked();
}

void MainWindow::on_apply_filter_clicked() {
    on_refresh_clicked();
}

void MainWindow::on_filter_currentTextChanged(const QString&) {
    on_refresh_clicked();
}

void MainWindow::on_window_currentTextChanged(const QString&) {
    on_refresh_clicked();
}

void MainWindow::on_window_length_valueChanged(int) {
    on_refresh_clicked();
}

void MainWindow::on_low_bound_valueChanged(double) {
    on_refresh_clicked();
}

void MainWindow::on_high_bound_valueChanged(double) {
    on_refresh_clicked();
}

void MainWindow::on_spline_clicked() {
    on_refresh_clicked();
}

void MainWindow::on_natural_clicked() {
    on_refresh_clicked();
}

void MainWindow::on_refresh_clicked() {
    on_step_size_valueChanged(ui->step_size->value());
}

void MainWindow::interpolate(arma::vec& result) {
    if(!ui->spline->isChecked()) {
        arma::interp1(source_data.col(0), source_data.col(1), time, result, "*linear", 0.);
        return;
    }

    CubicSpline interp(source_data.col(0), source_data.col(1), ui->natural->isChecked());

    result.set_size(time.n_elem);

    for(auto I = 0llu; I < time.n_elem; ++I)
        result(I) = interp.evaluate(time(I));
}

void MainWindow::on_step_size_valueChanged(const double step_size) {
    if(source_data.empty()) return;

    time = arma::regspace(0, step_size, source_data.col(0).max());

    if(ui->acceleration->isChecked())
        interpolate(acceleration);
    else if(ui->velocity->isChecked())
        interpolate(velocity);
    else
        interpolate(displacement);

    update_data();
}

void MainWindow::on_pa_valueChanged(double) {
    update_data();
}

void MainWindow::on_pb_valueChanged(double) {
    // if(ui->scheme->currentText() == "Newmark") ui->pa->setRange(.25 * std::pow(ui->pb->value(), 2.), 1.);
    update_data();
}

void MainWindow::on_pc_valueChanged(double) {
    update_data();
}

void MainWindow::on_v_offset_valueChanged(double) {
    update_data();
}

void MainWindow::on_u_offset_valueChanged(double) {
    update_data();
}

void MainWindow::load_data(const QString& name) {
    auto loaded = false;
    if(source_data.load(name.toStdString()) && source_data.n_cols == 2llu)
        loaded = true;
    else if(arma::Col<int> single_column; single_column.load(name.toStdString())) {
        loaded = true;
        source_data.zeros(single_column.n_elem - 2llu, 2);
        for(auto I = 0llu, J = 2llu; I < source_data.n_rows; ++I, ++J) {
            source_data(I, 0) = double(I);
            source_data(I, 1) = single_column(J);
        }
        source_data.col(0) *= 1E-3 * double(single_column(0));
    }

    if(!loaded) {
        QMessageBox::critical(this, tr("Error"), tr("Cannot read data in either one column format or two column format."));
        source_data.reset();
        return;
    }

    source_data.col(1) /= std::max(std::fabs(source_data.col(1).max()), std::fabs(source_data.col(1).min()));
}

void MainWindow::save_data(const QString& name) {
    const auto a_fft = perform_transform(acceleration, true);
    const auto v_fft = perform_transform(velocity, true);
    const auto u_fft = perform_transform(displacement, true);

    const arma::mat time_data_to_export = arma::join_rows(time, acceleration, velocity, displacement);

    arma::field<std::string> header(time_data_to_export.n_cols);
    header(0) = "time";
    header(1) = "acceleration";
    header(2) = "velocity";
    header(3) = "displacement";

    time_data_to_export.save(arma::csv_name((name + "_time").toStdString(), header));

    header(0) = "frequency";

    const arma::mat frequency_data_to_export = arma::join_rows(a_fft, v_fft.col(1), u_fft.col(1));

    frequency_data_to_export.save(arma::csv_name((name + "_frequency").toStdString(), header));
}

void MainWindow::update_data() {
    if(time.empty()) return;

    const auto time_diff = arma::diff(time);
    if(time_diff.is_empty()) return;

    const auto step_size = arma::diff(time).min();
    if(step_size == 0.) return;

    std::unique_ptr<IntegrationScheme> integrator;

    if(ui->scheme->currentText() == "Newmark") {
        const auto beta = ui->pa->value();
        const auto gamma = ui->pb->value();
        integrator = std::make_unique<Newmark>(step_size, beta, gamma);
    }
    else if(ui->scheme->currentText() == "BatheTwoStep")
        integrator = std::make_unique<BatheTwoStep>(step_size);
    else if(ui->scheme->currentText() == "GeneralizedAlpha")
        integrator = std::make_unique<GeneralizedAlpha>(step_size, ui->pa->value());
    else {
        arma::vec R{ui->pa->value(), ui->pb->value(), ui->pc->value()};
        if(ui->scheme->currentText() == "GSSSS-U0")
            integrator = std::make_unique<GSSSSU0>(step_size, std::move(R));
        else if(ui->scheme->currentText() == "GSSSS-V0")
            integrator = std::make_unique<GSSSSV0>(step_size, std::move(R));
    }

    arma::vec h;

    if(ui->apply_filter->isChecked()) {
        const auto fn = .5 / step_size;
        const auto fb = std::min(1., ui->high_bound->value() / fn);
        const auto fa = std::max(0., ui->low_bound->value() / fn);

        if(fb <= fa) return;

        if(ui->filter->currentText() == "Low Pass") {
            if(ui->window->currentText() == "Hamming")
                h = fir_low_pass<WindowType::Hamming>(2 * ui->window_length->value(), fb);
            else if(ui->window->currentText() == "Hann")
                h = fir_low_pass<WindowType::Hann>(2 * ui->window_length->value(), fb);
            else if(ui->window->currentText() == "Blackman")
                h = fir_low_pass<WindowType::Blackman>(2 * ui->window_length->value(), fb);
            else if(ui->window->currentText() == "BlackmanNuttall")
                h = fir_low_pass<WindowType::BlackmanNuttall>(2 * ui->window_length->value(), fb);
            else if(ui->window->currentText() == "BlackmanHarris")
                h = fir_low_pass<WindowType::BlackmanHarris>(2 * ui->window_length->value(), fb);
            else if(ui->window->currentText() == "FlatTop")
                h = fir_low_pass<WindowType::FlatTop>(2 * ui->window_length->value(), fb);
        }
        else if(ui->filter->currentText() == "High Pass") {
            if(ui->window->currentText() == "Hamming")
                h = fir_high_pass<WindowType::Hamming>(2 * ui->window_length->value(), fa);
            else if(ui->window->currentText() == "Hann")
                h = fir_high_pass<WindowType::Hann>(2 * ui->window_length->value(), fa);
            else if(ui->window->currentText() == "Blackman")
                h = fir_high_pass<WindowType::Blackman>(2 * ui->window_length->value(), fa);
            else if(ui->window->currentText() == "BlackmanNuttall")
                h = fir_high_pass<WindowType::BlackmanNuttall>(2 * ui->window_length->value(), fa);
            else if(ui->window->currentText() == "BlackmanHarris")
                h = fir_high_pass<WindowType::BlackmanHarris>(2 * ui->window_length->value(), fa);
            else if(ui->window->currentText() == "FlatTop")
                h = fir_high_pass<WindowType::FlatTop>(2 * ui->window_length->value(), fa);
        }
        else if(ui->filter->currentText() == "Band Pass") {
            if(ui->window->currentText() == "Hamming")
                h = fir_band_pass<WindowType::Hamming>(2 * ui->window_length->value(), fa, fb);
            else if(ui->window->currentText() == "Hann")
                h = fir_band_pass<WindowType::Hann>(2 * ui->window_length->value(), fa, fb);
            else if(ui->window->currentText() == "Blackman")
                h = fir_band_pass<WindowType::Blackman>(2 * ui->window_length->value(), fa, fb);
            else if(ui->window->currentText() == "BlackmanNuttall")
                h = fir_band_pass<WindowType::BlackmanNuttall>(2 * ui->window_length->value(), fa, fb);
            else if(ui->window->currentText() == "BlackmanHarris")
                h = fir_band_pass<WindowType::BlackmanHarris>(2 * ui->window_length->value(), fa, fb);
            else if(ui->window->currentText() == "FlatTop")
                h = fir_band_pass<WindowType::FlatTop>(2 * ui->window_length->value(), fa, fb);
        }
        else if(ui->filter->currentText() == "Band Stop") {
            if(ui->window->currentText() == "Hamming")
                h = fir_band_stop<WindowType::Hamming>(2 * ui->window_length->value(), fa, fb);
            else if(ui->window->currentText() == "Hann")
                h = fir_band_stop<WindowType::Hann>(2 * ui->window_length->value(), fa, fb);
            else if(ui->window->currentText() == "Blackman")
                h = fir_band_stop<WindowType::Blackman>(2 * ui->window_length->value(), fa, fb);
            else if(ui->window->currentText() == "BlackmanNuttall")
                h = fir_band_stop<WindowType::BlackmanNuttall>(2 * ui->window_length->value(), fa, fb);
            else if(ui->window->currentText() == "BlackmanHarris")
                h = fir_band_stop<WindowType::BlackmanHarris>(2 * ui->window_length->value(), fa, fb);
            else if(ui->window->currentText() == "FlatTop")
                h = fir_band_stop<WindowType::FlatTop>(2 * ui->window_length->value(), fa, fb);
        }
    }

    if(ui->acceleration->isChecked()) {
        velocity.zeros(time.n_elem);
        displacement.zeros(time.n_elem);
        acceleration = apply_filter(acceleration, h);

        velocity(0) = ui->v_offset->value();
        displacement(0) = ui->u_offset->value();

        integrator->update_from_acceleration(displacement, velocity, acceleration);
    }
    else if(ui->velocity->isChecked()) {
        acceleration.zeros(time.n_elem);
        displacement.zeros(time.n_elem);
        velocity = apply_filter(velocity, h);

        displacement(0) = ui->u_offset->value();

        integrator->update_from_velocity(displacement, velocity, acceleration);
    }
    else {
        acceleration.zeros(time.n_elem);
        velocity.zeros(time.n_elem);
        displacement = apply_filter(displacement, h);

        integrator->update_from_displacement(displacement, velocity, acceleration);
    }

    on_frequency_clicked(ui->frequency->checkState());
}

void MainWindow::set_label() {
    ui->pa->setRange(0,1);
    ui->pb->setRange(0,1);
    ui->pc->setRange(0,1);

    if(ui->scheme->currentText() == "Newmark") {
        ui->pa_label->setText("beta");
        ui->pa->setValue(.25);
        ui->pa_label->setDisabled(false);
        ui->pa->setDisabled(false);

        ui->pb_label->setText("gamma");
        ui->pb->setValue(.5);
        ui->pb_label->setDisabled(false);
        ui->pb->setDisabled(false);

        ui->pc_label->setText("");
        ui->pc->setValue(0.);
        ui->pc_label->setDisabled(true);
        ui->pc->setDisabled(true);
    }
    else if(ui->scheme->currentText() == "BatheTwoStep") {
        ui->pa_label->setText("");
        ui->pa->setValue(0);
        ui->pa_label->setDisabled(true);
        ui->pa->setDisabled(true);

        ui->pb_label->setText("");
        ui->pb->setValue(0);
        ui->pb_label->setDisabled(true);
        ui->pb->setDisabled(true);

        ui->pc_label->setText("");
        ui->pc->setValue(0);
        ui->pc_label->setDisabled(true);
        ui->pc->setDisabled(true);
    }
    else if(ui->scheme->currentText() == "GeneralizedAlpha") {
        ui->pa_label->setText("spectral radius");
        ui->pa->setValue(.8);
        ui->pa_label->setDisabled(false);
        ui->pa->setDisabled(false);

        ui->pb_label->setText("");
        ui->pb->setValue(0);
        ui->pb_label->setDisabled(true);
        ui->pb->setDisabled(true);

        ui->pc_label->setText("");
        ui->pc->setValue(0);
        ui->pc_label->setDisabled(true);
        ui->pc->setDisabled(true);
    }
    else if(ui->scheme->currentText() == "GSSSS-U0") {
        ui->pa_label->setText("spectral radius 1");
        ui->pa->setValue(1);
        ui->pa_label->setDisabled(false);
        ui->pa->setDisabled(false);

        ui->pb_label->setText("spectral radius 2");
        ui->pb->setValue(1);
        ui->pb_label->setDisabled(false);
        ui->pb->setDisabled(false);

        ui->pc_label->setText("spectral radius 3");
        ui->pc->setValue(0);
        ui->pc_label->setDisabled(false);
        ui->pc->setDisabled(true);
    }
    else if(ui->scheme->currentText() == "GSSSS-V0") {
        ui->pa_label->setText("spectral radius 1");
        ui->pa->setValue(1);
        ui->pa_label->setDisabled(false);
        ui->pa->setDisabled(true);

        ui->pb_label->setText("spectral radius 2");
        ui->pb->setValue(1);
        ui->pb_label->setDisabled(false);
        ui->pb->setDisabled(true);

        ui->pc_label->setText("spectral radius 3");
        ui->pc->setValue(1);
        ui->pc_label->setDisabled(false);
        ui->pc->setDisabled(false);
    }
}

void MainWindow::replot() {
    ui->source_canvas->setBackground(background_color);
    ui->target_canvas->setBackground(background_color);
    ui->source_canvas->replot();
    ui->target_canvas->replot();
}

void MainWindow::on_freq_a_clicked() {
    on_frequency_clicked(ui->frequency->checkState());
}

void MainWindow::on_freq_v_clicked() {
    on_frequency_clicked(ui->frequency->checkState());
}

void MainWindow::on_freq_u_clicked() {
    on_frequency_clicked(ui->frequency->checkState());
}

void MainWindow::on_frequency_clicked(const bool checked) {
    if(time.empty()) return;

    if(!checked) {
        if(ui->freq_a->isChecked())
            plot_time_curve(ui->target_canvas, arma::join_rows(time, acceleration), "Target Acceleration Amplitude");
        else if(ui->freq_v->isChecked())
            plot_time_curve(ui->target_canvas, arma::join_rows(time, velocity), "Target Velocity Amplitude");
        else if(ui->freq_u->isChecked())
            plot_time_curve(ui->target_canvas, arma::join_rows(time, displacement), "Target Displacement Amplitude");
    }
    else if(ui->freq_a->isChecked())
        plot_frequency_curve(ui->target_canvas, perform_transform(acceleration), "Target Acceleration Amplitude");
    else if(ui->freq_v->isChecked())
        plot_frequency_curve(ui->target_canvas, perform_transform(velocity), "Target Velocity Amplitude");
    else if(ui->freq_u->isChecked())
        plot_frequency_curve(ui->target_canvas, perform_transform(displacement), "Target Displacement Amplitude");

    replot();
}

arma::mat MainWindow::perform_transform(const arma::vec& data, const bool full) {
    const auto length = 2 << nextpow2(time.n_elem);
    const arma::vec time_diff = arma::diff(time);
    if(time_diff.empty()) return {};
    const auto step_size = time_diff.min();
    if(step_size == 0.) return {};

    const arma::vec fft_frequency = arma::regspace(0, 1, length - 1) / (step_size * double(length));
    arma::cx_vec fft_cx_magnitude = arma::fft(data, arma::uword(length));
    arma::vec fft_magnitude = 2. * arma::abs(fft_cx_magnitude) / double(data.n_elem);
    fft_magnitude(0) /= 2.;

    auto half_length = arma::uword(length) / 2;

    if(full)
        return arma::join_rows(fft_frequency.head(half_length), fft_magnitude.head(half_length));

    const auto cutoff = std::fabs(ui->high_bound->value());
    const arma::uvec head_length = arma::find(fft_frequency > cutoff);

    if(!head_length.empty() && head_length(0) + 1 < half_length) half_length = head_length(0) + 1;

    return arma::join_rows(fft_frequency.head(half_length), fft_magnitude.head(half_length));
}

void MainWindow::plot_time_curve(QCustomPlot* canvas, const arma::mat& data, const char* y_label) {
    if(data.empty()) return;

    initialise_canvas(canvas, "Time", y_label);

    plot_curve(canvas, data.col(0), data.col(1));
}

void MainWindow::plot_frequency_curve(QCustomPlot* canvas, const arma::mat& data, const char* y_label) {
    if(data.empty()) return;

    initialise_canvas(canvas, "Frequency", y_label);

    plot_curve(canvas, data.col(0), data.col(1));
}

void MainWindow::initialise_canvas(QCustomPlot* canvas, const char* x_label, const char* y_label) {
    canvas->setNoAntialiasingOnDrag(true);
    canvas->setInteractions(QCP::iNone);
    canvas->setPlottingHint(QCP::PlottingHint::phImmediateRefresh);
    canvas->setNotAntialiasedElement(QCP::AntialiasedElement::aePlottables);

    canvas->xAxis->setLabel(x_label);
    canvas->yAxis->setLabel(y_label);

    canvas->xAxis->grid()->setSubGridVisible(true);
    canvas->yAxis->grid()->setSubGridVisible(true);
}

void MainWindow::plot_curve(QCustomPlot* canvas, const arma::vec& x_data, const arma::vec& y_data) {
    const auto x_max = x_data.max();
    const auto y_max = y_data.max();
    const auto y_min = y_data.min();

    canvas->clearGraphs();
    canvas->xAxis->setRange(0., 1.05 * x_max);
    canvas->yAxis->setRange(1.1 * y_min, 1.1 * y_max);

    QPen pen;
    pen.setWidth(x_data.size() > 1000 ? 1 : 2);
    pen.setColor(canvas == ui->source_canvas ? source_color : target_color);

    canvas->addGraph();
    canvas->graph()->setPen(pen);
    canvas->graph()->setData(to_vector(x_data), to_vector(y_data));
    canvas->replot();
}

void MainWindow::on_light_clicked(const bool checked) {
    if(checked) {
        QFile file(":/utilities/stylesheet_francesco.qss");
        file.open(QFile::ReadOnly);
        setStyleSheet(QString::fromLatin1(file.readAll()));

        background_color = QColor(131, 134, 137);
    }
    else {
        setStyleSheet("");

        background_color = Qt::white;
    }

    replot();
}

void MainWindow::on_listen_clicked() {
    if(time.empty()) return;

    QAudioFormat audioFormat;
    audioFormat.setSampleRate(static_cast<int>(ui->upsample_ratio->value() / ui->step_size->value()));
    audioFormat.setChannelCount(1);
    audioFormat.setSampleSize(32);
    audioFormat.setCodec("audio/pcm");
    audioFormat.setByteOrder(QAudioFormat::LittleEndian);
    audioFormat.setSampleType(QAudioFormat::Float);

    QAudioDeviceInfo deviceInfo(QAudioDeviceInfo::defaultOutputDevice());
    if(!deviceInfo.isFormatSupported(audioFormat)) return;

    arma::fvec f_data;
    if(ui->freq_a->isChecked())
        f_data = arma::conv_to<arma::fvec>::from(acceleration);
    else if(ui->freq_v->isChecked())
        f_data = arma::conv_to<arma::fvec>::from(velocity);
    else
        f_data = arma::conv_to<arma::fvec>::from(displacement);

    f_data /= std::max(std::fabs(f_data.min()), std::fabs(f_data.max()));

    auto* byteBuffer = new QByteArray((char*)f_data.memptr(), sizeof(float) * acceleration.n_elem);

    auto* input = new QBuffer(byteBuffer);
    input->open(QIODevice::ReadOnly);

    auto* audio = new QAudioOutput(audioFormat, this);

    connect(audio, &QAudioOutput::stateChanged, [audio, input, byteBuffer](QAudio::State newState) {
        if(newState != QAudio::IdleState) return;

        delete audio;
        delete input;
        delete byteBuffer;
    });

    audio->start(input);
}
