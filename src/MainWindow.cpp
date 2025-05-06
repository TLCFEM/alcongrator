#include "MainWindow.h"
#include "resampling.h"
#include "ui_MainWindow.h"

#include <CubicSpline.h>
#include <IntegrationScheme.h>
#include <QAudioFormat>
#include <QAudioOutput>
#include <QMediaPlayer>

QVector<double> to_vector(const arma::vec& in) {
    QVector<double> out(int(in.n_elem), 0.);
    for(auto I = 0; I < int(in.n_elem); ++I) out[I] = in(I);
    return out;
}

arma::uword nextpow2(const arma::uword in) { return arma::uword(std::ceil(log2(double(in)))); }

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent), main_page(new Ui::MainWindow), filter_page(this) {
    main_page->setupUi(this);
    set_label();

    connect(&filter_page, &CustomFilter::accepted, this, &MainWindow::set_coefficient);
}

MainWindow::~MainWindow() {
    delete main_page;
}

void MainWindow::on_load_data_clicked() {
    const QString name = QFileDialog::getOpenFileName(this, tr("Load Data"), "", tr("All Files (*)"));

    if(name.isEmpty())
        return;

    load_data(name);

    if(source_data.empty()) return;

    plot_time_curve(main_page->source_canvas, source_data, "Source Amplitude");

    main_page->step_size->setRange(.0001, source_data.col(0).max() / std::min(source_data.n_elem, 2llu));
    main_page->step_size->setValue(arma::diff(source_data.col(0)).min());

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
    main_page->u_offset->setEnabled(checked);
    main_page->v_offset->setEnabled(checked);
    main_page->a_offset->setDisabled(checked);
    on_refresh_clicked();
}

void MainWindow::on_velocity_clicked(const bool checked) {
    main_page->u_offset->setEnabled(checked);
    main_page->v_offset->setDisabled(checked);
    main_page->a_offset->setEnabled(checked);
    on_refresh_clicked();
}

void MainWindow::on_displacement_clicked(const bool checked) {
    main_page->u_offset->setDisabled(checked);
    main_page->v_offset->setEnabled(checked);
    main_page->a_offset->setEnabled(checked);
    on_refresh_clicked();
}

void MainWindow::on_apply_filter_clicked() {
    on_refresh_clicked();
}

void MainWindow::on_filter_currentTextChanged(const QString& filter_type) {
    if(filter_type == "Custom") {
        main_page->window->setDisabled(true);
        main_page->window_length->setDisabled(true);
        main_page->low_bound->setDisabled(true);
        main_page->high_bound->setDisabled(true);
    }
    else {
        main_page->window->setEnabled(true);
        main_page->window_length->setEnabled(true);
        main_page->low_bound->setEnabled(true);
        main_page->high_bound->setEnabled(true);
    }

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
    on_step_size_valueChanged(main_page->step_size->value());
}

void MainWindow::interpolate(arma::vec& result) {
    if(!main_page->spline->isChecked()) {
        arma::interp1(source_data.col(0), source_data.col(1), time, result, "*linear", 0.);
        return;
    }

    CubicSpline interp(source_data.col(0), source_data.col(1), main_page->natural->isChecked());

    result.set_size(time.n_elem);

    for(auto I = 0llu; I < time.n_elem; ++I)
        result(I) = interp.evaluate(time(I));
}

void MainWindow::on_step_size_valueChanged(const double step_size) {
    if(source_data.empty()) return;

    time = arma::regspace(0, step_size, source_data.col(0).max());

    if(main_page->acceleration->isChecked())
        interpolate(acceleration);
    else if(main_page->velocity->isChecked())
        interpolate(velocity);
    else
        interpolate(displacement);

    update_data();
}

void MainWindow::on_pa_valueChanged(double) {
    on_refresh_clicked();
}

void MainWindow::on_pb_valueChanged(double) {
    on_refresh_clicked();
}

void MainWindow::on_pc_valueChanged(double) {
    on_refresh_clicked();
}

void MainWindow::on_v_offset_valueChanged(double) {
    on_refresh_clicked();
}

void MainWindow::on_u_offset_valueChanged(double) {
    on_refresh_clicked();
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
    const auto a_fft = perform_transform(acceleration);
    const auto v_fft = perform_transform(velocity);
    const auto u_fft = perform_transform(displacement);

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

    if(main_page->scheme->currentText() == "Newmark") {
        const auto beta = main_page->pa->value();
        const auto gamma = main_page->pb->value();
        integrator = std::make_unique<Newmark>(step_size, beta, gamma);
    }
    else if(main_page->scheme->currentText() == "BatheTwoStep")
        integrator = std::make_unique<BatheTwoStep>(step_size);
    else if(main_page->scheme->currentText() == "GeneralizedAlpha")
        integrator = std::make_unique<GeneralizedAlpha>(step_size, main_page->pa->value());
    else {
        arma::vec R{main_page->pa->value(), main_page->pb->value(), main_page->pc->value()};
        if(main_page->scheme->currentText() == "GSSSS-U0")
            integrator = std::make_unique<GSSSSU0>(step_size, std::move(R));
        else if(main_page->scheme->currentText() == "GSSSS-V0")
            integrator = std::make_unique<GSSSSV0>(step_size, std::move(R));
    }

    arma::vec h;

    if(main_page->apply_filter->isChecked()) {
        const auto fn = .5 / step_size;
        const auto fb = std::min(1., main_page->high_bound->value() / fn);
        const auto fa = std::max(0., main_page->low_bound->value() / fn);

        if(fb <= fa) return;

        if(main_page->filter->currentText() == "Low Pass") {
            if(main_page->window->currentText() == "Hamming")
                h = fir_low_pass<WindowType::Hamming>(2 * main_page->window_length->value(), fb);
            else if(main_page->window->currentText() == "Hann")
                h = fir_low_pass<WindowType::Hann>(2 * main_page->window_length->value(), fb);
            else if(main_page->window->currentText() == "Blackman")
                h = fir_low_pass<WindowType::Blackman>(2 * main_page->window_length->value(), fb);
            else if(main_page->window->currentText() == "BlackmanNuttall")
                h = fir_low_pass<WindowType::BlackmanNuttall>(2 * main_page->window_length->value(), fb);
            else if(main_page->window->currentText() == "BlackmanHarris")
                h = fir_low_pass<WindowType::BlackmanHarris>(2 * main_page->window_length->value(), fb);
            else if(main_page->window->currentText() == "FlatTop")
                h = fir_low_pass<WindowType::FlatTop>(2 * main_page->window_length->value(), fb);
        }
        else if(main_page->filter->currentText() == "High Pass") {
            if(main_page->window->currentText() == "Hamming")
                h = fir_high_pass<WindowType::Hamming>(2 * main_page->window_length->value(), fa);
            else if(main_page->window->currentText() == "Hann")
                h = fir_high_pass<WindowType::Hann>(2 * main_page->window_length->value(), fa);
            else if(main_page->window->currentText() == "Blackman")
                h = fir_high_pass<WindowType::Blackman>(2 * main_page->window_length->value(), fa);
            else if(main_page->window->currentText() == "BlackmanNuttall")
                h = fir_high_pass<WindowType::BlackmanNuttall>(2 * main_page->window_length->value(), fa);
            else if(main_page->window->currentText() == "BlackmanHarris")
                h = fir_high_pass<WindowType::BlackmanHarris>(2 * main_page->window_length->value(), fa);
            else if(main_page->window->currentText() == "FlatTop")
                h = fir_high_pass<WindowType::FlatTop>(2 * main_page->window_length->value(), fa);
        }
        else if(main_page->filter->currentText() == "Band Pass") {
            if(main_page->window->currentText() == "Hamming")
                h = fir_band_pass<WindowType::Hamming>(2 * main_page->window_length->value(), fa, fb);
            else if(main_page->window->currentText() == "Hann")
                h = fir_band_pass<WindowType::Hann>(2 * main_page->window_length->value(), fa, fb);
            else if(main_page->window->currentText() == "Blackman")
                h = fir_band_pass<WindowType::Blackman>(2 * main_page->window_length->value(), fa, fb);
            else if(main_page->window->currentText() == "BlackmanNuttall")
                h = fir_band_pass<WindowType::BlackmanNuttall>(2 * main_page->window_length->value(), fa, fb);
            else if(main_page->window->currentText() == "BlackmanHarris")
                h = fir_band_pass<WindowType::BlackmanHarris>(2 * main_page->window_length->value(), fa, fb);
            else if(main_page->window->currentText() == "FlatTop")
                h = fir_band_pass<WindowType::FlatTop>(2 * main_page->window_length->value(), fa, fb);
        }
        else if(main_page->filter->currentText() == "Band Stop") {
            if(main_page->window->currentText() == "Hamming")
                h = fir_band_stop<WindowType::Hamming>(2 * main_page->window_length->value(), fa, fb);
            else if(main_page->window->currentText() == "Hann")
                h = fir_band_stop<WindowType::Hann>(2 * main_page->window_length->value(), fa, fb);
            else if(main_page->window->currentText() == "Blackman")
                h = fir_band_stop<WindowType::Blackman>(2 * main_page->window_length->value(), fa, fb);
            else if(main_page->window->currentText() == "BlackmanNuttall")
                h = fir_band_stop<WindowType::BlackmanNuttall>(2 * main_page->window_length->value(), fa, fb);
            else if(main_page->window->currentText() == "BlackmanHarris")
                h = fir_band_stop<WindowType::BlackmanHarris>(2 * main_page->window_length->value(), fa, fb);
            else if(main_page->window->currentText() == "FlatTop")
                h = fir_band_stop<WindowType::FlatTop>(2 * main_page->window_length->value(), fa, fb);
        }
        else if(main_page->filter->currentText() == "Custom") {
            if(custom_filter.empty()) {
                QMessageBox::information(this, tr("Need coefficients"), tr("Please set filter coefficients first."));
                return;
            }
            h = custom_filter;
        }
    }

    if(main_page->acceleration->isChecked()) {
        velocity.zeros(time.n_elem);
        displacement.zeros(time.n_elem);
        acceleration = apply_filter(acceleration, h);

        velocity(0) = main_page->v_offset->value();
        displacement(0) = main_page->u_offset->value();

        integrator->update_from_acceleration(displacement, velocity, acceleration);
    }
    else if(main_page->velocity->isChecked()) {
        acceleration.zeros(time.n_elem);
        displacement.zeros(time.n_elem);
        velocity = apply_filter(velocity, h);

        acceleration(0) = main_page->a_offset->value();
        displacement(0) = main_page->u_offset->value();

        integrator->update_from_velocity(displacement, velocity, acceleration);
    }
    else {
        acceleration.zeros(time.n_elem);
        velocity.zeros(time.n_elem);
        displacement = apply_filter(displacement, h);

        acceleration(0) = main_page->a_offset->value();
        velocity(0) = main_page->v_offset->value();

        integrator->update_from_displacement(displacement, velocity, acceleration);
    }

    on_frequency_clicked(main_page->frequency->checkState());
}

void MainWindow::set_label() {
    main_page->pa->setRange(0, 1);
    main_page->pb->setRange(0, 1);
    main_page->pc->setRange(0, 1);

    if(main_page->scheme->currentText() == "Newmark") {
        main_page->pa_label->setText("beta");
        main_page->pa->setValue(.25);
        main_page->pa_label->setDisabled(false);
        main_page->pa->setDisabled(false);

        main_page->pb_label->setText("gamma");
        main_page->pb->setValue(.5);
        main_page->pb_label->setDisabled(false);
        main_page->pb->setDisabled(false);

        main_page->pc_label->setText("");
        main_page->pc->setValue(0.);
        main_page->pc_label->setDisabled(true);
        main_page->pc->setDisabled(true);
    }
    else if(main_page->scheme->currentText() == "BatheTwoStep") {
        main_page->pa_label->setText("");
        main_page->pa->setValue(0);
        main_page->pa_label->setDisabled(true);
        main_page->pa->setDisabled(true);

        main_page->pb_label->setText("");
        main_page->pb->setValue(0);
        main_page->pb_label->setDisabled(true);
        main_page->pb->setDisabled(true);

        main_page->pc_label->setText("");
        main_page->pc->setValue(0);
        main_page->pc_label->setDisabled(true);
        main_page->pc->setDisabled(true);
    }
    else if(main_page->scheme->currentText() == "GeneralizedAlpha") {
        main_page->pa_label->setText("spectral radius");
        main_page->pa->setValue(.8);
        main_page->pa_label->setDisabled(false);
        main_page->pa->setDisabled(false);

        main_page->pb_label->setText("");
        main_page->pb->setValue(0);
        main_page->pb_label->setDisabled(true);
        main_page->pb->setDisabled(true);

        main_page->pc_label->setText("");
        main_page->pc->setValue(0);
        main_page->pc_label->setDisabled(true);
        main_page->pc->setDisabled(true);
    }
    else if(main_page->scheme->currentText() == "GSSSS-U0") {
        main_page->pa_label->setText("spectral radius 1");
        main_page->pa->setValue(1);
        main_page->pa_label->setDisabled(false);
        main_page->pa->setDisabled(false);

        main_page->pb_label->setText("spectral radius 2");
        main_page->pb->setValue(1);
        main_page->pb_label->setDisabled(false);
        main_page->pb->setDisabled(false);

        main_page->pc_label->setText("spectral radius 3");
        main_page->pc->setValue(0);
        main_page->pc_label->setDisabled(false);
        main_page->pc->setDisabled(true);
    }
    else if(main_page->scheme->currentText() == "GSSSS-V0") {
        main_page->pa_label->setText("spectral radius 1");
        main_page->pa->setValue(1);
        main_page->pa_label->setDisabled(false);
        main_page->pa->setDisabled(true);

        main_page->pb_label->setText("spectral radius 2");
        main_page->pb->setValue(1);
        main_page->pb_label->setDisabled(false);
        main_page->pb->setDisabled(true);

        main_page->pc_label->setText("spectral radius 3");
        main_page->pc->setValue(1);
        main_page->pc_label->setDisabled(false);
        main_page->pc->setDisabled(false);
    }
}

void MainWindow::replot() {
    main_page->source_canvas->setBackground(background_color);
    main_page->target_canvas->setBackground(background_color);
    main_page->source_canvas->replot();
    main_page->target_canvas->replot();
}

void MainWindow::on_freq_a_clicked() {
    on_frequency_clicked(main_page->frequency->checkState());
}

void MainWindow::on_freq_v_clicked() {
    on_frequency_clicked(main_page->frequency->checkState());
}

void MainWindow::on_freq_u_clicked() {
    on_frequency_clicked(main_page->frequency->checkState());
}

void MainWindow::on_frequency_clicked(const bool checked) {
    if(time.empty()) return;

    if(!checked) {
        if(main_page->freq_a->isChecked())
            plot_time_curve(main_page->target_canvas, arma::join_rows(time, acceleration), "Target Acceleration Amplitude");
        else if(main_page->freq_v->isChecked())
            plot_time_curve(main_page->target_canvas, arma::join_rows(time, velocity), "Target Velocity Amplitude");
        else if(main_page->freq_u->isChecked())
            plot_time_curve(main_page->target_canvas, arma::join_rows(time, displacement), "Target Displacement Amplitude");
    }
    else if(main_page->freq_a->isChecked())
        plot_frequency_curve(main_page->target_canvas, f_a = perform_transform(acceleration), "Target Acceleration Amplitude");
    else if(main_page->freq_v->isChecked())
        plot_frequency_curve(main_page->target_canvas, f_v = perform_transform(velocity), "Target Velocity Amplitude");
    else if(main_page->freq_u->isChecked())
        plot_frequency_curve(main_page->target_canvas, f_u = perform_transform(displacement), "Target Displacement Amplitude");

    replot();
}

arma::mat MainWindow::perform_transform(const arma::vec& data) {
    const auto length = std::max(1024, 2 << nextpow2(time.n_elem));
    const arma::vec time_diff = arma::diff(time);
    if(time_diff.empty()) return {};
    const auto step_size = time_diff.min();
    if(step_size == 0.) return {};

    const arma::vec fft_frequency = arma::regspace(0, 1, length - 1) / (step_size * double(length));
    arma::cx_vec fft_cx_magnitude = arma::fft(data, arma::uword(length));
    arma::vec fft_magnitude = 2. * arma::abs(fft_cx_magnitude) / double(data.n_elem);

    const auto half_length = arma::uword(length) / 2;

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
    const auto x_min = x_data.min();
    const auto x_max = x_data.max();
    const auto y_max = y_data.max();
    const auto y_min = y_data.min();

    canvas->clearGraphs();

    QPen pen;
    pen.setWidth(x_data.size() > 1000 ? 1 : 2);
    pen.setColor(canvas == main_page->source_canvas ? source_color : target_color);

    canvas->addGraph();
    if(canvas == main_page->source_canvas) {
        if(x_data.size() <= 1000) canvas->graph()->setLineStyle(QCPGraph::lsImpulse);
        canvas->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 3));
    }
    canvas->graph()->setPen(pen);
    canvas->graph()->setData(to_vector(x_data), to_vector(y_data));
    canvas->graph()->rescaleAxes(true);

    canvas->xAxis->setRange(x_min, x_max);
    canvas->xAxis->scaleRange(1.02, canvas->xAxis->range().center());
    if(main_page->logarithmic->isChecked() && canvas == main_page->target_canvas) {
        canvas->yAxis->setScaleType(QCPAxis::stLogarithmic);
        QSharedPointer<QCPAxisTickerLog> log_ticker(new QCPAxisTickerLog);
        log_ticker->setTickCount(10);
        canvas->yAxis->setTicker(log_ticker);
        canvas->yAxis->setNumberFormat("eb");
        canvas->yAxis->setNumberPrecision(0);
        canvas->yAxis->setRange(std::pow(10., std::floor(std::log10(y_min))), std::pow(10., std::ceil(std::log10(y_max))));
    }
    else {
        canvas->yAxis->setScaleType(QCPAxis::stLinear);
        QSharedPointer<QCPAxisTickerFixed> linear_ticker(new QCPAxisTickerFixed);
        linear_ticker->setTickStep(std::ceil((y_max - y_min) * 2.) * .05);
        canvas->yAxis->setTicker(linear_ticker);
        canvas->yAxis->setTicks(true);
        canvas->yAxis->setTickLabels(true);
        canvas->yAxis->setNumberFormat("gb");
        canvas->yAxis->setNumberPrecision(2);
        canvas->yAxis->setRange(y_min, y_max);
        canvas->yAxis->scaleRange(1.05, canvas->yAxis->range().center());
    }

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

    QAudioFormat audio_format;
    audio_format.setSampleRate(static_cast<int>(main_page->upsample_ratio->value() / main_page->step_size->value()));
    audio_format.setChannelCount(1);
    audio_format.setSampleFormat(QAudioFormat::Float);

    arma::fvec f_data;
    if(main_page->freq_a->isChecked())
        f_data = arma::conv_to<arma::fvec>::from(acceleration);
    else if(main_page->freq_v->isChecked())
        f_data = arma::conv_to<arma::fvec>::from(velocity);
    else
        f_data = arma::conv_to<arma::fvec>::from(displacement);

    f_data /= std::max(std::fabs(f_data.min()), std::fabs(f_data.max()));

    QByteArray byte_buffer((char*)f_data.memptr(), sizeof(float) * acceleration.n_elem);

    QMediaPlayer player;
    QAudioOutput audio_output;
    player.setAudioOutput(&audio_output);
    player.setSource(QUrl::fromEncoded(byte_buffer));
    player.play();
}

void MainWindow::on_quantile_clicked() {
    arma::vec result;
    arma::vec q{.8, .85, .9, .95, .999};

    auto compute_quantile = [&](const arma::mat& input) {
        if(input.empty()) return arma::vec{};

        const arma::vec f_frequency = input.col(0);
        arma::vec f_magnitude = arma::cumsum(arma::square(input.col(1)));
        f_magnitude /= f_magnitude.back();

        result.set_size(q.n_elem);

        for(auto I = 0llu; I < q.n_elem; ++I)
            result(I) = f_frequency(arma::find(f_magnitude > q(I), 1).eval()(0));

        return result;
    };

    if(!main_page->frequency->isChecked()) {
        QMessageBox::information(this, tr("Quantile"), tr("Please toggle >Frequency< to compute."));
        return;
    }

    if(main_page->freq_a->isChecked())
        compute_quantile(f_a);
    else if(main_page->freq_v->isChecked())
        compute_quantile(f_v);
    else if(main_page->freq_u->isChecked())
        compute_quantile(f_u);

    std::ostringstream output;

    result.print(output, "Frequency @ 80% 85% 90% 95% 99.9%:");

    QMessageBox::information(this, tr("Quantile"), QString(output.str().c_str()));
}

void MainWindow::on_custom_coef_clicked() {
    filter_page.set_coefficient(custom_filter);
    filter_page.show();
}

void MainWindow::set_coefficient() {
    custom_filter = filter_page.get_coefficient();
}

void MainWindow::on_logarithmic_clicked(bool) {
    on_refresh_clicked();
}
