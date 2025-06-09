function plotFrequencySpectra(time, x1f, x1e, dx1f, dx1e, rotor_speed_rads, ...
                                  flapwise_nat_freq, edgewise_nat_freq)
    % Plots frequency spectra for displacements and velocities with harmonics marked
    % and includes vertical lines for natural frequencies.
    %
    % Inputs:
    % - time: time vector
    % - x1f, x1e: displacements
    % - dx1f, dx1e: velocities
    % - rotor_speed_rpm: rotor speed in radians per second
    % - flapwise_nat_freq: flapwise natural frequency (Hz)
    % - edgewise_nat_freq: edgewise natural frequency (Hz)

    % === Parameters ===
    rotor_speed_hz = rotor_speed_rads/2/pi;  % Convert to Hz
    num_harmonics = 5;                      % Number of harmonics to plot

    % === Setup FFT parameters ===
    N = length(time);
    dt = mean(diff(time));
    Fs = 1/dt;
    f = (0:N-1)*(Fs/N);
    f_oneSided = f(1:floor(N/2));

    % === Signals to process ===
    signals = {x1f, x1e, dx1f, dx1e};
    signal_names = {'x1f (Flapwise Disp)', ...
                    'x1e (Edgewise Disp)', ...
                    'dx1f (Flapwise Vel)', ...
                    'dx1e (Edgewise Vel)'};

    % === Loop through signals ===
    for i = 1:length(signals)
        signal = detrend(signals{i});       % Remove any linear trend
        windowed_signal = signal .* hann(N);
        Y = fft(windowed_signal);
        Y_mag = abs(Y)/(N/2);
        Y_mag_oneSided = Y_mag(1:floor(N/2));

        % Determine axis limits
        y_max = max(Y_mag_oneSided);
        max_plot_freq = max([num_harmonics+1, ...
                             flapwise_nat_freq*1.2, ...
                             edgewise_nat_freq*1.2]) * rotor_speed_hz;
        max_plot_freq = min(max_plot_freq, max(f_oneSided));  % don't exceed Nyquist

        % === Plot ===
        figure;
        plot(f_oneSided, Y_mag_oneSided, 'b', 'LineWidth', 1.5);
        xlabel('Frequency [Hz]');
        ylabel('Amplitude');
        title(['Frequency Spectrum of ', signal_names{i}]);
        grid on;
        hold on;

        % Plot harmonics
        for k = 1:num_harmonics
            x_harmonic = k * rotor_speed_hz;
            if x_harmonic <= max_plot_freq
                xline(x_harmonic, 'k--', sprintf('%dP', k), ...
                    'LabelOrientation', 'horizontal', ...
                    'LabelVerticalAlignment', 'bottom', ...
                    'LabelHorizontalAlignment', 'center', ...
                    'FontSize', 10);
            end
        end

        % Plot natural frequencies
        xline(flapwise_nat_freq, 'r-.', 'LineWidth', 1.5);
        text(flapwise_nat_freq, y_max * 0.9, 'Flapwise NF', ...
            'Rotation', 90, ...
            'Color', 'r', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 10);

        xline(edgewise_nat_freq, 'm-.', 'LineWidth', 1.5);
        text(edgewise_nat_freq, y_max * 0.9, 'Edgewise NF', ...
            'Rotation', 90, ...
            'Color', 'm', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 10);

        % Set axis limits
        xlim([0, max_plot_freq * 1.2]);
        ylim([0, y_max * 1.1]);

        legend('Spectrum');
    end
end
