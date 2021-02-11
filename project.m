
 % Summary of this function and detailed explanation goes here

 % First convert the record into matlab (creates recordm.mat):
 % wfdb2mat -r record


load('t100m.mat');
signal = val(1,:);

fs = 200;
signal = signal(:);
signal_mean = mean(signal);
signal = signal - signal_mean

delay = 0;
skip = 0;
m_selected_RR = 0;
mean_RR = 0;
ser_back = 0;
ax = zeros(1,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lowpass filter
%z = (1-z^-6)^2 / (1-z^-1)^2
%zakasnitev filtra = 2 * (6-1)/2 = 5

a_1 = [1,-1];
b_1 = [1,0,0,0,0,0,-1];

a = conv(a_1,a_1);
b = conv(b_1,b_1);
signal_f_l = filter(b,a,signal);

signal_f_l = signal_f_l / max(abs(signal_f_l)); %normalize

delay  = delay + 5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%highpass filter
%z = (-1+32z^-16 + z^-32)/(1+z^-1)
%zakasnitev filtra = 1 * (32-1)/2
a_2 = [1,1];
b_2 = zeros(1,33);
b_2(1) = -1;
b_2(17) = 32;
b_2(33) = 1;

signal_f_h = filter(b_2,a_2,signal_f_l);
signal_f_h = signal_f_h / max(abs(signal_f_h)); %normalize

delay = delay + round((32-1)/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%derivative filter
%z = (1/8T)(-z^-2 - 2z^-1 + 2z + z^2) 

const = (1/8) *fs;

a_3 = [1];
b_3 = const.*[1,2,0,-2,-1];

signal_f_d = filter(b_3,a_3,signal_f_h);
signal_f_d = signal_f_d / max(abs(signal_f_d));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%squaring function
%Y(nT) = x(nT)^2
signal_s = signal_f_d.^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%moving window integration
%Y(nT) = (1/N)[x(nT-(N-1)T)+x(nT-(N-2)T)+ ... +x(nT)]

N = 30; %as in paper

a_4 = [1];
b_4 = ones(1 ,N)/N;

signal_m_w_i = filter(b_4, a_4, signal_s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fiducial mark
miliseconds = 200; %QRS complex can not occur more closely than this physiologically
dist = round(miliseconds / 1000 * fs) ;

[pks,locs] = findpeaks(signal_m_w_i(1 : 2 * fs),'MINPEAKDISTANCE',dist); %value and location

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjusting the thresholds
seconds = 2; %2 seconds into the signal to determine the threshold
signal_2s = signal_m_w_i(1 : seconds * fs);

SPK = max(signal_2s);
NPK = 0;
T1 = 0.75 * NPK + 0.25 * SPK;
T2 = 0.5 * T1;

signal_f_h_2s = signal_f_h(1 : seconds * fs);

SPK_s = max(signal_f_h_2s);
NPK_s = 0;
T1_s = NPK_s + 0.25 * (SPK_s - NPK_s);
T2_s = 0.5 * T1_s;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%detecting QRS

% initialization
qrs_value = []; %value of qrs
qrs_position = []; %position
qrs_candidate = []; %candidates
qrs_raw = []; %final positions
T_check = 0;
THR_s = [];
THR2_s = [];

THR_1 = [];
THR_2 = [];
for i = 1 : length(pks)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %peaks in the filtered signal
    if (locs(i) - N >= 1)
        [y_i, x_i] = max(signal_f_h(locs(i) - N : locs(i))); %high_pass filter signal
    else
        if (i==1)
            [y_i,x_i]=max(signal_f_h(1:locs(i))); %high_pass filter signal
        end
    end
    x_i = x_i + locs(i) - N - 1;
    qrs_candidate = [qrs_candidate, x_i];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % adjusting the average RR interval and rate limits
    if (length(qrs_position) >= 9) %eight most recent beats
        RR_int = diff(qrs_position)
        mean_RR1 = mean(RR_int(end-7:end)); %average of the last eight beats having RR intervals
        %in RR_int2 we put beats that falls between low and high
        %RR-interval limits
        RR_int2 = RR_int(RR_int >= 0.92 * mean_RR1);
        RR_int2 = RR_int2(RR_int2 <= 1.16 * mean_RR1);
        if (length(RR_int2) >= 8) 
            mean_RR2 = mean(RR_int2(end-7:end));
        else
            mean_RR2 = mean(RR_int2);
        end       
    else
        mean_RR2 = mean(diff(qrs_position));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %qrs complex not found during the interval specified by RR MISSED LIMIT
    %threshold - there were no peaks in the time = 166% average time
    if ((length(qrs_position) ~= 0) && ((locs(i) - qrs_position(end)) >= round(1.66 * mean_RR2))) 

        %go back looking for peaks - max in that interval
        [pks_missed, locs_missed] = max(signal_m_w_i(qrs_position(end) + dist : locs(i) - dist));
        locs_missed = qrs_position(end)+ dist + locs_missed -1; 

        % check if it above the second threshold
        if (pks_missed > T2)

           % find the peak in the filtered signal
           [pks_filtered, index_filtered] = max(signal_f_h(locs_missed - N : locs_missed));

           % check if it above the second threshold -filtered signal
           if (pks_filtered > T2_s)
               qrs_value = [qrs_value, pks_missed];
               qrs_position = [qrs_position, locs_missed];
               qrs_raw=[qrs_raw, index_filtered + locs_missed - N - 1];
               % change the threshold
               SPK = 0.25 * pks_missed + 0.75 * SPK;
               SPK_s = 0.25 * pks_filtered + 0.75 * SPK;
               T1 = 0.75 * SPK + 0.25 * NPK;
               T2 = 0.5 * T1;
               T1_s = 0.75 * SPK_s + 0.25 * NPK_s;
               T2_s = 0.5 * T2;
           end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    upper_bound_seconds = 360 %miliseconds
    upper_bound = round(upper_bound_seconds / 1000 * fs)
    upper_bound_seconds_t = 75
    upper_bound_t = round(upper_bound_seconds_t*fs)

    % distinguish noise and QRS
    if (pks(i) > T1) 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %T-wave identification

        %condition for T-wave identification
        if ((length(qrs_position) ~= 0) && (locs(i) - qrs_position(end) <= upper_bound))
            % need to check if it is a T-wave

            signal_t = signal_f_h(qrs_raw(end) - upper_bound_t : qrs_raw(end));
            Slope_qrs = max(abs(diff(signal_t)));
            Slope = max(abs(diff(signal_f_h(x_i - upper_bound_t : x_i))));

            %if less than half -> T wave
            if (Slope <= 0.5*Slope_qrs)
                T_check=1;
            end
        end

        %~T wave -> QRS complex
        if (T_check==0)

            if (y_i > T1_s) 
                %it is a QRS
                qrs_position = [qrs_position, locs(i)];
                qrs_value = [qrs_value, pks(i)];
                qrs_raw = [qrs_raw, x_i];
                SPK = 0.125 * pks(i) + 0.875 * SPK;
                SPK_s = 0.125 *y_i + 0.875 * SPK_s;
            else
                NPK = 0.125 * pks(i) + 0.875 * NPK;
                NPK_s = 0.125 * y_i + 0.875 * NPK_s;
            end
        else
            % This is a T wave
            NPK = 0.125 * pks(i) + 0.875 * NPK;
            NPK_s = 0.125 * y_i + 0.875 * NPK_s;   
        end
    else
        % This is noise
        NPK = 0.125 * pks(i) + 0.875 * NPK;
        NPK_s = 0.125 * y_i + 0.875 * NPK_s;
    end
    % update the threshold
    T1 = 0.75 * NPK + 0.25 * SPK;
    T2 = 0.5 * T1;
    T1_s = 0.75 * NPK_s + 0.25 * SPK_s;
    T2_s = 0.5 * T1_s;
    THR_s = [THR_s, T1_s]
    THR2_s = [THR2_s, T2_s];

    THR_1 = [THR_1, T1]
    THR_2 = [THR_2, T2]
    % reset T_check
    T_check=0;



    % show the results
    %figure();
    %plot((1:length(signal)),signal /  max(abs(signal)));
    %for j=1:length(qrs_raw)
    %    hold on; xline(qrs_raw(j)-delay, '--r')
    %end
end


    

    
    
    
    
    
    

    figure;
    ax(1) = subplot(321);
    plot((1:length(signal)),signal /  max(abs(signal)));
    for j=1:length(qrs_raw)
        hold on; xline((qrs_raw(j)-delay), '--r')
    end
    axis tight;
    title('Raw signal');

    ax(2)=subplot(322);
    plot(signal_f_l);
    axis tight;
    title('Low pass filtered');

    ax(3)=subplot(323);
    plot(signal_f_h);
    axis tight;
    title('High pass filtered');

    ax(4)=subplot(324);
    plot(signal_f_d );
    axis tight;
    title('Derivative filtered')

    ax(5)=subplot(325);
    plot(signal_s );
    axis tight;
    title('Square function')

    ax(6)=subplot(326);
    plot(signal_m_w_i / max(signal_m_w_i));
    %hold on; plot(qrs_candidate, THR_1)
    %hold on; plot(qrs_candidate, THR_2)
    for j=1:length(qrs_raw)
        hold on; xline((qrs_raw(j) + 11), '--r')
        % all delay - delay low pass = 11
    end
    axis tight;
    title('Moving Window integration')


    %qrs_raw - delay
