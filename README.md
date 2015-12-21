# A communication system simulation
This program simulates a basic communication system using Matlab, and it plots BER curves in order to compare the performance of several codification algorithms. It includes the following components:

![Diagram](https://github.com/flandrade/communication-system-simulation/blob/master/images/diagram.png)

**main.m** is the main file.

## Components and options

### Input voices
**Voz.wav** is the audio file that contains the message. The program plots the input signal.

#### Output
**fundamental_frequency**: the fundamental frequency of the message in the audio file.

**x**: vector that contains the message.

### Quantization
This module plots the first input signal with the levels of quantization and the quantized signal.

#### Options
**level**: the number of levels for quantization.

**option_quantization**: quantization processes. Available options:
- 1 - Uniform
- 2 - Mu-Law
- 4 - A-Law

#### Output
**xq**: quantized message.

#### Output
**quantization_error**: quantization error.

## Codification
This simulation codifies the message using the following methods:
- Hamming (7,4)
- Convolutional codes: soft decision
- Convolutional codes: hard decision

## Modulation
This module modulates the message according the codification, and it plots the constellation for the selected modulation.

#### Options
**option_modulation**: modulation processes. Available options:
- 1 - BPSK
- 2 - QPSK
- 3 - BPSK and QPSK

#### Output
Message modulated according to the selected modulation:

**BPSK variables**: bitsm1 (no codification), bitsm2 (Hamming), bitsm3 (Convolutionl)

**QPSK variables**: bitsmqpsk1 (no codification), bitsmqpsk2 (Hamming), bitsmqpsk3 (Convolutional)

## BER Curves
This module uses a loop in order to simulate an AWGN channel with several Eb/N0 values. Eb/N0 is the energy per bit to noise power spectral density ratio. Value between 1 and 6, where 6 is for the least noisy channel.

Demodulation and decodification are performed in this loop. The program plots the BER curves of several codification algorithms.

#### Output
Probability of error gives the average rate of occurrence of decoding errors. Pe error for the codification algorithms according to the selected modulation:

**BPSK variables**: errorpe_bpsk_nocod (no codification), errorpe_bpsk_hamming (Hamming), errorpe_bpsk_hard (Convolutionl: Hard decision), errorpe_bpsk_soft (Convolutional: Soft decision)

**QPSK variables**: errorpe_qpsk_nocod (no codification), errorpe_qpsk_hamming (Hamming), errorpe_qpsk_hard (Convolutionl: Hard decision)

## Graphs

![Plots](https://github.com/flandrade/communication-system-simulation/blob/master/images/graphs.jpg)

### Acknowledgements
This program was developed during the communication course "Comunicación y codificación digital" at Universidad de las Fuerzas Armadas ESPE.
