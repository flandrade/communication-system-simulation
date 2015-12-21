# A communication system simulation
This program simulates a basic communication system using Matlab and plots the BER curve. It includes the following components:

![Diagram](https://github.com/flandrade/communication-system-simulation/blob/master/images/diagram.png)

**main.m** is the main file.

## Components and options

### Input voices
**Voz.wav** is the audio file that contains the message. The program will plot the input signal.

#### Output

**fundamental_frequency**: the fundamental frequency of the message in the audio file. 

### Quantization
This program will plot the first input signal with the levels for quantization and the quantized signal.

#### Options
**nivel**: the number of levels for quantization.

**opcion**: quantization processes. Available options:
- 1 - Uniform
- 2 - Mu-Law
- 4 - A-Law

### Output
**errorcuantizacion**: quantization error

### Modulation
