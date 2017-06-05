# Automatic_segmentation
Program performing automatic segmentation of sound records

Program segments the loaded WAVE audio file into distinct notes by searching the difference in pitches. The segmented samples are saved into distinct audio files. It can be used for building the database of sound samples for further use in concatenative synthesis.

This program is an implementation of Piszczalski and Galler's Component Frequency Ratios algorithm. Although it works fine on most 16-bit monophonic recordings of numerous musical instruments, it still has some problems with detecting transition points between pitches in several cases (especially when the difference in frequency is small - about half a tone, or due to existence of resonant frequencies or reverb).

Program is open for further development, e.g. tuning of error tolerance, developing better continuity conditions, smoothing frequency values, providing second measure value (amplitude RMS) to improve segmentation process.
