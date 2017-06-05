#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>


// Allegro5 headers
#include <allegro5/allegro.h>
#include <allegro5/allegro_primitives.h>

//FFTW3 Library
#include <fftw3.h>

using namespace std;

struct wavFile
{
    //WAVE file header descriptions
        char chunkID[5] = {'R','I','F','F'}; //"RIFF"
        unsigned int chunkSize; //size of the file
        char format[5] = {'W','A','V','E'}; //"WAVE"

        char subchunk1ID[5] = {'f','m','t',' '}; //"fmt "
        unsigned int subchunk1Size;
        short audioFormat; //PCM = 1 (i.e. linear quantization)
        short channels; //channel numbers: 1-MONO, 2-STEREO
        unsigned int sampleRate;
        unsigned int byteRate; // == SampleRate * NumChannels * BitsPerSample/8
        short blockAlign; //number of bytes for one sample including all channels
        short bitsPerSample;

        char subchunk2ID[4] = {'d','a','t','a'}; // "data"
        unsigned int dataSize; // == NumSamples * NumChannels * BitsPerSample/8
        short* data; // actual sound data
        unsigned int nsam;
        void loadData()
        {
//            nsam = 8*(dataSize/channels)/bitsPerSample;
//            nsam = dataSize/blockAlign;
            nsam = dataSize/blockAlign;
            data = new short[nsam];
        }
}wavefil;

struct harmonic
{
    double frequency;
    double f0;          //suggested fundamental
    double amplitude;

    harmonic(double a, double b, double c) //structure constructor
    {
        frequency = a;
        f0 = b;
        amplitude = c;
    }

    inline bool operator!=(const harmonic &harm)  //finding if compared frequency differs by more than 10%
    {
        return (frequency < 0.9*(harm.frequency) || frequency > 1.1*(harm.frequency));
    }
};

struct higher_amplitude //structure for sorting harmonic structures due to their FFT amplitudes
{
    inline bool operator() (const harmonic &harm1, const harmonic harm2)
    {
        return (harm1.amplitude > harm2.amplitude);
    }
};

struct segment_struct
{
    unsigned start; //the number of signal sample where the segment begins
    unsigned finish; //the number of signal sample where the segment ends
    segment_struct(unsigned a, unsigned b) //segment constructor
    {
        start = a;
        finish = b;
    }
};

//functions
bool loadWavFile(wavFile* wavef, string filePath) //loading .wav file data
{
    ifstream fw(filePath.c_str());
    if(!fw.good())
    {
        cout << "File couldn't be opened." << endl;
        return false;
    }

    //Searching for 'RIFF' chunk
    char chunkId[4] = {0};
    if(fw.read(chunkId, 4) &&
         (chunkId[0] == 'R' &&
          chunkId[1] == 'I' &&
          chunkId[2] == 'F' &&
          chunkId[3] == 'F'))
          			 cout<<wavef->chunkID<<endl;

    fw.read((char*)&wavef->chunkSize, sizeof(unsigned long));

    //Searching for 'WAVE' chunk
    int chunkSize = 0;
    while(fw.read(chunkId, 4) &&
         (chunkId[0] != 'W' &&
          chunkId[1] != 'A' &&
          chunkId[2] != 'V' &&
          chunkId[3] != 'E'))
    {
      fw.read(reinterpret_cast<char*>(&chunkSize), 4); // Read the chunk's size
      fw.seekg(chunkSize, ios_base::cur); // Skip the chunk
    }
    cout<<wavef->format<<endl;

    //Searching for 'fmt ' chunk
    while (fw.read(chunkId, 4) &&
         (chunkId[0] != 'f' ||
          chunkId[1] != 'm' ||
          chunkId[2] != 't' ||
          chunkId[3] != ' '))
    {
      fw.read(reinterpret_cast<char*>(&chunkSize), 4); // Read the chunk's size
      fw.seekg(chunkSize, ios_base::cur); // Skip the chunk
    }

    //Reading 'fmt ' components and saving them to wavFile structure
    fw.read((char*)&wavef->subchunk1Size, sizeof(unsigned int));
    fw.read((char*)&wavef->audioFormat, sizeof(short));
    fw.read((char*)&wavef->channels, sizeof(short));
    fw.read((char*)&wavef->sampleRate, sizeof(unsigned));
    fw.read((char*)&wavef->byteRate, sizeof(unsigned));
    fw.read((char*)&wavef->blockAlign, sizeof(short));
    fw.read((char*)&wavef->bitsPerSample, sizeof(short));

    //Searching for 'data' chunk
    while (fw.read(chunkId, 4) &&
         (chunkId[0] != 'd' ||
          chunkId[1] != 'a' ||
          chunkId[2] != 't' ||
          chunkId[3] != 'a'))
    {
      fw.read(reinterpret_cast<char*>(&chunkSize), 4); // Read the chunk's size
      fw.seekg(chunkSize, ios_base::cur); // Skip the chunk
    }

    fw.read((char*)&wavef->dataSize, sizeof(unsigned int));
    wavef->loadData(); //allocating memory to dynamic array for data reading

    short db[wavef->nsam]; //data buffer
    fw.read((char*)&db,wavef->dataSize);
    unsigned int j = 0;

    //assigning the WAVE data from read data buffer
    while(j<(wavef->nsam))
    {
        wavef->data[j] = db[j];
        j++;
    }
    fw.close();

    //Writing .wav file properties to the console
    if(wavef->channels == 1 ? cout << "MONO\n" : cout << "STEREO\n");
    cout << "Sample rate: " << wavef->sampleRate << "Hz" << endl;
    cout << wavef->byteRate << " B/s" << endl;
    cout << "Block alignment: " << wavef->blockAlign << endl;
    cout << wavef->bitsPerSample << "bit sound" << endl;
    cout << "data chunk size: " << wavef->dataSize << endl;
    cout << "number of samples: "<< wavef->nsam << endl;

    return true;
}

bool writeWavSegments(wavFile* wavef, vector <segment_struct> segments) //writing the output segments to .wav files
{
    cout<<endl;
    for(unsigned int seg=0; seg<segments.size(); seg++) //for each found segment
    {
        cout<<"Writing "<<seg+1<<" segment of .wav file... ";
        string name;
        stringstream ss;
        ss<<"out"<<seg+1<<".wav"; //constructing name of the output file with the number of consecutive segment
        ofstream wavSegment; //output .wav file
        wavSegment.open(ss.str().c_str());

        stringstream ss2;
        ss2<<"out"<<seg+1<<".txt";
        ofstream segmentval; //text file with segment values
        segmentval.open(ss2.str().c_str());

        //writing .wav file header
        char chunkID[4] = {'R','I','F','F'};
        wavSegment.write(chunkID, 4);
        wavSegment.write((char*)&wavef->chunkSize, sizeof(unsigned int));
        char format[4] = {'W','A','V','E'};
        wavSegment.write(format, 4);
        char subchunk1ID[4] = {'f','m','t',' '};
        wavSegment.write(subchunk1ID, 4);
        wavSegment.write((char*)&wavef->subchunk1Size, sizeof(unsigned int));
        wavSegment.write((char*)&wavef->audioFormat, sizeof(short));
        wavSegment.write((char*)&wavef->channels, sizeof(short));
        wavSegment.write((char*)&wavef->sampleRate, sizeof(unsigned int));
        wavSegment.write((char*)&wavef->byteRate, sizeof(unsigned int));
        wavSegment.write((char*)&wavef->blockAlign, sizeof(short));
        wavSegment.write((char*)&wavef->bitsPerSample, sizeof(short));
        char subchunk2ID[4] = {'d','a','t','a'};
        wavSegment.write(subchunk2ID, 4);
        unsigned int numOfSamples;
        numOfSamples = segments[seg].finish - segments[seg].start;
        unsigned long dataSize = numOfSamples*wavef->blockAlign;
        wavSegment.write((char*)&dataSize, sizeof(unsigned int));

        //writing values to .wav file sample by sample
        for(unsigned i=segments[seg].start; i<=segments[seg].finish; i++)
        {
            wavSegment.write((char*)&wavef->data[i], wavef->blockAlign);
            wavSegment.write((char*)&wavef->data[i], wavef->blockAlign);
            segmentval<<wavef->data[i]<<endl;
        }

        wavSegment.close();
        segmentval.close();
        cout<<"done."<<endl;
    }
    return true;
}

vector<segment_struct> analyser(short data[], unsigned long numOfSamples)
{
        cout<<"Analyser started..."<<endl;
        unsigned int N = floor((wavefil.sampleRate*32/1000)/2); //number of FFT bins (1/2 samples every 32ms)
        cout<<"Number of samples per interval : "<< N<<endl;
        double *in, *out;
        in = (double*) fftw_malloc(sizeof(double)*N);
        out = (double*) fftw_malloc(sizeof(double)*N);
        fftw_plan p;

        //FFT computing
        ofstream dataFile_FFT;
        dataFile_FFT.open("dataFile_FFT.txt");
        vector < vector <double> > samplesArray;
        int index = 0;
        while(N*index<numOfSamples)
        {
            //Computing FFT
            p = fftw_plan_r2r_1d(N, in, out, FFTW_REDFT00, FFTW_ESTIMATE);

            for(unsigned int i=0; i<N; i++)
            {
                if((i+N*index)<numOfSamples)
                    in[i] = (double) data[i+(N*index)];
            }
            fftw_execute(p);

            vector <double> tempArray;
            for(unsigned int i=0; i<N; i++)
            {
                tempArray.push_back(abs(out[i])); //saving absolute value of output in container
                dataFile_FFT<<tempArray[i]<<endl;
            }
            samplesArray.push_back(tempArray);
            dataFile_FFT<<endl;

            fftw_destroy_plan(p);
            index++;
        }
        cout<<index-1<<" FFTs computed. "<<endl;
        dataFile_FFT.close();

        fftw_free(in);
        fftw_free(out);

        //Find 6 greatest peaks in FFT arrays
        double peaks[samplesArray.size()][6];       // array storing peaks amplitude
        double frequencies[samplesArray.size()][6]; //array storing value of peaks frequencies
        ofstream peakFil;
        peakFil.open("peakFil.txt");
        for(unsigned int i=0; i<samplesArray.size(); i++)
        {
            for(int k=0; k<6; k++) //initializing peak values as zeros
            {
                peaks[i][k] = 0;//samplesArray[i][0];
            }
            for(unsigned int j=0; j<samplesArray[i].size(); j++)
            {
                for(int k=0; k<6; k++)
                {
                    if(samplesArray[i][j]>peaks[i][k])      //if value is greater than k's peak
                    {
                        for(int l=k+1; l<6; l++)        //move values
                        {
                            if(peaks[i][l] == peaks[i][l-1])
                                peaks[i][l] = peaks[i][l-1];
                            if(frequencies[i][l] == frequencies[i][l-1])
                                frequencies[i][l] = frequencies[i][l-1];
                        }
                        peaks[i][k] = samplesArray[i][j];               //assign peak amplitude
                        frequencies[i][k] = (j+1)*wavefil.sampleRate/N; //assign value of peak's frequency
                        break;
                    }
                }
            }
            for(int k=0; k<6; k++)
            {
                if(frequencies[i][k] != frequencies[i][k]) //checking if peak value is NaN
                    frequencies[i][k] = 0;
                peakFil<<frequencies[i][k]<<" "; //write peaks to text file
            }
            peakFil<<endl;
        }
        peakFil.close();
        cout<<"Peaks evaluated."<<endl;

        //Calculate peak ratios
        ofstream ratioFil;
        ratioFil.open("ratioFil.txt");
        double ratios[samplesArray.size()][6][6]; //6x6 matrix of ratios
        for(unsigned int i=0; i<samplesArray.size(); i++)
        {
            for(int j=0; j<6; j++)
            {
                for(int k=0; k<6; k++)
                {
                    if(frequencies[i][k]!=0)
                        ratios[i][j][k] = frequencies[i][j]/frequencies[i][k];
                    else
                        ratios[i][j][k] = 0;
                    ratioFil<<ratios[i][j][k]<<" ";
                }
                ratioFil<<endl;
            }
            ratioFil<<endl;
        }
        ratioFil.close();
        cout<<"Peak ratios calculated."<<endl;

        //estimate greatest suggestion of f0 by ocurrence of harmonics
        ofstream f0File;
        f0File.open("f0.txt");
        vector <vector <harmonic>> harmonics(samplesArray.size());
        for(unsigned int i=0; i<samplesArray.size(); i++)
        {
            for(int j=0; j<6; j++)
            {
                for(int k=0; k<6; k++)
                {
                    if(((ratios[i][j][k] > 0.9*(1.0/2.0)) && (ratios[i][j][k] < 1.1*(1.0/2.0))) ||
                       ((ratios[i][j][k] > 0.9*(1.0/3.0)) && (ratios[i][j][k] < 1.1*(1.0/3.0))) ||
                       ((ratios[i][j][k] > 0.9*(1.0/4.0)) && (ratios[i][j][k] < 1.1*(1.0/4.0))) ||
                       ((ratios[i][j][k] > 0.9*(1.0/5.0)) && (ratios[i][j][k] < 1.1*(1.0/5.0))) ||
                       ((ratios[i][j][k] > 0.9*(1.0/6.0)) && (ratios[i][j][k] < 1.1*(1.0/6.0))))
                    {
                        harmonics[i].push_back(harmonic(frequencies[i][j], frequencies[i][j], peaks[i][j]));
                    }

                    else if(((ratios[i][j][k] > 0.9*(2.0/3.0)) && (ratios[i][j][k] < 1.1*(2.0/3.0))) ||
                            ((ratios[i][j][k] > 0.9*(2.0/5.0)) && (ratios[i][j][k] < 1.1*(2.0/5.0))))
                    {
                        harmonics[i].push_back(harmonic(frequencies[i][j], frequencies[i][j]/2.0, peaks[i][j]));
                    }

                    else if(((ratios[i][j][k] > 0.9*(3.0/4.0)) && (ratios[i][j][k] < 1.1*(3.0/4.0))) ||
                            ((ratios[i][j][k] > 0.9*(3.0/5.0)) && (ratios[i][j][k] < 1.1*(3.0/5.0))))
                    {
                        harmonics[i].push_back(harmonic(frequencies[i][j], frequencies[i][j]/3.0, peaks[i][j]));
                    }

                    else if((ratios[i][j][k] > 0.9*(4.0/5.0)) && (ratios[i][j][k] < 1.1*(4.0/5.0)))
                    {
                        harmonics[i].push_back(harmonic(frequencies[i][j], frequencies[i][j]/4.0, peaks[i][j]));
                    }

                    else if((ratios[i][j][k] > 0.9*(5.0/6.0)) && (ratios[i][j][k] < 1.1*(5.0/6.0)))
                    {
                        harmonics[i].push_back(harmonic(frequencies[i][j], frequencies[i][j]/5.0, peaks[i][j]));
                    }
                }
            }
            if(harmonics[i].size() == 0)
                harmonics[i].push_back(harmonic(0,0,0));
            for(unsigned int j=0; j<harmonics[i].size(); j++)
            {
                f0File<<harmonics[i][j].frequency<<" ";
            }
            f0File<<endl;
        }
        f0File.close();
        cout<<"f0 extracted."<<endl;

        //estimation part 2
        //evaluate histogram
        vector <vector <harmonic>> hist(samplesArray.size());
        ofstream f0_hist;
        f0_hist.open("f0_hist.txt");
        for(unsigned int i=0; i<samplesArray.size(); i++)
        {
            hist[i].push_back(harmonics[i][0]);              //start with adding first frequency
            for(unsigned int j=1; j<harmonics[i].size(); j++)
            {
                for(unsigned int k=0; k<hist[i].size(); k++) //look for the same f0 and sum its amplitude
                {
                    if(harmonics[i][j].frequency > 0.9*(hist[i][k].frequency) && harmonics[i][j].frequency < 1.1*(hist[i][k].frequency))
                    {
                        hist[i][k].amplitude += peaks[i][j];
                        break;
                    }
                }
            }
            sort(hist[i].begin(), hist[i].end(), higher_amplitude());
            f0_hist<<hist[i][0].frequency<<'\t'<<hist[i][0].amplitude<<endl;
        }
        cout<<"Histogram evaluated."<<endl;
        f0_hist.close();

        //search for f0 pattern
        cout<<"Segmentation... ";
        ofstream f0_out;
        vector <segment_struct> segments;
        int count_segment = 4;
        int num_segment = 0;
        segment_struct temp_seg(0,0);
        f0_out.open("output.txt");
        f0_out<<hist[0][0].frequency<<endl<<hist[1][0].frequency<<endl<<hist[2][0].frequency<<endl<<hist[3][0].frequency<<endl;
        for(unsigned int i=3; i<hist.size()-4; i++)
        {
            //if the loop has reached the last iteration, end the last segment
            if(i==hist.size()-5)
            {
                temp_seg.finish = i*N-1;
                if(count_segment>10)
                {
                    cout<<num_segment+1<<": "<<temp_seg.start<<"-"<<temp_seg.finish<<endl;
                    segments.push_back(temp_seg);
                    num_segment++;
                }
            }
            //check if the checked frame falls out of the pattern, i.e. if frequency is in tolerance range with 3 previous and 3 next frames
            //meaning there is another segment
            if((hist[i][0]!=hist[i-1][0]) &&
               (hist[i][0]!=hist[i-2][0]) &&
               (hist[i][0]!=hist[i-3][0]) &&
               (hist[i-1][0]!=hist[i+1][0]) &&
               (hist[i-1][0]!=hist[i+2][0]) &&
               (hist[i-1][0]!=hist[i+3][0]))
            {
                temp_seg.finish = i*N-1;
                f0_out<<endl;

                //if segment has lasts for more than 10 frames, count it as actual segment
                if(count_segment>10)
                {
                    cout<<num_segment+1<<": "<<temp_seg.start<<"-"<<temp_seg.finish<<endl; //write first and last signal sample of the segment
                    segments.push_back(temp_seg);
                    num_segment++;
                }

                //segment has been found, so clear the counter for the next segment
                temp_seg.start = i*N;
                count_segment = 0;
            }
            f0_out<<hist[i][0].frequency<<endl;
            count_segment++;
        }
        f0_out.close();
        cout<<"done."<<endl;

        return segments;
}

void drawWaveForm(short data[], unsigned long numOfSamples)
{
    //Allegro initialization
    ALLEGRO_DISPLAY *display = NULL;
    int wWidth = 1280, wHeight = 680;
    al_init();
    al_init_primitives_addon();
    al_install_keyboard();
    ALLEGRO_KEYBOARD_STATE keyState;
    display = al_create_display(wWidth, wHeight);
    ALLEGRO_COLOR color_black = al_map_rgb(255, 255, 255);

    //Find the biggest absolute values for window resolution purposes
    int wndresW = numOfSamples/wWidth;
    double wndresH;
    short aMax=abs(data[0]);

    //write sound sample's values to txt file
    ofstream ampsFil;
    ampsFil.open("vals.txt");
    for(unsigned int i=0; i<numOfSamples; i++)
    {
        ampsFil<<data[i]<<endl;
    }
    ampsFil.close();

    //find the maximum value for calculation of the window resolution
    for(unsigned int i=1; i<numOfSamples; i++)
    {
        if(abs(data[i])>aMax)
            aMax = abs(data[i]);
    }
    wndresH = aMax/(wHeight*2);

    cout<<endl<<"UP/DOWN to switch resolution"<<endl
    <<"ESC to quit and perform segmentation."<<endl;

    int mpctor = 1;
    bool windowed = true;
    do
    {
        //Plotting fitted to window
        al_clear_to_color(al_map_rgb(0,0,0));
        if(windowed)
        {
            int sum1 = 0, sum2 = 0;
            for(int i=1; i<wWidth; i++)
            {
                al_put_pixel(i, wHeight/2, color_black); //Straight line of 0 value
                //algepraic sum of values to be plotted
                sum1 = 0;
                sum2=0;
                for(int j=0; j<wndresW; j++)
                {
                    sum1 += data[(i-1)*wndresW+j];
                    sum2 += data[i*wndresW+j];
                }
                sum1=sum1/wndresW;
                sum2=sum2/wndresW;
                al_draw_line(i-1, sum1/wndresH + wHeight/2, i, sum2/wndresH + wHeight/2, color_black, 1.0);
            }
            al_flip_display();
        }
        else
        {
            for(int i=1; i<wWidth; i++)
            {
                al_draw_line(i-1, wHeight/2, i, wHeight/2, color_black, 1.0);

                if((numOfSamples - mpctor*wWidth)>(unsigned int)i)
                    al_draw_line(i-1, data[i-1+mpctor*wWidth] + wHeight/2, i, data[i+mpctor*wWidth] + wHeight/2, color_black, 1.0);
                else
                    al_put_pixel(i, wHeight/2, color_black);
            }
            al_flip_display();
        }

        al_get_keyboard_state(&keyState);
        if(al_key_down(&keyState, ALLEGRO_KEY_RIGHT))
        {
            if(mpctor*wWidth<(int)numOfSamples)
                mpctor++;
            else
                mpctor=mpctor;
        }

        else if(al_key_down(&keyState, ALLEGRO_KEY_LEFT))
        {
            if(mpctor>1)
                mpctor--;
        }
        else if(al_key_down(&keyState, ALLEGRO_KEY_UP))
        {
            windowed = true;
        }
        else if(al_key_down(&keyState, ALLEGRO_KEY_DOWN))
        {
            windowed = false;
        }
    } while(!al_key_down(&keyState, ALLEGRO_KEY_ESCAPE));

    al_destroy_display(display);
}


int main()
{
    cout << "Enter the file path: ";
    string filePath;
    cin >> filePath;

    //loading .wav file sound sample
    if(loadWavFile(&wavefil, filePath))
    {
        cout << "File loaded successfully." << endl;
    }
    else
    {
        cout << "Problem occurred while loading the file." << endl;
        return 1;
    }

    //plotting waveform
    //not essential for program to work; can be commented out
    drawWaveForm(wavefil.data, wavefil.nsam);

    //Analyse and segment data
    vector <segment_struct> segments = analyser(wavefil.data, wavefil.nsam);

    //write segment values to .txt and .wav files
    writeWavSegments(&wavefil, segments);

    return 0;
}
