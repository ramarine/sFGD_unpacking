/* NOTE ON COMPILING THIS CODE:::::
     This .C file has ROOT plotting written into it, and has been written on a system where ROOT has been configured with
     C++ 14. In order to compile this on your machine, you need to run the following command (without quotations):
     " g++ -o data_decode_v2.exe data_decode_v2.c `root-config --cflags --libs` "

     In order to run this on your machine, you need to run the executable with the binary file you want to data check
     as a command line input and the output path with the file name [and output path]:
     " ./data_decode_v2.exe "${input_file}" ["${output_path}"]

     Depending on the machine, you may not have the rights to make the output folder and/or the files as this script executes, so you may need
     to FIRST create the folder and files to match the name of your .bin file, i.e.:
      Input: path/binaryFile.bin                         or             Input: path/binaryFile.bin /decoder_result_path/

      Output: path/binaryFile_LogErrors.txt                             Output: /decoder_result_path/binaryFile_LogErrors.txt
              path/binaryFile_LogWords.txt                                      /decoder_result_path/binaryFile_LogWords.txt
              path/binaryFile_histo.pdf                                         /decoder_result_path/binaryFile_histo.pdf
              path/binaryFile_histo.root                                        /decoder_result_path/binaryFile_histo.root

*/

#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <sys/ioctl.h>
#include <time.h>
#include <vector>
#include <array>
#include <fstream>

//ROOT Libraries
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TH3D.h>
#include <TFile.h>
#include <TTree.h>
#include <TLine.h>
#include <TPaveText.h>


const int N_FEB = 14;     // Number of FEBs per crate, Note that this is also used for the number of FEB error in the OCB data packet trailer
int nSlot_FEB = 0;        // Slot of the active FEB
int buflen = 4;
bool VERBOSE = false;   // Print out the Event data
bool READ_only = false; // Print out only *_LogErrors.txt


//words to read from file
uint32_t word;            // Word to read
uint32_t prev_word;       // Previous word
FILE *read_fp, *write_fp_words, *write_fp_errors; // Files to read the word from and to write the outputs
int n_line = 0;           // Line number in write_fp_words

bool is_checks_done = 0; // To keep track of the checks state w.r.t.current event, 0 = not done, 1 = done
int Nempty_words = 0;    // Total number of empty word
int are_prev_word_empty = 0; // Keep track of the number of empty words in a row
std::vector<std::array<int, 3>> where_empty_words; // Vector containing the location of the empty words (event, line, number in a row)
int Nevent = 0;          // Total number of event == nbre of OCB header
int Nevent_checked = 0;  // Total number of event checked == nbre of full event (no empty word)
int Nerrors = 0;         // Total number of errors
int Nskip_event = 0;     // Number of times there are a jump in the event number series
std::vector<std::array<int, 3>> where_skip_events;

//OCB header
std::vector<int> where_missing_OCB_header; // Event number with a missing OCB data packet header

//Gate header
bool prev_gate_word = 0; // Keep track of the previous gate (header = 0, trailer = 1), for the gate time

int N_missing_gate_headerA = 0; // Total number of gate header A (id0)
std::array<std::vector<std::vector<int>>, N_FEB> where_missing_gate_headerA; // Precise location of the missing gate header A: Event number, n_line (in *_LogWords)
int N_missing_gate_headerB = 0; // Missing Gate header B (id 0)
std::array<std::vector<std::vector<int>>, N_FEB> where_missing_gate_headerB;
int N_missing_gate_header_time = 0; // Missing Gate time (for header, id 7)
std::array<std::vector<std::vector<int>>, N_FEB> where_missing_gate_header_time;

int N_MmFEB_gate_number = 0; // Mismatch of the Gate number between header and trailer
std::array<std::vector<std::vector<int>>, N_FEB> where_MmFEB_gate_number;

int N_Mm_gate_type = 0; // Mismatch of the Gate type (header A) betweem FEBs
std::array<std::vector<std::vector<int>>, N_FEB> where_Mm_gate_type;
int N_Mm_gate_number = 0; // Mismatch of the Gate number (header A) betweem FEBs
std::array<std::vector<std::vector<int>>, N_FEB> where_Mm_gate_number;
int N_Mm_gate_time = 0; // Mismatch of the Gate time from GTS (header B) betweem FEBs
std::array<std::vector<std::vector<int>>, N_FEB> where_Mm_gate_time;
int N_Mm_gate_timeOpenGate = {0};
std::array<std::vector<std::vector<int>>, N_FEB> where_Mm_gate_timeOpenGate;

//GTS header
int N_missing_GTS_header = 0; // Missing GTS header (ID 1)
std::array<std::vector<std::vector<int>>, N_FEB> where_missing_GTS_header;

//Hit words
int Gain = 0;
int Channel_id = 0;

//GTS trailers
std::array<int, N_FEB> N_GTS_trailer1_checked = {0}; // Number of GTS trailer1 read during checks per FEB
int N_missing_GTS_trailer1 = 0; //Missing GTS trailer 1 (ID 4)
std::array<std::vector<std::vector<int>>, N_FEB> where_missing_GTS_trailer1;
int N_missing_GTS_trailer2 = 0; // MIssing GTS trailer 2 (ID 5)
std::array<std::vector<std::vector<int>>, N_FEB> where_missing_GTS_trailer2;
int N_Mm_GTS_tag = 0; // Mismatch of the first GTS tag (trailer1, ID 4) between FEBs (only one counted per event)
std::array<std::vector<std::vector<int>>, N_FEB> where_Mm_GTS_tag;
int N_Mm_GTS_time = 0; // Mismatch of the fist GTS time (trailer2, ID 5) between FEBs (only one counted per event)
std::array<std::vector<std::vector<int>>, N_FEB> where_Mm_GTS_time;
int N_MmFEB_GTS_tag = 0; // Mismatch of the GTS tag between header (ID 1) and trailer 1 (ID 5) in a given FEB
std::array<std::vector<std::vector<int>>, N_FEB> where_MmFEB_GTS_tag;
int N_W_GTS1 = 0; // Number of wrong first GTS word (First GTS word should be a trailer 1)
std::array<std::vector<std::vector<int>>, N_FEB> where_W_GTS1;
int N_W_GTS_tag = 0; // Number of wrong GTS tag in trailer 1 w.r.t the previous GTS trailer 1 in a FEB (GTS_tag = n+m for the mth gate trailer 1, with GTS_tag = n for the first trailer 1)
std::array<std::vector<std::vector<int>>, N_FEB> where_W_GTS_tag;
int N_W_GTS_time = 0; // Number of wrong GTS time in trailer 1 w.r.t the previous GTS trailer 1 in a FEB (GTS_time = n+m for the mth gate trailer 1, with GTS_time = n for the first trailer 1, only one error per FEB, per event is counted)
std::array<std::vector<std::vector<int>>, N_FEB> where_W_GTS_time;
int N_W_NGTSperGate = 0; // Number of wrong number of GTS word per Gate (betweem Gate header and gate trailer there should be 4 GTS word of each type = 12 GTS word, only one counted per event)
std::array<std::vector<std::vector<int>>, N_FEB> where_W_NGTSperGate;


//Gate trailers
int N_missing_gate_trailer = 0; // Missing Gate trailer (id 6)
std::array<std::vector<std::vector<int>>, N_FEB> where_missing_gate_trailer;
int N_missing_gate_trailer_time = 0; // Missing Gate time (for trailer, id7)
std::array<std::vector<std::vector<int>>, N_FEB> where_missing_gate_trailer_time;

int N_W_gate_trailer_type = 0; // Number of wrong Gate type (for trailer gate_type == 0)
std::array<std::vector<std::vector<int>>, N_FEB> where_W_gate_trailer_type;
int N_Mm_gate_trailer_number = 0; // Mismatch of the Gate number (trailer) betweem FEBs
std::array<std::vector<std::vector<int>>, N_FEB> where_Mm_gate_trailer_number;
int N_Mm_gate_timeClosedGate = 0; // Mismatch of the Gate time on closed gate (trailer) betweem FEBs
std::array<std::vector<std::vector<int>>, N_FEB> where_Mm_gate_timeClosedGate;


//FEB trailer
int Nmissing_FEB_trailer = 0; // Missing FEB trailer (id6)
std::array<int, N_FEB> slot_missing_FEB_trailer = {0};
std::array<int, N_FEB> slot_EventdoneTimeout = {0};
  //Errors contained in FEB trailer
int N_event_done_timeout = 0;
std::array<int, N_FEB> slot_D1_FIFO_full = {0};
int N_D1_FIFO_full = 0;
std::array<int, N_FEB> slot_D0_FIFO_full = {0};
int N_D0_FIFO_full = 0;
std::array<int, N_FEB> slot_Numb_decoder_errors = {0};
int N_decoder_errors = 0;


// Gate type vs gate number discrepency (last 2 bits of gate number should match gate type)
int N_gate_TagNumb_error = 0;
std::array<std::vector<std::vector<int>>, N_FEB> where_gate_TagNumb_error; //[nevent, n_line]


//OCB trailer
std::vector<std::vector<int>> where_missing_OCB_trailer; // Event number with a missing OCB data packet trailer
std::vector<int> where_gate_open_timeout; // Event number of gate open timeout error
std::vector<int> where_gate_close_timeout; // Event number of gate close timeout error
std::array<std::vector<int>, 14> where_FEB_error; //Number of the event where there is a FEB data packet i error (i is the FEB number)
std::array<int, 14> slot_FEB_error = {0}; //Number of FEB data packet i error (i is the FEB number)
int N_FEB_error = 0; // Total number of "FEB data packet i error" in OCB trailer





//Limit to the number of event read, to accelerate the code during debugging
//to read all event, set to 0
unsigned int to_read_Nevent = 0;
//  If the word-loop exit because it reaches the end of the file (not through to_read_Nevent) it reamins to run the checks on the last event.
// This variable keeps track of that
bool is_Nevent_limit = 0;
//std::vector<int> problematic_event = {};

//Root outputs
TCanvas * c1;

TH1I* hGatetimesOpen[N_FEB];

TH2D* hGateTypeMismatch[N_FEB];
TH2D* hGateNumberMismatch[N_FEB];
TH2D* hGatetimeMismatch[N_FEB];
TH2D* hGatetimeOpengateMismatch[N_FEB];

TH2D* hGateNumberTrailerMismatch[N_FEB];
TH2D* hGatetimeClosedgateMismatch[N_FEB];

TH2D* hGateNumberMismatchFEB;

TH2D* hGTSTagMismatch[N_FEB]; // Mismatch of the fist GTS tag (trailer1, ID 4) between FEBs
TH2D* hGTSTimeMismatch[N_FEB]; // Mismatch of the fist GTS time (trailer2, ID 5) between FEBs
TH2D* hNGTSperGate; // Number of GTS trailer1 per Gate (betweem Gate header and gate trailer there should be 4 GTS trailer 1)
TH2D* hIncrementGTSTime; // Increament between successive GTS time (trailer 2, ID 5)

TH1I* hLGAmplitude[N_FEB];
TH1I* hHGAmplitude[N_FEB];


//Structure containing one event
struct struct_OCB_header{
    int gate_type;
    int gate_tag;
    int nevent;
};

struct struct_gate_header {
    int A_board_id;
    int A_gate_type;
    int A_gate_number;
    int B_board_id;
    int B_gate_time;
    int gate_time;
};

struct struct_hold_time {
    int start_board_id;
    int start_time;
    int stop_board_id;
    int stop_time;
};

struct struct_GTS_header{
    int GTS_tag;
    bool is_inGate;
};

struct struct_hit_time {
    int channel_id;
    int hit_id;
    int tag_id;
    bool edge;
    int hit_time;
};

struct struct_hit_amplitude {
    int channel_id;
    int hit_id;
    int tag_id;
    int amplitude_id;
    int amplitude_meas;
};

struct struct_GTS_trailers{
    int GTS_tag;
    int Data;
    int GTS_time;
    std::array<bool, 2> are_inGate;
};

struct struct_gate_trailer {
    int board_id;
    int gate_type;
    int gate_number;
    int gate_time;
};

struct struct_FEB_trailer {
    bool Event_done_timeout;
    bool D1_FIFO_full;
    bool D0_FIFO_full;
    int N_decoder_error;
};

struct struct_FEB_packet {
    struct_gate_header gate_header;
    bool is_gate_headerA;
    bool is_gate_headerB;
    bool is_gate_header_time;

    struct_hold_time hold_time;
    bool is_hold_start_time;
    bool is_hold_stop_time;

    std::vector<struct_GTS_header> GTS_headers;
    std::vector<int> GTS_IDLog;

    std::vector<struct_hit_time> hits_time;
    std::vector<struct_hit_amplitude> hits_amplitude;

    std::vector<struct_GTS_trailers> GTS_trailers;

    struct_gate_trailer gate_trailer;
    bool is_gate_trailer;
    bool is_gate_trailer_time;

    struct_FEB_trailer FEB_trailer;
    bool is_FEB_trailer;

};

struct struct_OCB_trailer {
    bool gate_open_timeout;
    bool gate_close_timeout;
    std::array<bool, 14> FEB_error = {0};
};

struct struct_event {
    struct_OCB_header OCB_header;
    bool is_OCB_header;
    std::array<struct_FEB_packet, N_FEB> FEBs_packet;
    std::array<bool, N_FEB> is_FEBs_packet;
    struct_OCB_trailer OCB_trailer;
    bool is_OCB_trailer;
};

//vector that holds all the events
std::vector<struct_event> events;


//-----------------------------------------
// Function to return a given bits interval out of a 32 bits words
// input:  uint32_t word: 32 bits word (see Protocol_FEB and ocb_data_format for the details of the bits)
//         int start_bit: number of the first bit to read (0...31)
//         int n_bit: number of bits to read (1...(31-start_bit))
// output: uint32_t : desired bits
// E.g.: (1879048192, 28, 4) -> 7 (Word ID tag)
//-----------------------------------------
int get_bits(uint32_t word, int start_bit, int n_bits) {
    uint32_t mask = ((1 << n_bits) - 1);
    return (word >> start_bit) & mask;
}

// For the struct_FEB_packet initialize the bool values (i.e. set the is_*) of a given slot
void Init_struct_FEB_packet(struct_FEB_packet &FEB_packet){
  FEB_packet.is_gate_headerA = 0;
  FEB_packet.is_gate_headerB = 0;
  FEB_packet.is_gate_header_time = 0;
  FEB_packet.is_hold_start_time = 0;
  FEB_packet.is_hold_stop_time = 0;
  FEB_packet.is_gate_trailer = 0;
  FEB_packet.is_gate_trailer_time = 0;
  FEB_packet.is_FEB_trailer = 0;
}

// Add a event variable to the vector and initialize the bool values (i.e. set the is_* and are_* variables to 0)
void Add_struct_event(std::vector<struct_event> &events){
    struct_event event;
    event.is_OCB_header = 0;
    event.is_FEBs_packet.fill(0);
    event.is_OCB_trailer = 0;
    for (int slot = 0; slot < N_FEB; slot++){
        Init_struct_FEB_packet(event.FEBs_packet[slot]);
    }
    events.push_back(event);
}

// For GTS vectors, initialize the data and then push them in the vector
void Add_struct_GTS_header(struct struct_FEB_packet &FEB_packet){
    struct_GTS_header GTS_header;
    if (FEB_packet.is_gate_trailer == 0){ // No gate trailer yet = the GTS is in the gate
        GTS_header.is_inGate = 1;
    } else {
        GTS_header.is_inGate = 0;
    }
    FEB_packet.GTS_headers.push_back(GTS_header);
}
void Add_struct_GTS_trailers(struct struct_FEB_packet &FEB_packet){
    struct_GTS_trailers GTS_trailers;
    GTS_trailers.are_inGate = {0, 0};
    FEB_packet.GTS_trailers.push_back(GTS_trailers);
}

// Print out the information stored for an event
void printEvent(struct_event &event){
    //OCB header
    printf("OCB header, is_OCB_header: %d\n\tgate type: %d\tgate tag: %d\tnevent: %d\n", event.is_OCB_header, event.OCB_header.gate_type, event.OCB_header.gate_tag, event.OCB_header.nevent);
    //FEBS Packet
    for (int slot = 0; slot < N_FEB; slot++){
        printf("\tFEB %d, is_FEBs_packet: %d\n", slot, event.is_FEBs_packet[slot]);
        //Gate header
        printf("\t\tGate header, is_gate_headerA: %d, is_gate_headerB: %d, is_gate_header_time: %d\n", event.FEBs_packet[slot].is_gate_headerA, event.FEBs_packet[slot].is_gate_headerB, event.FEBs_packet[slot].is_gate_header_time);
        //Gate header A
        printf("\t\t\tGate header A\tBoard id: %d\tGate type: %d\tGate number: %d\n", event.FEBs_packet[slot].gate_header.A_board_id, event.FEBs_packet[slot].gate_header.A_gate_type, event.FEBs_packet[slot].gate_header.A_gate_number);
        //Gate header B
        printf("\t\t\tGate header B\tBoard id: %d\tGate time: %d\n", event.FEBs_packet[slot].gate_header.B_board_id, event.FEBs_packet[slot].gate_header.B_gate_time);
        //Gate header time
        printf("\t\t\tGate header time\tGate time: %d\n", event.FEBs_packet[slot].gate_header.gate_time);
        //Hold time
        printf("\t\tHold time, is_hold_start_time: %d, is_hold_stop_time: %d\n", event.FEBs_packet[slot].is_hold_start_time, event.FEBs_packet[slot].is_hold_stop_time);
        //Hold time start
        printf("\t\t\tHold time start\tBoard ID: %d\tStart time: %d\n", event.FEBs_packet[slot].hold_time.start_board_id, event.FEBs_packet[slot].hold_time.start_time);
        //Hold time stop
        printf("\t\t\tHold time stop\tBoard ID: %d\tStop time: %d\n", event.FEBs_packet[slot].hold_time.stop_board_id, event.FEBs_packet[slot].hold_time.stop_time);
        //Gate trailers
        printf("\t\tGate trailer, is_gate_trailer: %d, is_gate_trailer_time: %d\n", event.FEBs_packet[slot].is_gate_trailer, event.FEBs_packet[slot].is_gate_trailer_time);
        //Gate trailer
        printf("\t\t\tGate trailer\tBoard id: %d\tGate type: %d\tGate number: %d\n", event.FEBs_packet[slot].gate_trailer.board_id, event.FEBs_packet[slot].gate_trailer.gate_type, event.FEBs_packet[slot].gate_trailer.gate_number);
        //Gate trailer time
        printf("\t\t\tGate trailer time\tGate time: %d\n", event.FEBs_packet[slot].gate_trailer.gate_time);
        //FEB trailer
        printf("\t\tFEB trailer, is_FEB_trailer: %d\n", event.FEBs_packet[slot].is_FEB_trailer);
        printf("\t\t\tEvent done timeout: %d\tD1_FIFO_full: %d\tD0_FIFO_full: %d\t N decoder error: %d\n\n", (int)event.FEBs_packet[slot].FEB_trailer.Event_done_timeout, (int)event.FEBs_packet[slot].FEB_trailer.D1_FIFO_full, (int)event.FEBs_packet[slot].FEB_trailer.D0_FIFO_full, event.FEBs_packet[slot].FEB_trailer.N_decoder_error);
    }
    //OCB trailer
    printf("OCB trailer, is_OCB_trailer: %d\n\tgate opeb timeout: %d\tgate close timeout: %d", event.is_OCB_trailer, (int)event.OCB_trailer.gate_open_timeout, (int)event.OCB_trailer.gate_close_timeout);
    for (int slot = 0; slot < N_FEB; slot++){
        printf("\tFEB %d error: %d;", slot, (int)event.OCB_trailer.FEB_error[slot]);
        if ((slot == 3) || (slot == 9)){
          printf("\n");
        }
    }
    printf("\n\n\n");
}

//Controle the structure of the event and print out/save the errors
void Run_checks(std::vector<struct_event> &events){
    Nevent_checked++;
    // Check if the previous event was properly closed (OCB trailer was read)
    int ref_gate_type = -1;
    if(events.back().is_OCB_trailer == 0){
        n_line++;
        fprintf(write_fp_words, "%d\t>> ERROR: OCB trailer not found in Event %d <<\n", n_line, events.back().OCB_header.nevent);
        std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
        where_missing_OCB_trailer.push_back(location);
        Nerrors++; //add an entry to the error counter
    } else { // if there was a trailer count the errors contain in it
        if (events.back().OCB_trailer.gate_open_timeout == 1){
            where_gate_open_timeout.push_back(events.back().OCB_header.nevent);
            Nerrors++;
        }
        if (events.back().OCB_trailer.gate_close_timeout == 1){
            where_gate_close_timeout.push_back(events.back().OCB_header.nevent);
            Nerrors++;
        }
        //Loop over the FEB data packet error
        for (int i = 0; i < N_FEB; i++){
            if (events.back().OCB_trailer.FEB_error[i] == 1){
                where_FEB_error[i].push_back(events.back().OCB_header.nevent);
                N_FEB_error++;
            }
        }
    }

    // Then proceed with the checks over the last event:
    // Define some reference values
    int ref_GTS_tag = -1;
    int ref_GTS_time = -1;
    bool is_Mm_GTS_tag_counted = false; // bool to count only one GTS tag mismatch per event
    bool is_Mm_GTS_time_counted = false; // bool to count only one GTS time mismatch per event
    bool is_W_NGTSperGate_counted = false; // bool to count only one wrong GTS number per event

    // Loop over the FEB slots
    for (int slot = 0; slot < N_FEB; slot++){
        //Check if the gate headers were correctly read, if yes perform the check
        if (events.back().FEBs_packet[slot].is_gate_headerA != 1){
            n_line++;
            fprintf(write_fp_words, "%d\t>> ERROR: The Gate header A is missing for the FEB %d<<\n", n_line,  slot);
            std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
            where_missing_gate_headerA[slot].push_back(location);
            N_missing_gate_headerA++;
            Nerrors++;
        }else{
            // Gate tag vs gate number discrepency (last 2 bits of gate number (ID 0) should match gate type (of OCB header, ID 8))
            if (events.back().OCB_header.gate_tag != get_bits(events.back().FEBs_packet[slot].gate_header.A_gate_number, 0,2)){
                n_line++;
                fprintf(write_fp_words, "%d\t>> ERROR: The gate number (%d) of FEB %d does not match the gate tag (%d) of the OCB header <<\n", n_line, get_bits(events.back().FEBs_packet[slot].gate_header.A_gate_number, 0,2), slot, events.back().OCB_header.gate_tag);
                N_gate_TagNumb_error++ ;
                std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
                where_gate_TagNumb_error[slot].push_back(location);
                Nerrors++;
            }
            // Check for Mismatch Gate type and number (header A) between Febs and plot them
            int ref_gate_type = events.back().FEBs_packet[slot].gate_header.A_gate_type;
            int ref_gate_number = events.back().FEBs_packet[slot].gate_header.A_gate_number;
            if (slot == 0){
                bool is_Mm_gate_type = 0; // Keep track is the mismatch was already counted in the error count (we want only one mismatch per event counted)
                bool is_Mm_gate_number = 0;
                for (int j_FEB = 1; j_FEB < N_FEB; j_FEB++){
                    hGateTypeMismatch[slot]->Fill(j_FEB, ref_gate_type - events.back().FEBs_packet[j_FEB].gate_header.A_gate_type); // Fill the histo
                    hGateNumberMismatch[slot]->Fill(j_FEB, ref_gate_number - events.back().FEBs_packet[j_FEB].gate_header.A_gate_number);
                    if ((events.back().FEBs_packet[j_FEB].gate_header.A_gate_type != ref_gate_type) && (is_Mm_gate_type == 0)){
                        n_line++;
                        fprintf(write_fp_words, "%d\t>> ERROR: The gate type (header A, id0) are mismatched between FEBs (e.g. FEB %i (%i) and %i (%i))<<\n", n_line, slot, ref_gate_type, j_FEB,  events.back().FEBs_packet[j_FEB].gate_header.A_gate_type);
                        N_Mm_gate_type++;
                        Nerrors++;
                        std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
                        where_Mm_gate_type[j_FEB].push_back(location);
                        is_Mm_gate_type = 1;
                    }
                    if ((events.back().FEBs_packet[j_FEB].gate_header.A_gate_number != ref_gate_number) && (is_Mm_gate_number == 0)){
                        n_line++;
                        fprintf(write_fp_words, "%d\t>> ERROR: The gate number (header A, id0) are mismatched between FEBs (e.g. FEB %i (%i) and %i (%i))<<\n", n_line, slot, ref_gate_number, j_FEB,  events.back().FEBs_packet[j_FEB].gate_header.A_gate_number);
                        N_Mm_gate_number++;
                        Nerrors++;
                        std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
                        where_Mm_gate_number[j_FEB].push_back(location);
                        is_Mm_gate_number = 1;
                    }
                }
            } else { // Fill the other histo (the Mismatch w.r.t each of the FEBs)
                for (int j_FEB = 0; j_FEB < N_FEB; j_FEB++){
                    if (j_FEB == slot) continue;
                    hGateTypeMismatch[slot]->Fill(j_FEB, ref_gate_type - events.back().FEBs_packet[j_FEB].gate_header.A_gate_type);
                    hGateNumberMismatch[slot]->Fill(j_FEB, ref_gate_number - events.back().FEBs_packet[j_FEB].gate_header.A_gate_number);
                }
            } // Check for mismatch between the gate number of the header (A) and trailer
            if (events.back().FEBs_packet[slot].is_gate_trailer == 1){
                hGateNumberMismatchFEB->Fill(slot, events.back().FEBs_packet[slot].gate_header.A_gate_number - events.back().FEBs_packet[slot].gate_trailer.gate_number);
                if (events.back().FEBs_packet[slot].gate_header.A_gate_number != events.back().FEBs_packet[slot].gate_trailer.gate_number){
                    n_line++;
                    fprintf(write_fp_words, "%d\t>> ERROR: The Gate number of the header (header A, ID 0) %d does not match the one of the trailer (ID 6) %d for the FEB %d <<\n", n_line, events.back().FEBs_packet[slot].gate_header.A_gate_number, events.back().FEBs_packet[slot].gate_trailer.gate_number, slot);
                    Nerrors++;
                    N_MmFEB_gate_number++;
                    std::vector<int> location = {events.back().OCB_header.nevent, n_line};
                    where_MmFEB_gate_number[slot].push_back(location);
                }
            }
        }
        if (events.back().FEBs_packet[slot].is_gate_headerB == 0){
            n_line++;
            fprintf(write_fp_words, "%d\t>> ERROR: The Gate header B (ID 0) is missing for the FEB %d <<\n", n_line, slot);
            std::vector<int> location = {events.back().OCB_header.nevent, n_line};
            where_missing_gate_headerB[slot].push_back(location);
            N_missing_gate_headerB++;
            Nerrors++;
        } else {
          // Check for Mismatch Gate time from GTS (header B) between Febs and plot them
          int ref_gate_time = events.back().FEBs_packet[slot].gate_header.B_gate_time;
          if (slot == 0){
              bool is_Mm_gate_time = 0; // Keep track is the Mismatch was already counted in the error count (we want only one Mismatch per event counted)
              for (int j_FEB = 1; j_FEB < N_FEB; j_FEB++){
                  hGatetimeMismatch[slot]->Fill(j_FEB, ref_gate_time - events.back().FEBs_packet[j_FEB].gate_header.B_gate_time); // Fill the histo
                  if ((events.back().FEBs_packet[j_FEB].gate_header.B_gate_time != ref_gate_time) && (is_Mm_gate_time == 0)){
                      n_line++;
                      fprintf(write_fp_words, "%d\t>> ERROR: The gate times from GTS (header B, ID 0) are mismatched between FEBs (e.g. FEB %i (%i) and %i (%i))<<\n", n_line, slot, ref_gate_time, j_FEB,  events.back().FEBs_packet[j_FEB].gate_header.B_gate_time);
                      N_Mm_gate_time++;
                      Nerrors++;
                      std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
                      where_Mm_gate_time[j_FEB].push_back(location);
                      is_Mm_gate_time = 1;
                  }
              }
          } else { //Fill the other histo (the Mismatch w.r.t each of the FEBs)
              for (int j_FEB = 0; j_FEB < N_FEB; j_FEB++){
                  if (j_FEB == slot) continue;
                  hGatetimeMismatch[slot]->Fill(j_FEB, ref_gate_time - events.back().FEBs_packet[j_FEB].gate_header.B_gate_time);

              }
          }
        }
        if (events.back().FEBs_packet[slot].is_gate_header_time == 0){
            n_line++;
            fprintf(write_fp_words, "%d\t>> ERROR: The Gate time word of Gate header (ID 7) is missing for the FEB %d <<\n", n_line,  slot);
            std::vector<int> location = {events.back().OCB_header.nevent, n_line};
            where_missing_gate_header_time[slot].push_back(location);
            N_missing_gate_header_time++;
            Nerrors++;
        } else { //if there is an header, fill the hGatetimesOpen
            hGatetimesOpen[slot]->Fill(events.back().FEBs_packet[slot].gate_header.gate_time);
            int ref_gate_timeOpenGate = events.back().FEBs_packet[slot].gate_header.gate_time;
            if (slot == 0){
                bool is_Mm_gate_timeOpenGate = 0; // Keep track is the mismatch was already counted in the error count (we want only one mismatch per event counted)
                for (int j_FEB = 1; j_FEB < N_FEB; j_FEB++){
                    hGatetimeOpengateMismatch[slot]->Fill(j_FEB, ref_gate_timeOpenGate - events.back().FEBs_packet[j_FEB].gate_header.gate_time);
                    if ((events.back().FEBs_packet[j_FEB].gate_header.gate_time != ref_gate_timeOpenGate) && (is_Mm_gate_timeOpenGate == 0)){
                        n_line++;
                        fprintf(write_fp_words, "%d\t>> ERROR: The gate time on open gate (Gate time, ID 7) are mismatched between FEBs (e.g. FEB %i (%i) and %i (%i))<<\n", n_line, slot, ref_gate_timeOpenGate, j_FEB,  events.back().FEBs_packet[j_FEB].gate_header.gate_time);
                        N_Mm_gate_timeOpenGate++;
                        Nerrors++;
                        std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
                        where_Mm_gate_timeOpenGate[j_FEB].push_back(location);
                        is_Mm_gate_timeOpenGate = 1;
                    }
                }
            } else { //Fill the other histo (the Mismatch w.r.t each of the FEBs)
                for (int j_FEB = 0; j_FEB < N_FEB; j_FEB++){
                    if (j_FEB == slot) continue;
                    hGatetimeOpengateMismatch[slot]->Fill(j_FEB, ref_gate_timeOpenGate - events.back().FEBs_packet[j_FEB].gate_header.gate_time);
                }
            }
        }

        //Check if the GTS Log to see if it starts with a trailer 1 (ID 4) and if a header/trailer is missing (correct order is 4->5->1->4->...)
        if (events.back().FEBs_packet[slot].GTS_IDLog[0] != 4){
          n_line++;
          fprintf(write_fp_words, "%d\t>> ERROR: A the first GTS word (ID %d) is not a trailer 1 (ID 4) word for FEB %d <<\n", n_line, events.back().FEBs_packet[slot].GTS_IDLog[0] , slot);
          N_W_GTS1++;
          std::vector<int> location = {events.back().OCB_header.nevent, n_line};
          where_W_GTS1[slot].push_back(location);
          Nerrors++;
        }
        for (int i = 1; i < events.back().FEBs_packet[slot].GTS_IDLog.size(); i++) {
            if (events.back().FEBs_packet[slot].GTS_IDLog[i] == 5){
                if (events.back().FEBs_packet[slot].GTS_IDLog[i-1] != 4){
                    n_line++;
                    fprintf(write_fp_words, "%d\t>> ERROR: A GTS trailer 1 (ID 4) is missing for FEB %d (associated to GTS tag: %d) <<\n");
                    N_missing_GTS_trailer1++;
                    std::vector<int> location = {events.back().OCB_header.nevent, n_line};
                    where_missing_GTS_trailer1[slot].push_back(location);
                    Nerrors++;
                }
            }else if (events.back().FEBs_packet[slot].GTS_IDLog[i] == 1){
                if (events.back().FEBs_packet[slot].GTS_IDLog[i-1] != 5){
                    n_line++;
                    fprintf(write_fp_words, "%d\t>> ERROR: A GTS trailer 2 (ID 5) is missing for FEB %d (associated to GTS tag: %d) <<\n");
                    N_missing_GTS_trailer2++;
                    std::vector<int> location = {events.back().OCB_header.nevent, n_line};
                    where_missing_GTS_trailer2[slot].push_back(location);
                    Nerrors++;
                }
            }else if (events.back().FEBs_packet[slot].GTS_IDLog[i] == 4){
                if (events.back().FEBs_packet[slot].GTS_IDLog[i-1] != 1){
                    n_line++;
                    fprintf(write_fp_words, "%d\t>> ERROR: A GTS header (ID 1) is missing  for FEB %d (associated to GTS tag: %d) <<\n");
                    N_missing_GTS_header++;
                    std::vector<int> location = {events.back().OCB_header.nevent, n_line};
                    where_missing_GTS_header[slot].push_back(location);
                    Nerrors++;
                }
            }
        }
        // Check for mismatched first GTS tag and time between FEBs
        if ((ref_GTS_tag == -1) || (ref_GTS_time == -1)){
            ref_GTS_tag = events.back().FEBs_packet[slot].GTS_trailers[0].GTS_tag;
            ref_GTS_time = events.back().FEBs_packet[slot].GTS_trailers[0].GTS_time;
        } else {
            if ((is_Mm_GTS_tag_counted == false) && (events.back().FEBs_packet[slot].GTS_trailers[0].GTS_tag != ref_GTS_tag)){
                n_line++;
                fprintf(write_fp_words, "%d\t>> ERROR: The GTS tag (%d) of the first trailer1 (ID 4) are mismatched between FEBs (ref. GTS tag = %d) <<\n",n_line, events.back().FEBs_packet[slot].GTS_trailers[0].GTS_tag, ref_GTS_tag);
                N_Mm_GTS_tag++;
                std::vector<int> location = {events.back().OCB_header.nevent, n_line};
                where_Mm_GTS_tag[slot].push_back(location);
                Nerrors++;
                is_Mm_GTS_tag_counted = true; //Count only one error per event
            }
            if ((is_Mm_GTS_time_counted == false) && (events.back().FEBs_packet[slot].GTS_trailers[0].GTS_time != ref_GTS_time)){
                n_line++;
                fprintf(write_fp_words, "%d\t>> ERROR: The GTS time (%d) of the first trailer2 (ID 5) are mismatched between FEBs (ref. GTS time = %d) <<\n",n_line, events.back().FEBs_packet[slot].GTS_trailers[0].GTS_time, ref_GTS_time);
                N_Mm_GTS_time++;
                std::vector<int> location = {events.back().OCB_header.nevent, n_line};
                where_Mm_GTS_time[slot].push_back(location);
                Nerrors++;
                is_Mm_GTS_time_counted = true; //Count only one error per event
            }
        }
        //Fill the histos
        for (int i_FEB = 0; i_FEB < N_FEB; i_FEB++){
            if (i_FEB != slot){
                hGTSTagMismatch[slot]->Fill(i_FEB, events.back().FEBs_packet[slot].GTS_trailers[0].GTS_tag - events.back().FEBs_packet[i_FEB].GTS_trailers[0].GTS_tag);
                hGTSTimeMismatch[slot]->Fill(i_FEB, events.back().FEBs_packet[slot].GTS_trailers[0].GTS_time - events.back().FEBs_packet[i_FEB].GTS_trailers[0].GTS_time);
            }
        }
        // Check if there is 4 Gates of each type between the Gate header and trailer
        int N_GTStrailer1inGate = 0;
        int N_GTStrailer2inGate = 0;
        int N_GTSheaderinGate = 0;
        bool is_W_GTS_tag_counted = false; // bool to count only one wrong GTS tag per FEB and per event
        bool is_W_GTS_time_counted = false; // bool to count only one wrong GTS time per FEB and per event
        if (events.back().FEBs_packet[slot].GTS_trailers.size() == events.back().FEBs_packet[slot].GTS_headers.size()){ //check tha there are the same number of GTS header and trailers
            int N_GTS = events.back().FEBs_packet[slot].GTS_trailers.size();
            N_GTS_trailer1_checked[slot] += N_GTS;
            N_GTStrailer1inGate += events.back().FEBs_packet[slot].GTS_trailers[0].are_inGate[0];
            N_GTStrailer2inGate += events.back().FEBs_packet[slot].GTS_trailers[0].are_inGate[1];
            N_GTSheaderinGate += events.back().FEBs_packet[slot].GTS_headers[0].is_inGate;
            for (int i_GTS = 1; i_GTS < N_GTS; i_GTS++){
                N_GTStrailer1inGate += events.back().FEBs_packet[slot].GTS_trailers[i_GTS].are_inGate[0];
                N_GTStrailer2inGate += events.back().FEBs_packet[slot].GTS_trailers[i_GTS].are_inGate[1];
                N_GTSheaderinGate += events.back().FEBs_packet[slot].GTS_headers[i_GTS].is_inGate;
                // Check if the GTS tag and time increaments are correct (should be +1)
                if ((is_W_GTS_tag_counted == false) && (events.back().FEBs_packet[slot].GTS_trailers[i_GTS].GTS_tag != (events.back().FEBs_packet[slot].GTS_trailers[i_GTS-1].GTS_tag+1))) {
                    n_line++;
                    fprintf(write_fp_words, "%d\t>> ERROR: The GTS tag of trailer 1 (ID 4) went from %d to %d for the FEB %d <<\n", n_line, events.back().FEBs_packet[slot].GTS_trailers[i_GTS-1].GTS_tag, events.back().FEBs_packet[slot].GTS_trailers[i_GTS].GTS_tag, slot);
                    N_W_GTS_tag++;
                    std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
                    where_W_GTS_tag[slot].push_back(location);
                    Nerrors++;
                    is_W_GTS_tag_counted = true; //Count only one error per event
                }
                if ((is_W_GTS_time_counted == false) && (events.back().FEBs_packet[slot].GTS_trailers[i_GTS].GTS_time != (events.back().FEBs_packet[slot].GTS_trailers[i_GTS-1].GTS_time+1))) {
                    n_line++;
                    fprintf(write_fp_words, "%d\t>> ERROR: The GTS time of trailer 1 (ID 4) went from %d to %d for the FEB %d <<\n", n_line, events.back().FEBs_packet[slot].GTS_trailers[i_GTS-1].GTS_time, events.back().FEBs_packet[slot].GTS_trailers[i_GTS].GTS_time, slot);
                    N_W_GTS_time++;
                    std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
                    where_W_GTS_time[slot].push_back(location);
                    Nerrors++;
                    is_W_GTS_time_counted = true; //Count only one error per event
                }
                // Check if the GTS tag of the header match the GTS tag of the trailer 1 (Note that the first GTS word should be a trailer 1, not a header)
                if (events.back().FEBs_packet[slot].GTS_headers[i_GTS-1].GTS_tag != events.back().FEBs_packet[slot].GTS_trailers[i_GTS].GTS_tag){
                    n_line++;
                    fprintf(write_fp_words, "%d\t>> ERROR: The GTS tags are mismatch between GTS header (ID 1) and GTS trailer 1 (ID 4) (e.g. %d and %d) for the FEB %d <<\n", n_line, events.back().FEBs_packet[slot].GTS_headers[i_GTS-1], events.back().FEBs_packet[slot].GTS_trailers[i_GTS].GTS_tag, slot);
                    N_MmFEB_GTS_tag++;
                    std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
                    where_MmFEB_GTS_tag[slot].push_back(location);
                    Nerrors++;
                }
                //Fill histo, use first and last bin as overflow bins
                if (((events.back().FEBs_packet[slot].GTS_trailers[i_GTS].GTS_time - events.back().FEBs_packet[slot].GTS_trailers[i_GTS-1].GTS_time) <= 5) && ((events.back().FEBs_packet[slot].GTS_trailers[i_GTS].GTS_time - events.back().FEBs_packet[slot].GTS_trailers[i_GTS-1].GTS_time) >= -5)){
                    hIncrementGTSTime->Fill(slot, events.back().FEBs_packet[slot].GTS_trailers[i_GTS].GTS_time - events.back().FEBs_packet[slot].GTS_trailers[i_GTS-1].GTS_time);
                } else if ((events.back().FEBs_packet[slot].GTS_trailers[i_GTS].GTS_time - events.back().FEBs_packet[slot].GTS_trailers[i_GTS-1].GTS_time) < -5) {
                    hIncrementGTSTime->Fill(slot, -5);
                } else if ((events.back().FEBs_packet[slot].GTS_trailers[i_GTS].GTS_time - events.back().FEBs_packet[slot].GTS_trailers[i_GTS-1].GTS_time) > 5) {
                    hIncrementGTSTime->Fill(slot, 5);
                }
            }
            if ((is_W_NGTSperGate_counted == false) && ((N_GTStrailer1inGate != 4) || (N_GTStrailer2inGate != 4) || (N_GTSheaderinGate != 4))) {
                n_line++;
                fprintf(write_fp_words, "%d\t>> ERROR: The number of GTS word (ID 1, 4, 5) between the Gate header and trailer is wrong (trailer 1: %d, trailer 2: %d, header: %d, expect 4 for each) for the FEB %d <<\n", n_line, N_GTStrailer1inGate, N_GTStrailer2inGate, N_GTSheaderinGate, slot);
                N_W_NGTSperGate++;
                std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
                where_W_NGTSperGate[slot].push_back(location);
                Nerrors++;
                is_W_NGTSperGate_counted = true; //Count only one error per event
            }
            //Fill histo
            hNGTSperGate->Fill(slot, N_GTStrailer1inGate);
        } else { //If the number of trailer does not match the number of header, print out an error to notify that the checks on GTS were not performed
            n_line++;
            fprintf(write_fp_words, "%d\t>> ERROR: The number of GTS trailers does not match the number of header for the FEB %d, therefore no check was performed about GTS information <<\n", n_line, N_GTStrailer1inGate, N_GTStrailer2inGate, N_GTSheaderinGate, slot);
        }

        //Check the gate trailer
        if (events.back().FEBs_packet[slot].is_gate_trailer == 0){
            n_line++;
            fprintf(write_fp_words, "%d\t>> ERROR: The gate trailer (ID 6) is missing for the FEB %d <<\n", n_line, slot);
            N_missing_gate_trailer++;
            std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
            where_missing_gate_trailer[slot].push_back(location);
            Nerrors++;
        } else {
            // Check for wrong gate type (trailer gate type == 0)
            if (events.back().FEBs_packet[slot].gate_trailer.gate_type != 0){
              n_line++;
              fprintf(write_fp_words, "%d\t>> ERROR: The gate type of the gate trailer (ID 6) is %d and  not 0 for the FEB %d <<\n", n_line, events.back().FEBs_packet[slot].gate_trailer.gate_type, slot);
              N_W_gate_trailer_type++;
              std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
              where_W_gate_trailer_type[slot].push_back(location);
              Nerrors++;
            }
            // Check for Mismatch Gate type and number (trailer) between Febs and plot them
            int ref_gate_trailer_number = events.back().FEBs_packet[slot].gate_trailer.gate_number;
            if (slot == 0){
                bool is_Mm_gate_trailer_number = 0;
                for (int j_FEB = 1; j_FEB < N_FEB; j_FEB++){
                    hGateNumberTrailerMismatch[slot]->Fill(j_FEB, ref_gate_trailer_number - events.back().FEBs_packet[j_FEB].gate_trailer.gate_number);
                    if ((events.back().FEBs_packet[j_FEB].gate_trailer.gate_number != ref_gate_trailer_number) && (is_Mm_gate_trailer_number == 0)){
                        n_line++;
                        fprintf(write_fp_words, "%d\t>> ERROR: The gate number (trailer, ID 6) are mismatched between FEBs (e.g. FEB %i (%i) and %i (%i)) <<\n", n_line, slot, ref_gate_trailer_number, j_FEB,  events.back().FEBs_packet[j_FEB].gate_trailer.gate_number);
                        N_Mm_gate_trailer_number++;
                        Nerrors++;
                        std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
                        where_Mm_gate_trailer_number[j_FEB].push_back(location);
                        is_Mm_gate_trailer_number = 1;
                    }
                }
            } else { //Fill the other histo (the Mismatch w.r.t each of the FEBs)
                for (int j_FEB = 0; j_FEB < N_FEB; j_FEB++){
                    if (j_FEB == slot) continue;
                    hGateNumberTrailerMismatch[slot]->Fill(j_FEB, ref_gate_trailer_number - events.back().FEBs_packet[j_FEB].gate_trailer.gate_number);
                }
            }
        }
        if (events.back().FEBs_packet[slot].is_gate_trailer_time == 0){
            n_line++;
            fprintf(write_fp_words, "%d\t>> ERROR: The Gate time word of Gate trailer (ID 6) is missing for the FEB %d <<\n", n_line, slot);
            N_missing_gate_trailer_time++;
            std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
            where_missing_gate_trailer_time[slot].push_back(location);
            Nerrors++;
        } else {
            int ref_gate_timeClosedGate = events.back().FEBs_packet[slot].gate_trailer.gate_time;
            if (slot == 0){
                bool is_Mm_gate_timeClosedGate = 0; // Keep track is the mismatch was already counted in the error count (we want only one mismatch per event counted)
                for (int j_FEB = 1; j_FEB < N_FEB; j_FEB++){
                    hGatetimeClosedgateMismatch[slot]->Fill(j_FEB, ref_gate_timeClosedGate - events.back().FEBs_packet[j_FEB].gate_trailer.gate_time);
                    if ((events.back().FEBs_packet[j_FEB].gate_trailer.gate_time != ref_gate_timeClosedGate) && (is_Mm_gate_timeClosedGate == 0)){
                        n_line++;
                        fprintf(write_fp_words, "%d\t>> ERROR: The gate time on closed gate (Gate time, ID 7) are mismatched between FEBs (e.g. FEB %i (%i) and %i (%i))<<\n", n_line, slot, ref_gate_timeClosedGate, j_FEB,  events.back().FEBs_packet[j_FEB].gate_trailer.gate_type);
                        N_Mm_gate_timeClosedGate++;
                        Nerrors++;
                        std::vector<int> location = {events.back().OCB_header.nevent, n_line}; //[nevent, n_line]
                        where_Mm_gate_timeClosedGate[j_FEB].push_back(location);
                        is_Mm_gate_timeClosedGate = 1;
                    }
                }
            } else { //Fill the other histo (the Mismatch w.r.t each of the FEBs)
                for (int j_FEB = 0; j_FEB < N_FEB; j_FEB++){
                    if (j_FEB == slot) continue;
                    hGatetimeClosedgateMismatch[slot]->Fill(j_FEB, ref_gate_timeClosedGate - events.back().FEBs_packet[j_FEB].gate_trailer.gate_time);
                }
            }

        }
        //Check if the FEB trailer was correctly read
        if(events.back().FEBs_packet[slot].is_FEB_trailer == 0){ //Check if missing
            n_line++;
            fprintf(write_fp_words, "%d\t>> ERROR: The FEB trailer (ID 6) is missing for the FEB in slot %d <<\n", n_line, slot);
            Nmissing_FEB_trailer++;
            slot_missing_FEB_trailer[slot]++;
            Nerrors++;
        }else{ //if not missing, record the error it outputs
            if (events.back().FEBs_packet[slot].FEB_trailer.Event_done_timeout == 1){
                slot_EventdoneTimeout[slot]++;
                N_event_done_timeout++;
                Nerrors++;
            }
            if (events.back().FEBs_packet[slot].FEB_trailer.D1_FIFO_full == 1){
                slot_D1_FIFO_full[slot]++;
                N_D1_FIFO_full++;
                Nerrors++;
            }
            if (events.back().FEBs_packet[slot].FEB_trailer.D0_FIFO_full == 1){
                slot_D0_FIFO_full[slot]++;
                N_D0_FIFO_full++;
                Nerrors++;
            }
            slot_Numb_decoder_errors[slot] += events.back().FEBs_packet[slot].FEB_trailer.N_decoder_error;
            N_decoder_errors += events.back().FEBs_packet[slot].FEB_trailer.N_decoder_error;
            Nerrors += events.back().FEBs_packet[slot].FEB_trailer.N_decoder_error;
        }

        //loop over the GTS trailer
        // for (int i_GTStrailer = 0; i_GTStrailer < events.back().FEBs_packet[i_FEB].GTS_trailers.size(); i_GTStrailer++){
        //     //Ref. GTS time
        //     int ref_GTSTime = events.back().FEBs_packet[i_FEB].GTS_trailers[i_GTStrailer].GTS_time;
        //     int ref_GTSTag = events.back().FEBs_packet[i_FEB].GTS_trailers[i_GTStrailer].GTS_tag;
        //     //loop over the FEB a second time to compute the GTS time differences w.r.t every FEBS
        //     for (int j_FEB = 0; j_FEB < N_FEB; j_FEB++){
        //         events.back().FEBs_packet[i_FEB].GTS_trailers[i_GTStrailer].GTS_time;
        //     }
        // }
    }
    //Check if the event number increment is correct (if the nevent[i] = nevent[i-1]+1)
    if ((events.size() > 1) && (events[Nevent-2].OCB_header.nevent+1 != events[Nevent-1].OCB_header.nevent)){
        n_line++;
        fprintf(write_fp_words, "%d\t>> ERROR: Event number went from %d to %d <<\n", n_line, events[Nevent-2].OCB_header.nevent, events[Nevent-1].OCB_header.nevent);
        std::array<int, 3> where_skip_event = {events[Nevent-2].OCB_header.nevent, events[Nevent-1].OCB_header.nevent, n_line};
        where_skip_events.push_back(where_skip_event);
        Nskip_event++;
        Nerrors++;
    }
}


//print the array of type std::array<int, 14> (most of the slot_* variables) in write_fp
void printArray(FILE *write_fp, const std::array<int, 14> &arr){
    for(int k=0; k< 14; k++){
        if(k==0){
            fprintf(write_fp,"{%i, ", arr[k]);
        }else if(k==13){
            fprintf(write_fp, "%i}\n", arr[k]);
        }else{
            fprintf(write_fp, "%i, ", arr[k]);
        }
    }
}

//print the array of type std::array<std::vector<std::vector<int>>, 14> (most of the where_* variables) in write_fp
void printArrayofVector(FILE *write_fp, const std::array<std::vector<std::vector<int>>, 14> &arr){
    fprintf(write_fp_errors, "| Per slots: \n");
    for(int slot = 0; slot < 14; slot++){
        if(arr[slot].size() != 0){
            if (arr[slot].size() < 10) {
                fprintf(write_fp,"|\tFEB %i [event (line)]: ", slot);
                for (int n = 0; n < arr[slot].size()-1; n++){
                    fprintf(write_fp,"%i (%i), ", arr[slot][n][0], arr[slot][n][1]);
                }
            } else {
                fprintf(write_fp,"|\tFEB %i [event (line)] (%i errors, 10 are displayed): ", slot, arr[slot].size());
                for (int n = 0; n < 9; n++){
                    fprintf(write_fp,"%i (%i), ", arr[slot][n][0], arr[slot][n][1]);
                }
            }
            fprintf(write_fp,"%i (%i)\n", arr[slot].back()[0], arr[slot].back()[1]);
        }else{
            fprintf(write_fp,"|\tFEB %i : none\n", slot);
        }
    }
}

//print the array of type std::vector<std::array<int, 3>> in write_fp
void printVectorofArray3(FILE *write_fp, const std::vector<std::array<int, 3>> &vec){
    for (int i = 0; i < vec.size(); i++){
        if (i%10 == 0){
            fprintf(write_fp, "\n| [%d, %d, %d]", vec[i][0], vec[i][1], vec[i][2]);
        } else{
            fprintf(write_fp, " ,[%d, %d, %d]", vec[i][0], vec[i][1], vec[i][2]);
        }
    }
    fprintf(write_fp, "\n");
}


//Print out a std::vector<std::vector<int>> struct (e.g. where_missing_OCB_trailer (NO FEB info))
void printVectorofVector(FILE *write_fp, std::vector<std::vector<int>> &vecofvec){
    int limit = std::min(10, (int)vecofvec.size());
    if (limit != 0){
        if (limit == 10){
            fprintf(write_fp,"| Detailed: [event (line)] (%d errors, 10 are displayed): ", limit);
        } else {
            fprintf(write_fp,"| Detailed: [event (line)]: ");
        }
        fprintf(write_fp,"|\t");
        for (int i = 0; i < limit-1; i++){
            fprintf(write_fp,"%i (%i), ", vecofvec[i][0], vecofvec[i][1]);
        }
        fprintf(write_fp,"%i (%i)\n", vecofvec.back()[0], vecofvec.back()[1]);
    }
}



//-----------------------------------------
// Briefly: input a .bin file and check that its structure follow the designed model, output the error detected
//
// The .bin files is read from an OCB header to an OCN trailer, the information are stored in the dedicated structure.
// Once an event is completly read, some checks are run to output the error in the structure. This is repeted until
// the end of the .bin file or untile the umber of event read match to_read_Nevent parameter. Then the errors are outputed
//through a .txt and some plots (.root and .pdf files)
// input: <file_name>.bin: binary file to decode and check
//        [output-path] : directory to put the result in it (optional)
// output: <file_name>_LogWords.txt: Log with the words decoded
//         <file_name>_LogErrors.txt: Log with the errors (summary and detailed)
//         <file_name>_histo.root: Root file with the plots
//         <file_name>_histo.pdf: pdf file with the plots
//-----------------------------------------
int main (int argc, char *argv[]) {
    char input_file[100]; //binary file (with its path)
    char file_name[100]; // file name (with extension)
    char output_path[100]; //path + name_file (without extension) for the output

    // Check the passed arguments and if the first one is a binary
    if (argc < 2) {
        printf("Binary file to decode not given!\n");
        return -1;
    } else if (argc == 2) { //output path not given, output in the same folder as the binary one
        if (strstr(argv[1], ".bin") == NULL) {
                printf("Input file must be a .bin file!\n");
                return -1;
        }
        strcpy(input_file, argv[1]);
        strcpy(output_path, argv[1]);

    } else if (argc == 3) { //path given, add the file name to the output path
        if (strstr(argv[1], ".bin") == NULL) {
            printf("Input file must be a .bin file!\n");
            return -1;
        }
        strcpy(input_file, argv[1]);
        strcpy(file_name, argv[1]);
        strcpy(output_path, argv[2]);
        char * pLastSlash; //pointer to the last slash in path
        pLastSlash = strrchr(output_path,'/');
        if (pLastSlash+1 != NULL){ //check if the last slash is after the output folder, if not add a slash
            if (*(pLastSlash+1) != '\0'){
                strcat(output_path, "/");
            }
        }
        char * pFileName;
        //Split the path into token w.r.t '/' and keep the last one (== file name)
        pFileName = strtok(file_name,"/");
        while (pFileName != NULL) {
            strcpy(file_name, pFileName);
            pFileName = strtok(NULL, "/");
        }
        strcat(output_path, file_name);
    } else {
        printf("Invalid number of arguments!\n");
        return -1;
    }
    //Remove the .bin extension
    char * pExt;
    pExt = strrchr(output_path,'.');
    strncpy(pExt,"\0",4); //get rid of the .bin

    // Set the different file name
    char output_LogWords[100];
    strcpy(output_LogWords, output_path);
    strcat(output_LogWords,"_LogWords.txt");

    char output_LogErrors[100];
    strcpy(output_LogErrors, output_path);
    strcat(output_LogErrors,"_LogErrors.txt");

    char output_root[100];
    strcpy(output_root, output_path);
    strcat(output_root,"_histo.root");

    char output_pdf[100];
    strcpy(output_pdf, output_path);
    strcat(output_pdf,"_histo.pdf");

    //Create the output text files
    std::ofstream{output_LogWords};
    std::ofstream{output_LogErrors};

    // Open file
    unsigned char currline[buflen];
    read_fp = fopen(input_file, "rb");
    write_fp_words = fopen(output_LogWords, "w");
    write_fp_errors = fopen(output_LogErrors, "w");

    if (read_fp == NULL) {
        printf("Error opening read file %s\n", input_file);
        return -1;
    }
    if (write_fp_words == NULL) {
        printf("Error opening write file %s\n", output_LogWords);
        return -1;
    }
    if (write_fp_errors == NULL) {
        printf("Error opening write file %s\n", output_LogErrors);
        return -1;
    }

    // Create the root file
    TFile *file = new TFile(output_root, "RECREATE");
    file->cd();

    // Initialize the canvas and the histograms
    c1 = new TCanvas("c1", "Output Histograms", 900, 600); //1920 x 1080
    int xlow = 0;
    int xup = 10000000; //60000; //Gate time is a 28-bits variable (268435455)
    int bin_size = 100;
    int nbinsx = (xup-xlow)/bin_size;
    for (int i_hist = 0; i_hist < N_FEB; i_hist++){
        hGatetimesOpen[i_hist] = new TH1I(Form("hGatetimesOpen_Feb%d", i_hist), Form("Gate time on open gate of FEB #%d", i_hist), nbinsx, xlow, xup);
        hGatetimesOpen[i_hist]->GetXaxis()->SetTitle("Gate time [resol. 10 ms]");
        hGatetimesOpen[i_hist]->GetYaxis()->SetTitle("Number of Entries");

        hGateTypeMismatch[i_hist] = new TH2D(Form("hGateTypeMismatchFeb%d", i_hist), Form("Gate type mismatch w.r.t FEB #%d", i_hist), 14, -0.5, 13.5, 11, -5.5, 5.5);
        hGateTypeMismatch[i_hist]->GetXaxis()->SetNdivisions(15);
        hGateTypeMismatch[i_hist]->GetYaxis()->SetNdivisions(11);
        hGateTypeMismatch[i_hist]->GetYaxis()->SetTickLength(0.005);
        hGateTypeMismatch[i_hist]->GetXaxis()->SetTitle("FEB Number");
        hGateTypeMismatch[i_hist]->GetYaxis()->SetTitle("Type difference");
        hGateTypeMismatch[i_hist]->GetZaxis()->SetTitle("Number of Entries");

        hGateNumberMismatch[i_hist] = new TH2D(Form("hGateNumberMismatchFEB%d", i_hist), Form("Gate number differences w.r.t FEB #%d", i_hist), 14, -0.5, 13.5, 11, -5.5, 5.5);
        hGateNumberMismatch[i_hist]->GetXaxis()->SetNdivisions(15);
        hGateNumberMismatch[i_hist]->GetYaxis()->SetNdivisions(11);
        hGateNumberMismatch[i_hist]->GetYaxis()->SetTickLength(0.005);
        hGateNumberMismatch[i_hist]->GetXaxis()->SetTitle("FEB Number");
        hGateNumberMismatch[i_hist]->GetYaxis()->SetTitle("Number difference");
        hGateNumberMismatch[i_hist]->GetZaxis()->SetTitle("Number of Entries");

        hGatetimeMismatch[i_hist] = new TH2D(Form("hGatetimeMismatchFeb%d", i_hist), Form("Gate time from GTS differences w.r.t FEB #%d", i_hist), 14, -0.5, 13.5, 11, -5.5, 5.5);
        hGatetimeMismatch[i_hist]->GetXaxis()->SetNdivisions(15);
        hGatetimeMismatch[i_hist]->GetYaxis()->SetNdivisions(11);
        hGatetimeMismatch[i_hist]->GetYaxis()->SetTickLength(0.005);
        hGatetimeMismatch[i_hist]->GetXaxis()->SetTitle("FEB Number");
        hGatetimeMismatch[i_hist]->GetYaxis()->SetTitle("Time difference [resol. 10us]");
        hGatetimeMismatch[i_hist]->GetZaxis()->SetTitle("Number of Entries");

        hGatetimeOpengateMismatch[i_hist] = new TH2D(Form("hGatetimeOpengateMismatchFeb%d", i_hist), Form("Gate time on open gate differences w.r.t FEB #%d", i_hist), 14, -0.5, 13.5, 11, -5.5, 5.5);
        hGatetimeOpengateMismatch[i_hist]->GetXaxis()->SetNdivisions(15);
        hGatetimeOpengateMismatch[i_hist]->GetYaxis()->SetNdivisions(11);
        hGatetimeOpengateMismatch[i_hist]->GetYaxis()->SetTickLength(0.005);
        hGatetimeOpengateMismatch[i_hist]->GetXaxis()->SetTitle("FEB Number");
        hGatetimeOpengateMismatch[i_hist]->GetYaxis()->SetTitle("Time difference [resol. 10ms]");
        hGatetimeOpengateMismatch[i_hist]->GetZaxis()->SetTitle("Number of Entries");

        hGateNumberTrailerMismatch[i_hist] = new TH2D(Form("hGateNumberTrailerMismatchFEB%d", i_hist), Form("Gate number (trailer) differences w.r.t FEB #%d", i_hist), 14, -0.5, 13.5, 11, -5.5, 5.5);
        hGateNumberTrailerMismatch[i_hist]->GetXaxis()->SetNdivisions(15);
        hGateNumberTrailerMismatch[i_hist]->GetYaxis()->SetNdivisions(11);
        hGateNumberTrailerMismatch[i_hist]->GetYaxis()->SetTickLength(0.005);
        hGateNumberTrailerMismatch[i_hist]->GetXaxis()->SetTitle("FEB Number");
        hGateNumberTrailerMismatch[i_hist]->GetYaxis()->SetTitle("Number difference");
        hGateNumberTrailerMismatch[i_hist]->GetZaxis()->SetTitle("Number of Entries");

        hGatetimeClosedgateMismatch[i_hist] = new TH2D(Form("hGatetimeClosedgateMismatchFeb%d", i_hist), Form("Gate time on closed gate (trailer) differences w.r.t FEB #%d", i_hist), 14, -0.5, 13.5, 11, -5.5, 5.5);
        hGatetimeClosedgateMismatch[i_hist]->GetXaxis()->SetNdivisions(15);
        hGatetimeClosedgateMismatch[i_hist]->GetYaxis()->SetNdivisions(11);
        hGatetimeClosedgateMismatch[i_hist]->GetYaxis()->SetTickLength(0.005);
        hGatetimeClosedgateMismatch[i_hist]->GetXaxis()->SetTitle("FEB Number");
        hGatetimeClosedgateMismatch[i_hist]->GetYaxis()->SetTitle("Time difference [resol. 10ms]");
        hGatetimeClosedgateMismatch[i_hist]->GetZaxis()->SetTitle("Number of Entries");

        hGTSTagMismatch[i_hist] = new TH2D(Form("hGTSTagMismatchFeb%d", i_hist), Form("GTS tag of the first GTS trailer 1 (ID 4) differences w.r.t FEB #%d", i_hist), 14, -0.5, 13.5, 11, -5.5, 5.5);
        hGTSTagMismatch[i_hist]->GetXaxis()->SetNdivisions(15);
        hGTSTagMismatch[i_hist]->GetYaxis()->SetNdivisions(11);
        hGTSTagMismatch[i_hist]->GetYaxis()->SetTickLength(0.005);
        hGTSTagMismatch[i_hist]->GetXaxis()->SetTitle("FEB Number");
        hGTSTagMismatch[i_hist]->GetYaxis()->SetTitle("Tag difference");
        hGTSTagMismatch[i_hist]->GetZaxis()->SetTitle("Number of Entries");

        hGTSTimeMismatch[i_hist] = new TH2D(Form("hGTSTimeMismatchFeb%d", i_hist), Form("GTS time of the first GTS trailer 2 (ID 5) differences w.r.t FEB #%d", i_hist), 14, -0.5, 13.5, 11, -5.5, 5.5);
        hGTSTimeMismatch[i_hist]->GetXaxis()->SetNdivisions(15);
        hGTSTimeMismatch[i_hist]->GetYaxis()->SetNdivisions(11);
        hGTSTimeMismatch[i_hist]->GetYaxis()->SetTickLength(0.005);
        hGTSTimeMismatch[i_hist]->GetXaxis()->SetTitle("FEB Number");
        hGTSTimeMismatch[i_hist]->GetYaxis()->SetTitle("Time difference [resol. 10us]");
        hGTSTimeMismatch[i_hist]->GetZaxis()->SetTitle("Number of Entries");

        hLGAmplitude[i_hist] = new TH1I(Form("hLGAmplitudeFeb%d", i_hist), Form("LG Amplitude of FEB #%d", i_hist), 256, 0, 256);
        hLGAmplitude[i_hist]->GetXaxis()->SetTitle("Channel Number");
        hLGAmplitude[i_hist]->GetYaxis()->SetTitle("Number of Entries");
        hHGAmplitude[i_hist] = new TH1I(Form("hHGAmplitudeFeb%d", i_hist), Form("HG Amplitude of FEB #%d", i_hist), 256, 0, 256);
        hHGAmplitude[i_hist]->GetXaxis()->SetTitle("Channel Number");
        hHGAmplitude[i_hist]->GetYaxis()->SetTitle("Number of Entries");

    }

    hGateNumberMismatchFEB = new TH2D("hGateNumberMismatchinFEB", "Gate number differences between header and trailer in FEBs" , 14, -0.5, 13.5, 11, -5.5, 5.5);
    hGateNumberMismatchFEB->GetXaxis()->SetNdivisions(15);
    hGateNumberMismatchFEB->GetYaxis()->SetNdivisions(11);
    hGateNumberMismatchFEB->GetYaxis()->SetTickLength(0.005);
    hGateNumberMismatchFEB->GetXaxis()->SetTitle("FEB Number");
    hGateNumberMismatchFEB->GetYaxis()->SetTitle("Number difference (head. - tail.)");
    hGateNumberMismatchFEB->GetZaxis()->SetTitle("Number of Entries");

    hIncrementGTSTime = new TH2D("hIncrementGTSTime", "Increment between two GTS times (trailer 2, ID 5): GTS_time[n] - GTS_time[n-1]" , 14, -0.5, 13.5, 11, -5.5, 5.5);
    hIncrementGTSTime->GetXaxis()->SetNdivisions(15);
    hIncrementGTSTime->GetYaxis()->SetNdivisions(11);
    hIncrementGTSTime->GetYaxis()->SetTickLength(0.005);
    hIncrementGTSTime->GetXaxis()->SetTitle("FEB Number");
    hIncrementGTSTime->GetYaxis()->SetTitle("Increment");
    hIncrementGTSTime->GetZaxis()->SetTitle("Number of Entries");

    hNGTSperGate = new TH2D("hNGTSperGate", "Number of GTS trailer 1 (ID 4) per Gate" , 14, -0.5, 13.5, 10, -0.5, 9.5);
    hNGTSperGate->GetXaxis()->SetNdivisions(15);
    hNGTSperGate->GetYaxis()->SetNdivisions(10);
    hNGTSperGate->GetYaxis()->SetTickLength(0.005);
    hNGTSperGate->GetXaxis()->SetTitle("FEB Number");
    hNGTSperGate->GetYaxis()->SetTitle("Number of GTS trailer 1 per Gate");
    hNGTSperGate->GetZaxis()->SetTitle("Number of Entries");



    // Loop over the words until there is no more word in the file input_file
    while (fread(&word, sizeof(word), 1, read_fp)) {
        // Detect the empty words
        if(get_bits(word, 0, 31) == 0){
            are_prev_word_empty++;
            Nempty_words++;
            continue;
        }
        // If the previous word was empty but not this one, print out the empty
        // word error and set the is_checks_done to 1, such that no checks are run on this incomplete event
        if (are_prev_word_empty != 0 && get_bits(word, 0, 31) != 0){
            n_line++;
            fprintf(write_fp_words,"%d\t>> ERROR: Detected %d empty word(s) after Event Number: %d\n", n_line, are_prev_word_empty, events.back().OCB_header.nevent);
            std::array<int, 3> where_empty_word = {events.back().OCB_header.nevent, n_line, are_prev_word_empty};
            where_empty_words.push_back(where_empty_word);
            are_prev_word_empty = 0;
            is_checks_done = 1;
            Nerrors++;
        }

        // Read the Word ID
        int word_id = get_bits(word, 28, 4);

        // OCB header
        if (word_id == 8){
            // New OCB header == new event, therefore run the checks on the previous event
            if(READ_only == false && events.size() != 0 && is_checks_done == 0){
                Run_checks(events);
            } else {
                is_checks_done = 0;
            }

            if (to_read_Nevent !=0 && Nevent >= to_read_Nevent) {
                is_Nevent_limit = 1;
                break;
            }

            //Print out the progress
            if (Nevent % 1000 == 0) {
                if (Nevent != 0){
                    fflush(stdout);
                }
                printf("The program is running, %d events were read\n", Nevent);
            }

            // Add an element to the events vector and read the word
            Add_struct_event(events);
            //Increase the counter of event
            Nevent++;
            events[Nevent-1].OCB_header.gate_type = get_bits(word, 25, 3);
            events[Nevent-1].OCB_header.gate_tag = get_bits(word, 23, 2);
            events[Nevent-1].OCB_header.nevent = get_bits(word, 0, 23);
            events[Nevent-1].is_OCB_header = 1;

            n_line++;
            fprintf(write_fp_words,"%d\tID  8: OCB header;\tGate type: %d;\tGate tag: %d;\tEvent number:%d\n", n_line, events[Nevent-1].OCB_header.gate_type, events[Nevent-1].OCB_header.gate_tag, events[Nevent-1].OCB_header.nevent);

        // Gate headers
        } else if (word_id == 0){
            nSlot_FEB = get_bits(word, 20, 4); //lower 4 bits of the gate number is the slot number;
            events[Nevent-1].is_FEBs_packet[nSlot_FEB] = 1;
            // Gate header A
            if (get_bits(word, 19, 1) == 0){
                events[Nevent-1].FEBs_packet[nSlot_FEB].gate_header.A_board_id = get_bits(word, 20, 8);
                events[Nevent-1].FEBs_packet[nSlot_FEB].gate_header.A_gate_type = get_bits(word, 16, 3);
                events[Nevent-1].FEBs_packet[nSlot_FEB].gate_header.A_gate_number = get_bits(word, 0, 16);
                events[Nevent-1].FEBs_packet[nSlot_FEB].is_gate_headerA = 1;
                n_line++;
                fprintf(write_fp_words,"%d\tID  0: Gate header A;\tBoard ID: %d (FEB %d);\tGate type: %d;\tGate number: %d\n", n_line, events[Nevent-1].FEBs_packet[nSlot_FEB].gate_header.A_board_id, get_bits(word, 20, 4), events[Nevent-1].FEBs_packet[nSlot_FEB].gate_header.A_gate_type, events[Nevent-1].FEBs_packet[nSlot_FEB].gate_header.A_gate_number);
            // Gate header B
            }else if (get_bits(word, 19, 1) == 1){
                events[Nevent-1].FEBs_packet[nSlot_FEB].gate_header.B_board_id = get_bits(word, 20, 8);
                events[Nevent-1].FEBs_packet[nSlot_FEB].gate_header.B_gate_time = get_bits(word, 0, 11);
                events[Nevent-1].FEBs_packet[nSlot_FEB].is_gate_headerB = 1;
                n_line++;
                fprintf(write_fp_words,"%d\tID  0: Gate header B;\tBoard ID: %d (FEB %d);\tGate time: %d;\n", n_line, events[Nevent-1].FEBs_packet[nSlot_FEB].gate_header.B_board_id, get_bits(word, 20, 4), events[Nevent-1].FEBs_packet[nSlot_FEB].gate_header.B_gate_time);
            }
            prev_gate_word = 0;
        }
        // Gate time
        else if (word_id == 7 && prev_gate_word == 0){
            events[Nevent-1].FEBs_packet[nSlot_FEB].gate_header.gate_time = get_bits(word, 0, 28);
            events[Nevent-1].FEBs_packet[nSlot_FEB].is_gate_header_time = 1;
            n_line++;
            fprintf(write_fp_words,"%d\tID  7: Gate time on open gate\t\tGate time: %d;\n", n_line, events[Nevent-1].FEBs_packet[nSlot_FEB].gate_header.gate_time);

        //HOLD time
        } else if (word_id == 11){
            if (get_bits(word, 19, 1)==0){
                events[Nevent-1].FEBs_packet[nSlot_FEB].hold_time.start_board_id = get_bits(word, 20, 8);
                events[Nevent-1].FEBs_packet[nSlot_FEB].hold_time.start_time = get_bits(word, 0, 28);
                events[Nevent-1].FEBs_packet[nSlot_FEB].is_hold_start_time = 1;
                n_line++;
                fprintf(write_fp_words,"%d\tID  11: HOLD time;\tBoard id: %d (FEB %d);\tHold Start time %d\n", n_line, events[Nevent-1].FEBs_packet[nSlot_FEB].hold_time.start_board_id, get_bits(word, 20, 4), events[Nevent-1].FEBs_packet[nSlot_FEB].hold_time.start_time);
            }else if (get_bits(word, 19, 1)==1){
                events[Nevent-1].FEBs_packet[nSlot_FEB].hold_time.stop_board_id = get_bits(word, 20, 8);
                events[Nevent-1].FEBs_packet[nSlot_FEB].hold_time.stop_time = get_bits(word, 0, 28);
                events[Nevent-1].FEBs_packet[nSlot_FEB].is_hold_stop_time = 1;
                n_line++;
                fprintf(write_fp_words,"%d\tID  11: HOLD time;\tBoard id: %d (FEB %d);\tHold Stop time %d\n", n_line, events[Nevent-1].FEBs_packet[nSlot_FEB].hold_time.stop_board_id, get_bits(word, 20, 4), events[Nevent-1].FEBs_packet[nSlot_FEB].hold_time.stop_time);
            }

        //GTS headers
        } else if (word_id == 1){
            Add_struct_GTS_header(events[Nevent-1].FEBs_packet[nSlot_FEB]); //.GTS.headers.is_inGate ckecked in the add function directlz\y
            events[Nevent-1].FEBs_packet[nSlot_FEB].GTS_headers.back().GTS_tag = get_bits(word, 0, 28);
            events[Nevent-1].FEBs_packet[nSlot_FEB].GTS_IDLog.push_back(word_id);
            n_line++;
            fprintf(write_fp_words,"%d\tID  1: GTS header;\tGTS tag: %d\n", n_line, events[Nevent-1].FEBs_packet[nSlot_FEB].GTS_headers.back().GTS_tag);

        //Hit time
        } else if (VERBOSE == true && word_id == 2){
            n_line++;
            fprintf(write_fp_words,"%d\tID  2: Hit time;\tChannel ID: %d\tHit ID: %d\tTag ID: %d\tEdge: %d\tHit time:%d\n", n_line, get_bits(word, 20, 8), get_bits(word, 17, 3), get_bits(word, 15, 2), get_bits(word, 14, 1), get_bits(word, 0, 13));

        //Hit Amplitude
        } else if (word_id == 3){
             Gain = get_bits(word, 12, 3);
             Channel_id = get_bits(word, 20, 8);
            if (Gain == 2){
                hHGAmplitude[nSlot_FEB]->Fill(Channel_id);
            } else if (Gain == 3){
                hLGAmplitude[nSlot_FEB]->Fill(Channel_id);
            }
            if (VERBOSE == true){
              n_line++;
              fprintf(write_fp_words,"%d\tID  3: Hit Amplitude;\tChannel ID: %d\tHit ID: %d\tTag ID: %d\tAmplitude ID: %d\tAmplitude measurement:%d\n", n_line, get_bits(word, 20, 8), get_bits(word, 17, 3), get_bits(word, 15, 2), get_bits(word, 12, 3), get_bits(word, 0, 12));
            }

        //GTS trailers
        } else if (word_id == 4){
            Add_struct_GTS_trailers(events[Nevent-1].FEBs_packet[nSlot_FEB]);
            events[Nevent-1].FEBs_packet[nSlot_FEB].GTS_trailers.back().GTS_tag = get_bits(word, 0, 28);
            if (events[Nevent-1].FEBs_packet[nSlot_FEB].is_gate_trailer == 0){ // No gate trailer yet = the GTS is in the gate
                events[Nevent-1].FEBs_packet[nSlot_FEB].GTS_trailers.back().are_inGate[0] = 1;
            }
            events[Nevent-1].FEBs_packet[nSlot_FEB].GTS_IDLog.push_back(word_id);
            n_line++;
            fprintf(write_fp_words,"%d\tID  4: GTS trailer1;\tGTS tag: %d\n", n_line, events[Nevent-1].FEBs_packet[nSlot_FEB].GTS_trailers.back().GTS_tag);

        } else if (word_id == 5){
            if (get_bits(prev_word, 28, 4) != 4){ // If the previous word was not the GTS trailer 1 (id 4), add a new struct_GTS_trailers
                Add_struct_GTS_trailers(events[Nevent-1].FEBs_packet[nSlot_FEB]);
            }
            events[Nevent-1].FEBs_packet[nSlot_FEB].GTS_trailers.back().Data = get_bits(word, 27, 1);
            events[Nevent-1].FEBs_packet[nSlot_FEB].GTS_trailers.back().GTS_time = get_bits(word, 0, 20);
            if (events[Nevent-1].FEBs_packet[nSlot_FEB].is_gate_trailer == 0){ // No gate trailer yet = the GTS is in the gate
                events[Nevent-1].FEBs_packet[nSlot_FEB].GTS_trailers.back().are_inGate[1] = 1;
            }
            events[Nevent-1].FEBs_packet[nSlot_FEB].GTS_IDLog.push_back(word_id);
            n_line++;
            fprintf(write_fp_words,"%d\tID  5: GTS trailer2;\tData: %d\t\tGTS time: %d\n", n_line, events[Nevent-1].FEBs_packet[nSlot_FEB].GTS_trailers.back().Data, events[Nevent-1].FEBs_packet[nSlot_FEB].GTS_trailers.back().GTS_time);

        //Gate trailer
        } else if (word_id == 6){
            events[Nevent-1].FEBs_packet[nSlot_FEB].gate_trailer.board_id = get_bits(word, 20, 8);
            events[Nevent-1].FEBs_packet[nSlot_FEB].gate_trailer.gate_type = get_bits(word, 16, 3);
            events[Nevent-1].FEBs_packet[nSlot_FEB].gate_trailer.gate_number = get_bits(word, 0, 16);
            events[Nevent-1].FEBs_packet[nSlot_FEB].is_gate_trailer = 1;
            prev_gate_word = 1;
            n_line++;
            fprintf(write_fp_words,"%d\tID  6: Gate trailer;\tBoard ID: %d (FEB %d)\tGate type: %d;\tGate number: %d\n", n_line, events[Nevent-1].FEBs_packet[nSlot_FEB].gate_trailer.board_id, get_bits(word, 20, 4), events[Nevent-1].FEBs_packet[nSlot_FEB].gate_trailer.gate_type, events[Nevent-1].FEBs_packet[nSlot_FEB].gate_trailer.gate_number);

        } else if (word_id == 7 && prev_gate_word == 1){
            events[Nevent-1].FEBs_packet[nSlot_FEB].gate_trailer.gate_time = get_bits(word, 0, 28);
            events[Nevent-1].FEBs_packet[nSlot_FEB].is_gate_trailer_time = 1;
            n_line++;
            fprintf(write_fp_words,"%d\tID  7: Gate time on closed gate\t\tGate time: %d;\n", n_line, events[Nevent-1].FEBs_packet[nSlot_FEB].gate_trailer.gate_time);

        //FEB trailer
        } else if (word_id == 13){
            events[Nevent-1].FEBs_packet[nSlot_FEB].FEB_trailer.Event_done_timeout = get_bits(word, 18, 1);
            events[Nevent-1].FEBs_packet[nSlot_FEB].FEB_trailer.D1_FIFO_full = get_bits(word, 17, 1);
            events[Nevent-1].FEBs_packet[nSlot_FEB].FEB_trailer.D0_FIFO_full = get_bits(word, 16, 1);
            events[Nevent-1].FEBs_packet[nSlot_FEB].FEB_trailer.N_decoder_error = get_bits(word, 0, 16);
            events[Nevent-1].FEBs_packet[nSlot_FEB].is_FEB_trailer = 1;
            n_line++;
            fprintf(write_fp_words,"%d\tID  13: FEB trailer;\t Event done timeout: %d;\tD1 FIFO full: %d;\tD0 FIFO full: %d;\tNumber of decoder errors: %d\n", n_line, events[Nevent-1].FEBs_packet[nSlot_FEB].FEB_trailer.Event_done_timeout, events[Nevent-1].FEBs_packet[nSlot_FEB].FEB_trailer.D1_FIFO_full, events[Nevent-1].FEBs_packet[nSlot_FEB].FEB_trailer.D0_FIFO_full, events[Nevent-1].FEBs_packet[nSlot_FEB].FEB_trailer.N_decoder_error);

        //OCB trailer
        } else if (word_id == 9){
            events[Nevent-1].OCB_trailer.gate_open_timeout = get_bits(word, 15, 1);
            events[Nevent-1].OCB_trailer.gate_close_timeout = get_bits(word, 14, 1);
            for (int i = 0; i < N_FEB; i++){
                events[Nevent-1].OCB_trailer.FEB_error[i] = get_bits(word, i, 1);
            }
            events[Nevent-1].is_OCB_trailer = 1;
            n_line++;
            fprintf(write_fp_words,"%d\tID  9: OCB trailer\n", n_line);
        }
        // Keep track of the previous word
        prev_word = word;
    }
    if (READ_only == true){
        printf("The program has ended, All the words (%d events) were read but not checked (READ only mode)\n", Nevent);
        // close files
        fclose(read_fp);
        fclose(write_fp_words);
        fclose(write_fp_errors);
        file->Close();
        return 0;
    }
    // Run the checks on the last event if it was not done (if it was not break by to_read_Nevent limit) and if the last event is complete (it has an OCB trailer)
    if (is_Nevent_limit == 0){
        if (events.back().is_OCB_trailer == 1){
            Run_checks(events);
        }
        printf("The program has ended, All the words (%d events) were read\n", Nevent);
    }else {
        printf("The program has ended, limit of the number of events to read reached, %d events read\n", Nevent);
    }
    printf("checks done/n");
    // Print out the probelmatic events
    // if (problematic_event.size() != 0){
    //     for (int i = 0; i < problematic_event.size(); i++){
    //         // Find the index of the event problematic_event[i]:
    //         int i_probEvent = problematic_event[i]-1;
    //         if (events[i_probEvent].OCB_header.nevent != problematic_event[i]){
    //           while (events[i_probEvent].OCB_header.nevent != problematic_event[i]){
    //               i_probEvent--;
    //           }
    //         }
    //         // Print the event
    //         printEvent(events[i_probEvent]);
    //     }
    // }

    //Print out the Summary
    if(Nerrors == 0) {
        fprintf(write_fp_errors, "No Errors for all! :)\n\n");
        fprintf(write_fp_errors, "No Errors for all! :)\n\n");
        fprintf(write_fp_errors, "No Errors for all! :)\n\n");
    }

    fprintf(write_fp_errors, "\n");
    fprintf(write_fp_errors, "====================== SUMMARY of the ERRORS ======================================================\n");
    fprintf(write_fp_errors, "| Total number of events: %d                      \n", Nevent);
    fprintf(write_fp_errors, "| Total number of events checked (full event): %d \n", Nevent_checked);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of errors: %d                      \n", Nerrors);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of missing Gate header words: %d   \n", N_missing_gate_headerA+N_missing_gate_headerB+N_missing_gate_header_time);
    fprintf(write_fp_errors, "| Total number of missing GTS words: %d           \n", N_missing_GTS_trailer1+N_missing_GTS_trailer2+N_missing_GTS_header);
    fprintf(write_fp_errors, "| Total number of missing Gate Trailer : %d       \n", N_missing_gate_trailer + N_missing_gate_trailer_time);
    fprintf(write_fp_errors, "| Total number of missing FEB trailer: %d         \n", Nmissing_FEB_trailer);
    fprintf(write_fp_errors, "| Total number of missing OCB trailer: %d           \n", where_missing_OCB_trailer.size());
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total FEB Trailer errors: %d                    \n", N_event_done_timeout + N_D1_FIFO_full + N_D0_FIFO_full);
    fprintf(write_fp_errors, "| Total OCB Trailer errors: %d                    \n", where_gate_open_timeout.size()+ where_gate_close_timeout.size() + N_FEB_error);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total Gate Tag vs Gate Number mismatch: %d \n", N_gate_TagNumb_error);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of information mismatch between Gate header and trailer: %d \n", N_MmFEB_gate_number);
    fprintf(write_fp_errors, "| Total number of error about GTS word information w.r.t a single FEB: %d \n", N_W_NGTSperGate+N_W_GTS1+N_W_GTS_tag+N_W_GTS_time+N_MmFEB_GTS_tag);
    fprintf(write_fp_errors, "| Total number of error about GTS word information mismatch between FEB: %d \n", N_Mm_GTS_tag+N_Mm_GTS_time);
    fprintf(write_fp_errors, "| Total number of event with Gate headers information mismatch (w.r.t. FEB 0): %d \n", N_Mm_gate_type+N_Mm_gate_number+N_Mm_gate_time+N_Mm_gate_timeOpenGate);
    fprintf(write_fp_errors, "| Total number of event with Gate trailer information mismatch (w.r.t. FEB 0): %d \n", N_Mm_gate_trailer_number+N_Mm_gate_timeClosedGate);
    fprintf(write_fp_errors, "| Total number of event with a wrong Gate trailer type: %d \n", N_W_gate_trailer_type);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Number of skipped events: %d                    \n", Nskip_event); // Number of number of event skipped
    fprintf(write_fp_errors, "| Number of empty word: %d                        \n", Nempty_words); // Number of number of event skipped
    fprintf(write_fp_errors, "===================================================================================================\n");
    fprintf(write_fp_errors, "\n");
    fprintf(write_fp_errors, "\n");
    fprintf(write_fp_errors, "\n");
    fprintf(write_fp_errors, "===================================================================================================\n");
    fprintf(write_fp_errors, "====================== Detailed Errors ============================================================\n");
    fprintf(write_fp_errors, "===================================================================================================\n");
    fprintf(write_fp_errors, "\n");
    fprintf(write_fp_errors, "====================== Gate Errors ================================================================\n");
    fprintf(write_fp_errors, "| Total number of missing Gate header words: %d   \n", N_missing_gate_headerA+N_missing_gate_headerB+N_missing_gate_header_time);
    fprintf(write_fp_errors, "|     Number of missing Gate header type A: %d  \n", N_missing_gate_headerA);
    printArrayofVector(write_fp_errors, where_missing_gate_headerA);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "|     Number of missing Gate header type B: %d  \n", N_missing_gate_headerB);
    printArrayofVector(write_fp_errors, where_missing_gate_headerB);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "|     Number of missing Gate time: %d  \n", N_missing_gate_header_time);
    printArrayofVector(write_fp_errors, where_missing_gate_header_time);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of missing Gate trailer words: %d   \n", N_missing_gate_trailer+N_missing_gate_trailer_time);
    fprintf(write_fp_errors, "|     Number of missing Gate trailer: %d          \n", N_missing_gate_trailer);
    printArrayofVector(write_fp_errors, where_missing_gate_trailer);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "|     Number of missing Gate trailer time: %d      \n", N_missing_gate_trailer_time);
    printArrayofVector(write_fp_errors, where_missing_gate_trailer_time);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of event with Gate number mismatch between header (header A, ID 0) and trailer (id 6): %d \n", N_MmFEB_gate_number);
    printArrayofVector(write_fp_errors, where_MmFEB_gate_number);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of event with Gate type (header A, ID 0) mismatch (w.r.t. FEB 0): %d \n", N_Mm_gate_type);
    printArrayofVector(write_fp_errors, where_Mm_gate_type);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of event with Gate number (header A, ID 0) mismatch (w.r.t. FEB 0): %d \n", N_Mm_gate_number);
    printArrayofVector(write_fp_errors, where_Mm_gate_number);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of event with Gate Time from GTS (header B, ID 0) mismatch (w.r.t. FEB 0): %d \n", N_Mm_gate_time);
    printArrayofVector(write_fp_errors, where_Mm_gate_time);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of event with Gate Time on open gate (ID 7) mismatch (w.r.t. FEB 0): %d \n", N_Mm_gate_timeOpenGate);
    printArrayofVector(write_fp_errors, where_Mm_gate_timeOpenGate);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of event with Gate number (trailer, ID 6) mismatch (w.r.t. FEB 0): %d \n", N_Mm_gate_trailer_number);
    printArrayofVector(write_fp_errors, where_Mm_gate_trailer_number);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of event with a wrong Gate type (trailer, ID 6): %d \n", N_W_gate_trailer_type);
    printArrayofVector(write_fp_errors, where_W_gate_trailer_type);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of event with Gate Time on closed gate (ID 7) mismatch (w.r.t. FEB 0): %d \n", N_Mm_gate_timeClosedGate);
    printArrayofVector(write_fp_errors, where_Mm_gate_timeClosedGate);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of mismatch between Gate Tag and Gate Number(header A ID 0): %d                                                 \n", N_gate_TagNumb_error);
    printArrayofVector(write_fp_errors, where_gate_TagNumb_error);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "===================================================================================================\n");
    fprintf(write_fp_errors, "\n");
    fprintf(write_fp_errors, "\n");
    fprintf(write_fp_errors, "====================== GTS Errors =================================================================\n");
    fprintf(write_fp_errors, "| Total number of missing GTS words: %d           \n", N_missing_GTS_trailer1+N_missing_GTS_trailer2+N_missing_GTS_header);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "|     Number of missing GTS trailer 1: %d         \n", N_missing_GTS_trailer1);
    printArrayofVector(write_fp_errors, where_missing_GTS_trailer1);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "|     Number of missing GTS trailer 2: %d         \n", N_missing_GTS_trailer2);
    printArrayofVector(write_fp_errors, where_missing_GTS_trailer2);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "|     Number of missing GTS header: %d            \n", N_missing_GTS_header);
    printArrayofVector(write_fp_errors, where_missing_GTS_header);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of event with GTS Tag (trailer 1, ID 4) mismatch (w.r.t. the first GTS tag recorded): %d\n| (Note: The FEB associated to the error is the first one to have shown the error) \n", N_Mm_GTS_tag);
    printArrayofVector(write_fp_errors, where_Mm_GTS_tag);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of event with GTS Time (trailer 2, ID 5) mismatch (w.r.t. the first GTS time recorded): %d \n| (Note: The FEB associated to the error is the first one to have shown the error) \n", N_Mm_GTS_time);
    printArrayofVector(write_fp_errors, where_Mm_GTS_time);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of mismatch GTS tag betweem header (ID 1) and trailer1 (ID 4): %d \n", N_MmFEB_GTS_tag);
    printArrayofVector(write_fp_errors, where_MmFEB_GTS_tag);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of Gate with wrong first GTS words: %d \n", N_W_GTS1);
    printArrayofVector(write_fp_errors, where_W_GTS1);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of Gate with wrong number of GTS words inside: %d \n| (Note: The FEB associated to the error is the first one to have shown the error) \n", N_W_NGTSperGate);
    printArrayofVector(write_fp_errors, where_W_NGTSperGate);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of Gate with wrong GTS tag increment: %d \n", N_W_GTS_tag);
    printArrayofVector(write_fp_errors, where_W_GTS_tag);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of Gate with wrong GTS time increment: %d \n| (Note: Only one error per FEB, per event is counted) \n", N_W_GTS_time);
    printArrayofVector(write_fp_errors, where_W_GTS_time);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "===================================================================================================\n");
    fprintf(write_fp_errors, "\n");
    fprintf(write_fp_errors, "\n");
    fprintf(write_fp_errors, "====================== FEB Errors =================================================================\n");
    fprintf(write_fp_errors, "|     Number of missing FEB trailer: %d           \n", Nmissing_FEB_trailer);
    fprintf(write_fp_errors, "| Per slots: ");
    printArray(write_fp_errors, slot_missing_FEB_trailer);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Number of Event done timeout errors: %d         \n", N_event_done_timeout);
    fprintf(write_fp_errors, "| Per slots: ");
    printArray(write_fp_errors, slot_EventdoneTimeout);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Number of D1 FIFO full errors: %d               \n", N_D1_FIFO_full);
    fprintf(write_fp_errors, "| Per slots: ");
    printArray(write_fp_errors, slot_D1_FIFO_full);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Number of D0 FIFO full errors: %d               \n", N_D0_FIFO_full);
    fprintf(write_fp_errors, "| Per slots: ");
    printArray(write_fp_errors, slot_D0_FIFO_full);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Total number of decoder errors: %d              \n", N_decoder_errors);
    fprintf(write_fp_errors, "| Per slots: ");
    printArray(write_fp_errors, slot_Numb_decoder_errors);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "===================================================================================================\n");
    fprintf(write_fp_errors, "\n");
    fprintf(write_fp_errors, "\n");
    fprintf(write_fp_errors, "====================== OCB Errors =================================================================\n");
    fprintf(write_fp_errors, "|     Number of missing OCB trailer: %d           \n", where_missing_OCB_trailer.size());
    printVectorofVector(write_fp_errors, where_missing_OCB_trailer);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "|     Number of Gate open timeout: %d             \n", where_gate_open_timeout.size());
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "|     Number of Gate close error: %d              \n", where_gate_close_timeout.size());
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "|     Number of FEB data packet i error: %d       \n", N_FEB_error);
    fprintf(write_fp_errors, "| Per FEB: ");
    printArray(write_fp_errors, slot_FEB_error);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "===================================================================================================\n");
    fprintf(write_fp_errors, "\n");
    fprintf(write_fp_errors, "\n");
    fprintf(write_fp_errors, "====================== General Errors =============================================================\n");
    fprintf(write_fp_errors, "| Number of jump in events number : %d            \n", Nskip_event); // Number of event number jump
    fprintf(write_fp_errors, "| where (jump: [from, to, line]):");
    printVectorofArray3(write_fp_errors, where_skip_events);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "| Number of empty word: %d                        \n", Nempty_words); // Number of empty word
    fprintf(write_fp_errors, "| where ([after/during event number, line, number of empty word in a row]):");
    printVectorofArray3(write_fp_errors, where_empty_words);
    fprintf(write_fp_errors, "|\n");
    fprintf(write_fp_errors, "===================================================================================================\n");
    fprintf(write_fp_errors, "\n");
    fprintf(write_fp_errors, "\n");
    fprintf(write_fp_errors, "\n");


    //Print the histos in .root file and .pdf file
    c1->Print(Form("%s[", output_pdf));     // Opens output file for multipage
    c1->Clear();

    //For hGatetimesOpen, set the correct range and print it
    c1->Divide(5, 3);
    for(int j_feb = 0; j_feb < N_FEB; j_feb++){
        c1->cd(j_feb+1);
        int LastBin = hGatetimesOpen[j_feb]->FindLastBinAbove(0, 1, 1, -1);
        int FirstBin = hGatetimesOpen[j_feb]->FindFirstBinAbove(0, 1, 1, -1);
        int margin = (LastBin - FirstBin)/10;
        if (FirstBin > margin){
            FirstBin -= margin;
        }else{
            FirstBin = 0;
        }
        LastBin += margin;
        hGatetimesOpen[j_feb]->GetXaxis()->SetRange(FirstBin, LastBin);
        hGatetimesOpen[j_feb]->Draw();
    }
    c1->Print(output_pdf);
    c1->Write();
    c1->Clear();

    c1->Divide(5, 3);
    gPad->SetRightMargin(0.18); //to be able to see the z axis title
    for(int j_feb = 0; j_feb < N_FEB; j_feb++){
        c1->cd(j_feb+1);
        hGateTypeMismatch[j_feb]->GetZaxis()->SetTitleOffset(1.3);
        hGateTypeMismatch[j_feb]->SetStats(0);
        hGateTypeMismatch[j_feb]->Draw("Text COLZ");
        TPaveText *text_hGateTypeMismatch = new TPaveText(0.55,0.84,0.8,0.89, "tblrNDC");
        text_hGateTypeMismatch->AddText(Form("Total Nbre of Entries per FEB: %d", Nevent_checked));
        // text_hGateTypeMismatch->SetTextSize(0.005);
        text_hGateTypeMismatch->SetBorderSize(1);
        text_hGateTypeMismatch->Draw("same");

    }
    c1->Print(output_pdf);
    c1->Write();
    c1->Clear();

    c1->Divide(5, 3);
    for(int j_feb = 0; j_feb < N_FEB; j_feb++){
        c1->cd(j_feb+1);
        hGateNumberMismatch[j_feb]->GetZaxis()->SetTitleOffset(1.3);
        hGateNumberMismatch[j_feb]->SetStats(0);
        hGateNumberMismatch[j_feb]->Draw("Text COLZ");
        TPaveText *text_hGateNumberMismatch = new TPaveText(0.55,0.84,0.8,0.89, "tblrNDC");
        text_hGateNumberMismatch->AddText(Form("Total Nbre of Entries per FEB: %d", Nevent_checked));
        // text_hGateNumberMismatch->SetTextSize(0.005);
        text_hGateNumberMismatch->SetBorderSize(1);
        text_hGateNumberMismatch->Draw("same");

    }
    c1->Print(output_pdf);
    c1->Write();
    c1->Clear();

    c1->Divide(5, 3);
    for(int j_feb = 0; j_feb < N_FEB; j_feb++){
        c1->cd(j_feb+1);
        hGatetimeMismatch[j_feb]->GetZaxis()->SetTitleOffset(1.3);
        hGatetimeMismatch[j_feb]->SetStats(0);
        hGatetimeMismatch[j_feb]->Draw("Text COLZ");
        TPaveText *text_hGatetimeMismatch = new TPaveText(0.55,0.84,0.8,0.89, "tblrNDC");
        text_hGatetimeMismatch->AddText(Form("Total Nbre of Entries per FEB: %d", Nevent_checked));
        // text_hGatetimeMismatch->SetTextSize(0.005);
        text_hGatetimeMismatch->SetBorderSize(1);
        text_hGatetimeMismatch->Draw("same");
    }
    c1->Print(output_pdf);
    c1->Write();
    c1->Clear();

    c1->Divide(5, 3);
    for(int j_feb = 0; j_feb < N_FEB; j_feb++){
        c1->cd(j_feb+1);
        hGatetimeOpengateMismatch[j_feb]->GetZaxis()->SetTitleOffset(1.3);
        hGatetimeOpengateMismatch[j_feb]->SetStats(0);
        hGatetimeOpengateMismatch[j_feb]->Draw("Text COLZ");
        TPaveText *text_hGatetimeOpengateMismatch = new TPaveText(0.55,0.84,0.8,0.89, "tblrNDC");
        text_hGatetimeOpengateMismatch->AddText(Form("Total Nbre of Entries per FEB: %d", Nevent_checked));
        text_hGatetimeOpengateMismatch->SetBorderSize(1);
        text_hGatetimeOpengateMismatch->Draw("same");
    }
    c1->Print(output_pdf);
    c1->Write();
    c1->Clear();

    c1->Divide(5, 3);
    for(int j_feb = 0; j_feb < N_FEB; j_feb++){
        c1->cd(j_feb+1);
        hGTSTagMismatch[j_feb]->GetZaxis()->SetTitleOffset(1.3);
        hGTSTagMismatch[j_feb]->SetStats(0);
        hGTSTagMismatch[j_feb]->Draw("Text COLZ");
        TPaveText *text_hGTSTagMismatch = new TPaveText(0.55,0.84,0.8,0.89, "tblrNDC");
        text_hGTSTagMismatch->AddText(Form("Total Nbre of Entries per FEB: %d", Nevent_checked));
        text_hGTSTagMismatch->SetBorderSize(1);
        text_hGTSTagMismatch->Draw("same");
    }
    c1->Print(output_pdf);
    c1->Write();
    c1->Clear();

    c1->Divide(5, 3);
    for(int j_feb = 0; j_feb < N_FEB; j_feb++){
        c1->cd(j_feb+1);
        hGTSTimeMismatch[j_feb]->GetZaxis()->SetTitleOffset(1.3);
        hGTSTimeMismatch[j_feb]->SetStats(0);
        hGTSTimeMismatch[j_feb]->Draw("Text COLZ");
        TPaveText *text_hGTSTimeMismatch = new TPaveText(0.55,0.84,0.8,0.89, "tblrNDC");
        text_hGTSTimeMismatch->AddText(Form("Total Nbre of Entries per FEB: %d", Nevent_checked));
        text_hGTSTimeMismatch->SetBorderSize(1);
        text_hGTSTimeMismatch->Draw("same");
    }
    c1->Print(output_pdf);
    c1->Write();
    c1->Clear();

    // c1->Divide(5, 3);
    // for(int j_feb = 0; j_feb < N_FEB; j_feb++){
    //     c1->cd(j_feb+1);
    //     [j_feb]->GetZaxis()->SetTitleOffset(1.3);
    //     [j_feb]->SetStats(0);
    //     [j_feb]->Draw("Text COLZ");
    //     TPaveText *text_ = new TPaveText(0.55,0.84,0.8,0.89, "tblrNDC");
    //     text_->AddText(Form("Total Nbre of Entries per FEB: %d", Nevent_checked));
    //     text_->SetBorderSize(1);
    //     text_->Draw("same");
    // }
    // c1->Print(output_pdf);
    // c1->Write();
    // c1->Clear();

    c1->Divide(5, 3);
    for(int j_feb = 0; j_feb < N_FEB; j_feb++){
        c1->cd(j_feb+1);
        hGateNumberTrailerMismatch[j_feb]->GetZaxis()->SetTitleOffset(1.3);
        hGateNumberTrailerMismatch[j_feb]->SetStats(0);
        hGateNumberTrailerMismatch[j_feb]->Draw("Text COLZ");
        TPaveText *text_hGateNumberTrailerMismatch = new TPaveText(0.55,0.84,0.8,0.89, "tblrNDC");
        text_hGateNumberTrailerMismatch->AddText(Form("Total Nbre of Entries per FEB: %d", Nevent_checked));
        text_hGateNumberTrailerMismatch->SetBorderSize(1);
        text_hGateNumberTrailerMismatch->Draw("same");
    }
    c1->Print(output_pdf);
    c1->Write();
    c1->Clear();

    c1->Divide(5, 3);
    for(int j_feb = 0; j_feb < N_FEB; j_feb++){
        c1->cd(j_feb+1);
        hGatetimeClosedgateMismatch[j_feb]->GetZaxis()->SetTitleOffset(1.3);
        hGatetimeClosedgateMismatch[j_feb]->SetStats(0);
        hGatetimeClosedgateMismatch[j_feb]->Draw("Text COLZ");
        TPaveText *text_hGatetimeClosedgateMismatch = new TPaveText(0.55,0.84,0.8,0.89, "tblrNDC");
        text_hGatetimeClosedgateMismatch->AddText(Form("Total Nbre of Entries per FEB: %d", Nevent_checked));
        text_hGatetimeClosedgateMismatch->SetBorderSize(1);
        text_hGatetimeClosedgateMismatch->Draw("same");
    }
    c1->Print(output_pdf);
    c1->Write();
    c1->Clear();

    // c1->Divide(5, 3);
    // for(int j_feb = 0; j_feb < N_FEB; j_feb++){
    //     c1->cd(j_feb+1);
    //     [j_feb]->GetZaxis()->SetTitleOffset(1.3);
    //     [j_feb]->SetStats(0);
    //     [j_feb]->Draw("Text COLZ");
    //     TPaveText *text_ = new TPaveText(0.55,0.84,0.8,0.89, "tblrNDC");
    //     text_->AddText(Form("Total Nbre of Entries per FEB: %d", Nevent_checked));
    //     text_->SetBorderSize(1);
    //     text_->Draw("same");
    // }
    // c1->Print(output_pdf);
    // c1->Write();
    // c1->Clear();



    hGateNumberMismatchFEB->GetZaxis()->SetTitleOffset(1.3);
    hGateNumberMismatchFEB->SetStats(0);
    hGateNumberMismatchFEB->Draw("Text COLZ");
    TPaveText *text_hGateNumberMismatchFEB = new TPaveText(0.55,0.84,0.8,0.89, "tblrNDC");
    text_hGateNumberMismatchFEB->AddText(Form("Total Nbre of Entries per FEB: %d", Nevent_checked));
    text_hGateNumberMismatchFEB->SetBorderSize(1);
    text_hGateNumberMismatchFEB->Draw("same");

    c1->Print(output_pdf);
    c1->Write();
    c1->Clear();

    hNGTSperGate->GetZaxis()->SetTitleOffset(1.3);
    hNGTSperGate->SetStats(0);
    hNGTSperGate->Draw("Text COLZ");
    TPaveText *text_hNGTSperGate = new TPaveText(0.55,0.84,0.8,0.89, "tblrNDC");
    text_hNGTSperGate->AddText(Form("Total Nbre of Entries per FEB: %d", Nevent_checked));
    text_hNGTSperGate->SetBorderSize(1);
    text_hNGTSperGate->Draw("same");

    c1->Print(output_pdf);
    c1->Write();
    c1->Clear();
    hIncrementGTSTime->GetZaxis()->SetTitleOffset(1.3);
    hIncrementGTSTime->SetStats(0);
    hIncrementGTSTime->GetYaxis()->ChangeLabel(1,-1.,-1.,-1,-1,-1, "< -4");
    hIncrementGTSTime->GetYaxis()->ChangeLabel(-1,-1.,-1.,-1,-1,-1, "> 4");
    hIncrementGTSTime->Draw("Text COLZ");
    TPaveText *text_hIncrementGTSTime = new TPaveText(0.3,0.76,0.8,0.89, "tblrNDC");
    text_hIncrementGTSTime->AddText("Number of GTS trailer per FEB:");
    text_hIncrementGTSTime->AddText(Form("0:%d; 1: %d, 2:%d, 3:%d, 4:%d, 5:%d, 6:%d,", N_GTS_trailer1_checked[0],N_GTS_trailer1_checked[1], N_GTS_trailer1_checked[2],
     N_GTS_trailer1_checked[3], N_GTS_trailer1_checked[4], N_GTS_trailer1_checked[5], N_GTS_trailer1_checked[6]));
    text_hIncrementGTSTime->AddText(Form("7:%d, 8:%d, 9:%d, 10:%d, 11:%d, 12:%d, 13:%d", N_GTS_trailer1_checked[7], N_GTS_trailer1_checked[8], N_GTS_trailer1_checked[9],
     N_GTS_trailer1_checked[10], N_GTS_trailer1_checked[11], N_GTS_trailer1_checked[12], N_GTS_trailer1_checked[13]));
    text_hIncrementGTSTime->AddText(Form("Total Nbre of Entries per FEB: %d", Nevent_checked));
    text_hIncrementGTSTime->SetBorderSize(1);
    text_hIncrementGTSTime->Draw("same");

    c1->Print(output_pdf);
    c1->Write();
    c1->Clear();

    gPad->SetRightMargin(0.1);
    c1->Divide(5, 3);
    for (int j_feb = 0; j_feb < N_FEB; j_feb++){
        c1->cd(j_feb+1);
        hHGAmplitude[j_feb]->Draw();
    }
    c1->Print(output_pdf);
    c1->Write();
    c1->Clear();

    c1->Divide(5, 3);
    for (int j_feb = 0; j_feb < N_FEB; j_feb++){
        c1->cd(j_feb+1);
        hLGAmplitude[j_feb]->Draw();
    }

    c1->Print(output_pdf);
    c1->Write();
    c1->Clear();

    c1->Print(Form("%s]", output_pdf));   // Finished; close the multipage file
    c1->Close();

    // close files
    fclose(read_fp);
    fclose(write_fp_words);
    fclose(write_fp_errors);

    file->cd();
    file->Write("");
    file->Close();

    // write to new file if we find gtime errors
    printf("%s decoded!\n", input_file);

    return 0;
}
