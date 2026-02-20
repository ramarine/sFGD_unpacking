
/*
 * data_decode_v3.cpp
 *
 * Standalone single-FEB unpacker for .bin/.daq files.
 *
 * Word assumptions (matching your earlier decoder style):
 *   - word ID: bits [31:28]
 *   - ID 0: Gate Header A/B, bit19 selects B (1) vs A (0)
 *   - ID 7: Gate time word, appears after header A (open) and after gate trailer (close)
 *   - ID 1: GTS header
 *   - ID 4: GTS trailer 1
 *   - ID 5: GTS trailer 2 (time/data)
 *   - ID 2: Hit time
 *   - ID 3: Hit amplitude
 *   - ID 6: Gate trailer
 *   - ID 11: Hold time, bit19 selects stop (1) vs start (0)
 *   - ID 13: FEB trailer (ends gate)
 *
 * Event (gate) boundary for single-FEB 
 *   Start: Gate Header A (ID0, bit19=0)
 *   End  : FEB Trailer (ID13)
 *
 * Outputs (prefix = derived from args):
 *   <prefix>LogWords.txt   (only detailed if VERBOSE=true)
 *   <prefix>LogErrors.txt
 *   <prefix>histo.root
 *   <prefix>histo.pdf
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <array>

#include <sys/stat.h>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1I.h"
#include "TH2D.h"

static bool VERBOSE  = false;   // set true for per-word logging
static bool READonly = false;   // if true, skip checks/histos (kept for parity)

static inline int getbits(uint32_t w, int startbit, int nbits) {
    uint32_t mask = (nbits == 32) ? 0xFFFFFFFFu : ((1u << nbits) - 1u);
    return (w >> startbit) & mask;
}

struct HitTime {
    int channelid=0; // [27:20]
    int hitid=0;     // [19:17]
    int tagid=0;     // [16:15]
    bool edge=false; // [14]
    int hittime=0;   // [12:0]
};

struct HitAmpl {
    int channelid=0; // [27:20]
    int hitid=0;     // [19:17]
    int tagid=0;     // [16:15]
    int gain=0;      // [14:12] (2=HG, 3=LG in your older code convention)
    int adc=0;       // [11:0]
};

struct GateState {
    bool started=false;
    int gatenumber=-1;

    // gate header A/B
    bool hasHeaderA=false, hasHeaderB=false;
    int Aboardid=0, Agatetype=0, Agatenumber=0;
    int Bboardid=0, Bgatetime=0;

    // gate time open/close (ID7)
    bool hasGateTimeOpen=false, hasGateTimeClose=false;
    int gateTimeOpen=0, gateTimeClose=0;

    // gate trailer
    bool hasGateTrailer=false;
    int Tboardid=0, Tgatetype=0, Tgatenumber=0;

    // FEB trailer
    bool hasFEBTrailer=false;
    bool Eventdonetimeout=false, D1FIFOfull=false, D0FIFOfull=false;
    int  Ndecodererror=0;

    // hits
    std::vector<HitTime> hit_times;
    std::vector<HitAmpl> hit_ampl;

    void reset() { *this = GateState{}; }
};

static bool is_directory(const std::string &path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return false;
    return S_ISDIR(st.st_mode);
}

static std::string basename_only(const std::string &p) {
    auto pos = p.find_last_of('/');
    if (pos == std::string::npos) return p;
    return p.substr(pos+1);
}

static std::string strip_extension(const std::string &p) {
    auto dot = p.find_last_of('.');
    if (dot == std::string::npos) return p;
    return p.substr(0, dot);
}

static std::string make_prefix(int argc, char **argv) {
    std::string in = argv[1];
    std::string base = strip_extension(basename_only(in));

    if (argc < 3) {
        // default: alongside input, as prefix = input path without extension
        return strip_extension(in);
    }

    std::string out = argv[2];

    // If argv[2] is an existing directory, use it + basename(input)
    if (is_directory(out)) {
        if (!out.empty() && out.back() != '/') out += "/";
        return out + base;
    }

    // Else treat argv[2] as prefix (can include path)
    return strip_extension(out);
}

static void check_gate_and_log(FILE *log, int &lineCounter, GateState &g, int &Nerrors,
                               int &NmissingA, int &NmissingB, int &NmissingOpen, int &NmissingTrl,
                               int &NmissingClose, int &NmissingFEB, int &NunmatchT, int &NunmatchA) {
    if (!g.started) return;

    auto err = [&](const char *fmt, ...) {
        va_list ap;
        va_start(ap, fmt);
        std::fprintf(log, "%d ERROR: ", ++lineCounter);
        std::vfprintf(log, fmt, ap);
        std::fprintf(log, "\n");
        va_end(ap);
        Nerrors++;
    };

    if (!g.hasHeaderA) { NmissingA++; err("Missing Gate Header A (ID0 A) for gate %d", g.gatenumber); }
    if (!g.hasHeaderB) { NmissingB++; err("Missing Gate Header B (ID0 B) for gate %d", g.gatenumber); }
    if (!g.hasGateTimeOpen) { NmissingOpen++; err("Missing Gate Time open (ID7) for gate %d", g.gatenumber); }
    if (!g.hasGateTrailer) { NmissingTrl++; err("Missing Gate Trailer (ID6) for gate %d", g.gatenumber); }
    if (!g.hasGateTimeClose) { NmissingClose++; err("Missing Gate Time close (ID7) for gate %d", g.gatenumber); }
    if (!g.hasFEBTrailer) { NmissingFEB++; err("Missing FEB Trailer (ID13) for gate %d", g.gatenumber); }

    // Pair hit time <-> amplitude by (channelid, hitid, tagid)
    for (const auto &t : g.hit_times) {
        bool found=false;
        for (const auto &a : g.hit_ampl) {
            if (a.channelid==t.channelid && a.hitid==t.hitid && a.tagid==t.tagid) { found=true; break; }
        }
        if (!found) { NunmatchT++; err("HitTime ch=%d hit=%d tag=%d has no matching amplitude [gate %d]",
                                       t.channelid, t.hitid, t.tagid, g.gatenumber); }
    }
    for (const auto &a : g.hit_ampl) {
        bool found=false;
        for (const auto &t : g.hit_times) {
            if (a.channelid==t.channelid && a.hitid==t.hitid && a.tagid==t.tagid) { found=true; break; }
        }
        if (!found) { NunmatchA++; err("HitAmpl ch=%d hit=%d tag=%d has no matching time [gate %d]",
                                       a.channelid, a.hitid, a.tagid, g.gatenumber); }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::printf("Usage: %s input.bin|input.daq [output_dir|output_prefix]\n", argv[0]);
        return 1;
    }

    const char *inputfile = argv[1];
    if (!(std::strstr(inputfile, ".bin") || std::strstr(inputfile, ".daq"))) {
        std::printf("Input file must be .bin or .daq\n");
        return 1;
    }

    std::string prefix = make_prefix(argc, argv);
    std::string logWords  = prefix + "LogWords.txt";
    std::string logErrors = prefix + "LogErrors.txt";
    std::string outRoot   = prefix + "histo.root";
    std::string outPdf    = prefix + "histo.pdf";

    FILE *in = std::fopen(inputfile, "rb");
    if (!in) { std::perror("Cannot open input"); return 1; }

    FILE *fw = std::fopen(logWords.c_str(), "w");
    FILE *fe = std::fopen(logErrors.c_str(), "w");
    if (!fw || !fe) { std::perror("Cannot open output logs"); return 1; }

    // ROOT
    TFile *rootfile = new TFile(outRoot.c_str(), "RECREATE");
    TCanvas *c1 = new TCanvas("c1", "Single FEB unpacker", 900, 600);

    TH1I *hGateOpen  = new TH1I("hGateOpen",  "Gate time open (ID7 after header A)",  2000, 0, 10000000);
    TH1I *hGateClose = new TH1I("hGateClose", "Gate time close (ID7 after trailer)",  2000, 0, 10000000);
    TH1I *hHitTime   = new TH1I("hHitTime",   "Hit time (leading edges)",             8192, 0, 8192);
    TH2D *hHitVsCh   = new TH2D("hHitVsCh",   "Hit time vs channel",                  256, -0.5, 255.5, 8192, 0, 8192);
    TH1I *hHG        = new TH1I("hHG",        "HG amplitude (ADC)",                   4096, 0, 4096);
    TH1I *hLG        = new TH1I("hLG",        "LG amplitude (ADC)",                   4096, 0, 4096);

    // Counters
    int nline_words = 0;
    int Nevent=0, Nerrors=0, Nempty=0;
    int NmissingA=0, NmissingB=0, NmissingOpen=0, NmissingTrl=0, NmissingClose=0, NmissingFEB=0;
    int NunmatchT=0, NunmatchA=0;

    GateState gate;
    bool prevWasGateTrailer = false; // used to interpret ID7 as open vs close

    unsigned char buf[4];

    while (std::fread(buf, 1, 4, in) == 4) {
        uint32_t w = ((uint32_t)buf[3] << 24) | ((uint32_t)buf[2] << 16) | ((uint32_t)buf[1] << 8) | (uint32_t)buf[0];
        if (w == 0x00000000u) { Nempty++; continue; }

        int id = getbits(w, 28, 4);

        if (VERBOSE) {
            std::fprintf(fw, "%d [ID %2d] 0x%08X\n", ++nline_words, id, w);
        }

        switch (id) {
            case 0: {
                bool isB = getbits(w, 19, 1);
                if (!isB) {
                    // Starting a new gate -> if previous gate existed, check it
                    if (gate.started && !READonly) {
                        check_gate_and_log(fw, nline_words, gate, Nerrors,
                                           NmissingA, NmissingB, NmissingOpen, NmissingTrl,
                                           NmissingClose, NmissingFEB, NunmatchT, NunmatchA);
                    }

                    gate.reset();
                    gate.started = true;

                    gate.Aboardid = getbits(w, 20, 8);
                    gate.Agatetype = getbits(w, 16, 3);
                    gate.Agatenumber = getbits(w, 0, 16);
                    gate.gatenumber = gate.Agatenumber;
                    gate.hasHeaderA = true;

                    prevWasGateTrailer = false;
                    Nevent++;

                    if (VERBOSE) std::fprintf(fw, "  GateHeaderA board=%d type=%d num=%d\n",
                                              gate.Aboardid, gate.Agatetype, gate.Agatenumber);
                } else {
                    if (!gate.started) break;
                    gate.Bboardid = getbits(w, 20, 8);
                    gate.Bgatetime = getbits(w, 0, 11);
                    gate.hasHeaderB = true;
                    if (VERBOSE) std::fprintf(fw, "  GateHeaderB board=%d time=%d\n", gate.Bboardid, gate.Bgatetime);
                }
                break;
            }

            case 7: {
                if (!gate.started) break;
                int t = getbits(w, 0, 28);
                if (!prevWasGateTrailer) {
                    gate.gateTimeOpen = t;
                    gate.hasGateTimeOpen = true;
                    hGateOpen->Fill(t);
                    if (VERBOSE) std::fprintf(fw, "  GateTimeOpen=%d\n", t);
                } else {
                    gate.gateTimeClose = t;
                    gate.hasGateTimeClose = true;
                    hGateClose->Fill(t);
                    if (VERBOSE) std::fprintf(fw, "  GateTimeClose=%d\n", t);
                }
                break;
            }

            case 6: {
                if (!gate.started) break;
                gate.Tboardid = getbits(w, 20, 8);
                gate.Tgatetype = getbits(w, 16, 3);
                gate.Tgatenumber = getbits(w, 0, 16);
                gate.hasGateTrailer = true;
                prevWasGateTrailer = true;
                if (VERBOSE) std::fprintf(fw, "  GateTrailer board=%d type=%d num=%d\n",
                                          gate.Tboardid, gate.Tgatetype, gate.Tgatenumber);
                break;
            }

            case 11: {
                // Hold time (not used further, but decoded if needed)
                // bit19: 0 start, 1 stop; board [27:20]; time [18:0]
                // You can extend this section if you want to log/store it.
                break;
            }

            case 2: {
                if (!gate.started) break;
                HitTime ht;
                ht.channelid = getbits(w, 20, 8);
                ht.hitid     = getbits(w, 17, 3);
                ht.tagid     = getbits(w, 15, 2);
                ht.edge      = getbits(w, 14, 1);
                ht.hittime   = getbits(w, 0, 13);
                gate.hit_times.push_back(ht);
                if (!ht.edge) {
                    hHitTime->Fill(ht.hittime);
                    hHitVsCh->Fill(ht.channelid, ht.hittime);
                }
                break;
            }

            case 3: {
                if (!gate.started) break;
                HitAmpl ha;
                ha.channelid = getbits(w, 20, 8);
                ha.hitid     = getbits(w, 17, 3);
                ha.tagid     = getbits(w, 15, 2);
                ha.gain      = getbits(w, 12, 3);
                ha.adc       = getbits(w, 0, 12);
                gate.hit_ampl.push_back(ha);
                if (ha.gain == 2) hHG->Fill(ha.adc);
                if (ha.gain == 3) hLG->Fill(ha.adc);
                break;
            }

            case 13: {
                if (!gate.started) break;

                // FEB trailer bit mapping (following your older datadecodev2 description):
                // Eventdonetimeout [18], D1FIFOfull [17], D0FIFOfull [16], Ndecodererror [15:0]
                gate.Eventdonetimeout = getbits(w, 18, 1);
                gate.D1FIFOfull       = getbits(w, 17, 1);
                gate.D0FIFOfull       = getbits(w, 16, 1);
                gate.Ndecodererror    = getbits(w, 0, 16);
                gate.hasFEBTrailer    = true;

                // End of gate -> checks
                if (!READonly) {
                    check_gate_and_log(fw, nline_words, gate, Nerrors,
                                       NmissingA, NmissingB, NmissingOpen, NmissingTrl,
                                       NmissingClose, NmissingFEB, NunmatchT, NunmatchA);
                }

                // Prepare for next gate
                prevWasGateTrailer = false;
                break;
            }

            default:
                // ignore other IDs for now (GTS etc. can be added later if needed)
                break;
        }
    }

    // Final: if file ended mid-gate, still check what we have
    if (gate.started && !READonly) {
        check_gate_and_log(fw, nline_words, gate, Nerrors,
                           NmissingA, NmissingB, NmissingOpen, NmissingTrl,
                           NmissingClose, NmissingFEB, NunmatchT, NunmatchA);
    }

    // Error summary
    std::fprintf(fe, "========================================\n");
    std::fprintf(fe, "Single-FEB unpacker summary\n");
    std::fprintf(fe, "========================================\n");
    std::fprintf(fe, "Gates read              : %d\n", Nevent);
    std::fprintf(fe, "Total errors            : %d\n", Nerrors);
    std::fprintf(fe, "Empty words             : %d\n", Nempty);
    std::fprintf(fe, "Missing Gate Header A   : %d\n", NmissingA);
    std::fprintf(fe, "Missing Gate Header B   : %d\n", NmissingB);
    std::fprintf(fe, "Missing Gate Time open  : %d\n", NmissingOpen);
    std::fprintf(fe, "Missing Gate Trailer    : %d\n", NmissingTrl);
    std::fprintf(fe, "Missing Gate Time close : %d\n", NmissingClose);
    std::fprintf(fe, "Missing FEB Trailer     : %d\n", NmissingFEB);
    std::fprintf(fe, "Unmatched HitTime       : %d\n", NunmatchT);
    std::fprintf(fe, "Unmatched HitAmplitude  : %d\n", NunmatchA);

    // Save histograms
    c1->Print((outPdf + "[").c_str());
    hGateOpen->Draw();  c1->Print(outPdf.c_str());
    hGateClose->Draw(); c1->Print(outPdf.c_str());
    hHitTime->Draw();   c1->Print(outPdf.c_str());
    hHitVsCh->Draw("COLZ"); c1->Print(outPdf.c_str());
    hHG->Draw();        c1->Print(outPdf.c_str());
    hLG->Draw();        c1->Print(outPdf.c_str());
    c1->Print((outPdf + "]").c_str());

    rootfile->Write();
    rootfile->Close();

    std::fclose(in);
    std::fclose(fw);
    std::fclose(fe);

    std::printf("Done.\n");
    std::printf("  LogWords : %s\n", logWords.c_str());
    std::printf("  LogErrors: %s\n", logErrors.c_str());
    std::printf("  ROOT     : %s\n", outRoot.c_str());
    std::printf("  PDF      : %s\n", outPdf.c_str());

    return 0;
}

