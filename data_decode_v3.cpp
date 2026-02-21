#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include <limits>
#include <sys/stat.h>

#include "TFile.h"
#include "TTree.h"

static inline uint32_t u32_from_le4(const unsigned char b[4]) {
  return (uint32_t)b[0] | ((uint32_t)b[1] << 8) | ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

static inline int getbits(uint32_t w, int startbit, int nbits) {
  uint32_t mask = (nbits == 32) ? 0xFFFFFFFFu : ((1u << nbits) - 1u);
  return (w >> startbit) & mask;
}

static bool is_directory(const std::string &path) {
  struct stat st;
  if (stat(path.c_str(), &st) != 0) return false;
  return S_ISDIR(st.st_mode);
}

static std::string basename_only(const std::string &p) {
  auto pos = p.find_last_of('/');
  if (pos == std::string::npos) return p;
  return p.substr(pos + 1);
}

static std::string strip_extension(const std::string &p) {
  auto dot = p.find_last_of('.');
  if (dot == std::string::npos) return p;
  return p.substr(0, dot);
}

static std::string make_prefix(int argc, char **argv) {
  std::string in = argv[1];
  std::string base = strip_extension(basename_only(in));

  if (argc < 3) return strip_extension(in);

  std::string out = argv[2];
  if (is_directory(out)) {
    if (!out.empty() && out.back() != '/') out += "/";
    return out + base;
  }
  return strip_extension(out);
}

// -------------------- Hit aggregation --------------------
struct HitKey {
  int ch;
  int hit;
  int tag;
  bool operator==(const HitKey &o) const { return ch==o.ch && hit==o.hit && tag==o.tag; }
};

struct HitKeyHash {
  std::size_t operator()(const HitKey &k) const noexcept {
    // small ints; pack into 32b then hash
    uint32_t x = (uint32_t)(k.ch & 0xFF) | ((uint32_t)(k.hit & 0x7) << 8) | ((uint32_t)(k.tag & 0x3) << 11);
    // mix
    x ^= x >> 16; x *= 0x7feb352dU; x ^= x >> 15; x *= 0x846ca68bU; x ^= x >> 16;
    return (std::size_t)x;
  }
};

struct HitAccum {
  // times
  bool has_le=false, has_te=false;
  int  t_le_tick=0, t_te_tick=0;
  int  n_le=0, n_te=0;

  // amplitudes
  bool has_hg=false, has_lg=false;
  int  adc_hg=0, adc_lg=0;
  int  n_hg=0, n_lg=0;
};

// -------------------- Gate state --------------------
struct GateState {
  bool started=false;
  int gate_number=-1;

  // Gate header A / B (word ID 0)
  bool hasHeaderA=false, hasHeaderB=false;
  int Aboardid=0, Agatetype=0, Agatenumber=0;
  int Bboardid=0, Bgatetime=0; // 11-bit time in header B

  // Gate time (word ID 7)
  bool hasGateTimeOpen=false, hasGateTimeClose=false;
  uint32_t gateTimeOpen=0, gateTimeClose=0;

  // Gate trailer (word ID 6) + trailer time (word ID 7)
  bool hasGateTrailer=false;
  int Tboardid=0, Tgatetype=0, Tgatenumber=0;

  // FEB trailer (word ID 13)
  bool hasFEBTrailer=false;
  int Eventdonetimeout=0, D1FIFOfull=0, D0FIFOfull=0;
  int Ndecodererror=0;

  // GTS anchors (we store first-seen inside gate + last-seen overall)
  bool hasGTS_first_inGate=false;
  uint32_t gts_tag_first_inGate=0;
  uint32_t gts_time_first_inGate=0;
  int gts_data_first_inGate=0;

  bool hasGTS_last=false;
  uint32_t gts_tag_last=0;
  uint32_t gts_time_last=0;
  int gts_data_last=0;

  // hit aggregation
  std::unordered_map<HitKey, HitAccum, HitKeyHash> hits;
  uint32_t n_hit_time_words=0;
  uint32_t n_hit_amp_words=0;

  // gate-trailer boundary flag (for ID7 open/close and for GTS inGate classification)
  bool afterGateTrailer=false;

  void reset() { *this = GateState{}; }
};

// -------------------- main --------------------
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
  std::string outRoot = prefix + "trees.root";

  FILE *in = std::fopen(inputfile, "rb");
  if (!in) { std::perror("Cannot open input"); return 1; }

  // Constants (your confirmed hit time unit)
  constexpr double HIT_TICK_NS = 2.5;

  // These are kept as *convenience* derived columns; raw ticks are always stored.
  // Reference code labels: gate time ~10 ms resolution, GTS time diff ~10 us resolution.
  constexpr double GATE_TICK_MS_ASSUMED = 10.0;
  constexpr double GTS_TICK_US_ASSUMED  = 10.0;

  // ROOT file + TTrees
  TFile *f = new TFile(outRoot.c_str(), "RECREATE");

  // ---------------- Gate tree ----------------
  TTree *tGate = new TTree("Gate", "One entry per gate (single FEB)");
  Int_t   gate_number=0, gate_type=0;
  Int_t   Aboardid=0, Bboardid=0, Tboardid=0;
  Int_t   Bgatetime=0;
  UInt_t  gate_open_tick28=0, gate_close_tick28=0;
  Double_t gate_open_ms_assumed=0, gate_close_ms_assumed=0;
  Int_t   Tgatetype=0, Tgatenumber=0;
  Int_t   hasHeaderA=0, hasHeaderB=0, hasGateOpen=0, hasGateTrailer=0, hasGateClose=0, hasFEBTrailer=0;

  // FEB trailer flags
  Int_t Eventdonetimeout=0, D1FIFOfull=0, D0FIFOfull=0, Ndecodererror=0;

  // GTS anchors
  Int_t hasGTS_first_inGate=0;
  UInt_t gts_tag_first_inGate=0, gts_time_first_inGate=0;
  Int_t gts_data_first_inGate=0;

  Int_t hasGTS_last=0;
  UInt_t gts_tag_last=0, gts_time_last=0;
  Int_t gts_data_last=0;

  // counts
  UInt_t n_hit_keys=0, n_hit_time_words=0, n_hit_amp_words=0;

  tGate->Branch("gate_number", &gate_number, "gate_number/I");
  tGate->Branch("gate_type",   &gate_type,   "gate_type/I");
  tGate->Branch("Aboardid",    &Aboardid,    "Aboardid/I");
  tGate->Branch("Bboardid",    &Bboardid,    "Bboardid/I");
  tGate->Branch("Bgatetime",   &Bgatetime,   "Bgatetime/I");
  tGate->Branch("Tboardid",    &Tboardid,    "Tboardid/I");
  tGate->Branch("Tgatetype",   &Tgatetype,   "Tgatetype/I");
  tGate->Branch("Tgatenumber", &Tgatenumber, "Tgatenumber/I");

  tGate->Branch("gate_open_tick28",  &gate_open_tick28,  "gate_open_tick28/i");
  tGate->Branch("gate_close_tick28", &gate_close_tick28, "gate_close_tick28/i");
  tGate->Branch("gate_open_ms_assumed",  &gate_open_ms_assumed,  "gate_open_ms_assumed/D");
  tGate->Branch("gate_close_ms_assumed", &gate_close_ms_assumed, "gate_close_ms_assumed/D");

  tGate->Branch("hasHeaderA",    &hasHeaderA,    "hasHeaderA/I");
  tGate->Branch("hasHeaderB",    &hasHeaderB,    "hasHeaderB/I");
  tGate->Branch("hasGateOpen",   &hasGateOpen,   "hasGateOpen/I");
  tGate->Branch("hasGateTrailer",&hasGateTrailer,"hasGateTrailer/I");
  tGate->Branch("hasGateClose",  &hasGateClose,  "hasGateClose/I");
  tGate->Branch("hasFEBTrailer", &hasFEBTrailer, "hasFEBTrailer/I");

  tGate->Branch("Eventdonetimeout", &Eventdonetimeout, "Eventdonetimeout/I");
  tGate->Branch("D1FIFOfull",       &D1FIFOfull,       "D1FIFOfull/I");
  tGate->Branch("D0FIFOfull",       &D0FIFOfull,       "D0FIFOfull/I");
  tGate->Branch("Ndecodererror",    &Ndecodererror,    "Ndecodererror/I");

  tGate->Branch("hasGTS_first_inGate", &hasGTS_first_inGate, "hasGTS_first_inGate/I");
  tGate->Branch("gts_tag_first_inGate", &gts_tag_first_inGate, "gts_tag_first_inGate/i");
  tGate->Branch("gts_time_first_inGate",&gts_time_first_inGate,"gts_time_first_inGate/i");
  tGate->Branch("gts_data_first_inGate",&gts_data_first_inGate,"gts_data_first_inGate/I");

  tGate->Branch("hasGTS_last", &hasGTS_last, "hasGTS_last/I");
  tGate->Branch("gts_tag_last", &gts_tag_last, "gts_tag_last/i");
  tGate->Branch("gts_time_last",&gts_time_last,"gts_time_last/i");
  tGate->Branch("gts_data_last",&gts_data_last,"gts_data_last/I");

  tGate->Branch("n_hit_keys",      &n_hit_keys,      "n_hit_keys/i");
  tGate->Branch("n_hit_time_words",&n_hit_time_words,"n_hit_time_words/i");
  tGate->Branch("n_hit_amp_words", &n_hit_amp_words, "n_hit_amp_words/i");

  // ---------------- Hit tree ----------------
  TTree *tHit = new TTree("Hit", "One entry per hit key (ch,hitid,tagid) per gate with LE/TE and HG/LG");

  Int_t   hit_gate_number=0;
  Int_t   channelid=0, hitid=0, tagid=0;
  Int_t   has_le=0, has_te=0;
  Int_t   t_le_tick=0, t_te_tick=0;
  Double_t t_le_ns=0, t_te_ns=0;
  Int_t   tot_tick=0;
  Double_t tot_ns=0;

  Int_t   has_hg=0, has_lg=0;
  Int_t   adc_hg=0, adc_lg=0;

  Int_t   n_le=0, n_te=0, n_hg=0, n_lg=0;

  // Helpful anchors for “absolute-ish” timelines without joining trees
  UInt_t  hit_gate_open_tick28=0;
  Double_t hit_gate_open_ms_assumed=0;

  Int_t   hit_hasGTS_first_inGate=0;
  UInt_t  hit_gts_tag_first_inGate=0, hit_gts_time_first_inGate=0;
  Double_t hit_gts_time_us_assumed=0;

  tHit->Branch("gate_number", &hit_gate_number, "gate_number/I");
  tHit->Branch("channelid", &channelid, "channelid/I");
  tHit->Branch("hitid",     &hitid,     "hitid/I");
  tHit->Branch("tagid",     &tagid,     "tagid/I");

  tHit->Branch("has_le", &has_le, "has_le/I");
  tHit->Branch("has_te", &has_te, "has_te/I");
  tHit->Branch("t_le_tick", &t_le_tick, "t_le_tick/I");
  tHit->Branch("t_te_tick", &t_te_tick, "t_te_tick/I");
  tHit->Branch("t_le_ns", &t_le_ns, "t_le_ns/D");
  tHit->Branch("t_te_ns", &t_te_ns, "t_te_ns/D");
  tHit->Branch("tot_tick", &tot_tick, "tot_tick/I");
  tHit->Branch("tot_ns", &tot_ns, "tot_ns/D");

  tHit->Branch("has_hg", &has_hg, "has_hg/I");
  tHit->Branch("has_lg", &has_lg, "has_lg/I");
  tHit->Branch("adc_hg", &adc_hg, "adc_hg/I");
  tHit->Branch("adc_lg", &adc_lg, "adc_lg/I");

  tHit->Branch("n_le", &n_le, "n_le/I");
  tHit->Branch("n_te", &n_te, "n_te/I");
  tHit->Branch("n_hg", &n_hg, "n_hg/I");
  tHit->Branch("n_lg", &n_lg, "n_lg/I");

  tHit->Branch("gate_open_tick28", &hit_gate_open_tick28, "gate_open_tick28/i");
  tHit->Branch("gate_open_ms_assumed", &hit_gate_open_ms_assumed, "gate_open_ms_assumed/D");

  tHit->Branch("hasGTS_first_inGate", &hit_hasGTS_first_inGate, "hasGTS_first_inGate/I");
  tHit->Branch("gts_tag_first_inGate", &hit_gts_tag_first_inGate, "gts_tag_first_inGate/i");
  tHit->Branch("gts_time_first_inGate",&hit_gts_time_first_inGate,"gts_time_first_inGate/i");
  tHit->Branch("gts_time_us_assumed",  &hit_gts_time_us_assumed,  "gts_time_us_assumed/D");

  // ---------------- GTS tree ----------------
  TTree *tGTS = new TTree("GTS", "GTS words stream (tag/time/header) with gate association");
  Int_t   gts_gate_number=0;
  Int_t   gts_wordid=0;       // 1,4,5
  Int_t   gts_inGate=0;       // 1 if before gate trailer
  UInt_t  gts_tag=0;          // for ID 1 and ID 4
  UInt_t  gts_time20=0;       // for ID 5
  Int_t   gts_data=0;         // for ID 5
  Double_t gts_time_us_assumed=0;

  tGTS->Branch("gate_number", &gts_gate_number, "gate_number/I");
  tGTS->Branch("wordid",      &gts_wordid,      "wordid/I");
  tGTS->Branch("inGate",      &gts_inGate,      "inGate/I");
  tGTS->Branch("GTStag",      &gts_tag,         "GTStag/i");
  tGTS->Branch("GTStime20",   &gts_time20,      "GTStime20/i");
  tGTS->Branch("Data",        &gts_data,        "Data/I");
  tGTS->Branch("GTStime_us_assumed", &gts_time_us_assumed, "GTStime_us_assumed/D");

  // ---------------- Decode loop ----------------
  GateState gate;

  auto flush_gate_to_trees = [&]() {
    if (!gate.started) return;

    // Fill Gate tree
    gate_number = gate.gate_number;
    gate_type   = gate.Agatetype;

    Aboardid = gate.Aboardid;
    Bboardid = gate.Bboardid;
    Bgatetime = gate.Bgatetime;

    gate_open_tick28  = gate.gateTimeOpen;
    gate_close_tick28 = gate.gateTimeClose;
    gate_open_ms_assumed  = (double)gate_open_tick28  * GATE_TICK_MS_ASSUMED;
    gate_close_ms_assumed = (double)gate_close_tick28 * GATE_TICK_MS_ASSUMED;

    Tboardid = gate.Tboardid;
    Tgatetype = gate.Tgatetype;
    Tgatenumber = gate.Tgatenumber;

    hasHeaderA = gate.hasHeaderA ? 1 : 0;
    hasHeaderB = gate.hasHeaderB ? 1 : 0;
    hasGateOpen = gate.hasGateTimeOpen ? 1 : 0;
    hasGateTrailer = gate.hasGateTrailer ? 1 : 0;
    hasGateClose = gate.hasGateTimeClose ? 1 : 0;
    hasFEBTrailer = gate.hasFEBTrailer ? 1 : 0;

    Eventdonetimeout = gate.Eventdonetimeout;
    D1FIFOfull       = gate.D1FIFOfull;
    D0FIFOfull       = gate.D0FIFOfull;
    Ndecodererror    = gate.Ndecodererror;

    hasGTS_first_inGate = gate.hasGTS_first_inGate ? 1 : 0;
    gts_tag_first_inGate = gate.gts_tag_first_inGate;
    gts_time_first_inGate = gate.gts_time_first_inGate;
    gts_data_first_inGate = gate.gts_data_first_inGate;

    hasGTS_last = gate.hasGTS_last ? 1 : 0;
    gts_tag_last = gate.gts_tag_last;
    gts_time_last = gate.gts_time_last;
    gts_data_last = gate.gts_data_last;

    n_hit_keys = (UInt_t)gate.hits.size();
    n_hit_time_words = gate.n_hit_time_words;
    n_hit_amp_words  = gate.n_hit_amp_words;

    tGate->Fill();

    // Fill Hit tree (one row per key)
    for (const auto &kv : gate.hits) {
      const HitKey &k = kv.first;
      const HitAccum &h = kv.second;

      hit_gate_number = gate.gate_number;
      channelid = k.ch;
      hitid = k.hit;
      tagid = k.tag;

      has_le = h.has_le ? 1 : 0;
      has_te = h.has_te ? 1 : 0;
      t_le_tick = h.t_le_tick;
      t_te_tick = h.t_te_tick;
      t_le_ns = h.has_le ? (double)h.t_le_tick * HIT_TICK_NS : std::numeric_limits<double>::quiet_NaN();
      t_te_ns = h.has_te ? (double)h.t_te_tick * HIT_TICK_NS : std::numeric_limits<double>::quiet_NaN();

      if (h.has_le && h.has_te) {
        tot_tick = h.t_te_tick - h.t_le_tick;
        tot_ns = (double)tot_tick * HIT_TICK_NS;
      } else {
        tot_tick = 0;
        tot_ns = std::numeric_limits<double>::quiet_NaN();
      }

      has_hg = h.has_hg ? 1 : 0;
      has_lg = h.has_lg ? 1 : 0;
      adc_hg = h.adc_hg;
      adc_lg = h.adc_lg;

      n_le = h.n_le;
      n_te = h.n_te;
      n_hg = h.n_hg;
      n_lg = h.n_lg;

      hit_gate_open_tick28 = gate.gateTimeOpen;
      hit_gate_open_ms_assumed = (double)hit_gate_open_tick28 * GATE_TICK_MS_ASSUMED;

      hit_hasGTS_first_inGate = gate.hasGTS_first_inGate ? 1 : 0;
      hit_gts_tag_first_inGate  = gate.gts_tag_first_inGate;
      hit_gts_time_first_inGate = gate.gts_time_first_inGate;
      hit_gts_time_us_assumed   = (double)hit_gts_time_first_inGate * GTS_TICK_US_ASSUMED;

      tHit->Fill();
    }
  };

  unsigned char buf[4];
  uint32_t prevword = 0;

  while (std::fread(buf, 1, 4, in) == 4) {
    uint32_t w = u32_from_le4(buf);
    if (w == 0x00000000u) { prevword = w; continue; }

    int id = getbits(w, 28, 4);

    switch (id) {
      case 0: {
        // Gate Header A/B; bit19 selects B (1) vs A (0)
        bool isB = getbits(w, 19, 1);

        if (!isB) {
          // New gate starts: flush previous gate if any
          if (gate.started) flush_gate_to_trees();

          gate.reset();
          gate.started = true;

          gate.Aboardid = getbits(w, 20, 8);
          gate.Agatetype = getbits(w, 16, 3);
          gate.Agatenumber = getbits(w, 0, 16);
          gate.gate_number = gate.Agatenumber;
          gate.hasHeaderA = true;

          gate.afterGateTrailer = false;
        } else {
          if (!gate.started) break;
          gate.Bboardid = getbits(w, 20, 8);
          gate.Bgatetime = getbits(w, 0, 11);
          gate.hasHeaderB = true;
        }
        break;
      }

      case 7: {
        // Gate time word: interpret as open vs close using afterGateTrailer flag
        if (!gate.started) break;
        uint32_t t = (uint32_t)getbits(w, 0, 28);

        if (!gate.afterGateTrailer) {
          gate.gateTimeOpen = t;
          gate.hasGateTimeOpen = true;
        } else {
          gate.gateTimeClose = t;
          gate.hasGateTimeClose = true;
        }
        break;
      }

      case 6: {
        // Gate trailer
        if (!gate.started) break;
        gate.Tboardid = getbits(w, 20, 8);
        gate.Tgatetype = getbits(w, 16, 3);
        gate.Tgatenumber = getbits(w, 0, 16);
        gate.hasGateTrailer = true;
        gate.afterGateTrailer = true;
        break;
      }

      case 2: {
        // Hit time: (channelid, hitid, tagid, edge, hittime)
        if (!gate.started) break;
        gate.n_hit_time_words++;

        int ch  = getbits(w, 20, 8);
        int hid = getbits(w, 17, 3);
        int tid = getbits(w, 15, 2);
        int edge = getbits(w, 14, 1);     // convention: 0=LE, 1=TE
        int ht  = getbits(w, 0, 13);

        HitKey key{ch, hid, tid};
        HitAccum &acc = gate.hits[key];

        if (edge == 0) {
          acc.n_le++;
          if (!acc.has_le) { acc.has_le = true; acc.t_le_tick = ht; }
          // if duplicates happen, keep the earliest LE
          else if (ht < acc.t_le_tick) { acc.t_le_tick = ht; }
        } else {
          acc.n_te++;
          if (!acc.has_te) { acc.has_te = true; acc.t_te_tick = ht; }
          // if duplicates happen, keep the latest TE
          else if (ht > acc.t_te_tick) { acc.t_te_tick = ht; }
        }
        break;
      }

      case 3: {
        // Hit amplitude: (channelid, hitid, tagid, amplitudeid, amplitudemeas)
        if (!gate.started) break;
        gate.n_hit_amp_words++;

        int ch  = getbits(w, 20, 8);
        int hid = getbits(w, 17, 3);
        int tid = getbits(w, 15, 2);
        int ampId = getbits(w, 12, 3);
        int adc   = getbits(w, 0, 12);

        HitKey key{ch, hid, tid};
        HitAccum &acc = gate.hits[key];

        // Convention used in your reference: 2=HG, 3=LG
        if (ampId == 2) {
          acc.n_hg++;
          acc.has_hg = true;
          acc.adc_hg = adc;
        } else if (ampId == 3) {
          acc.n_lg++;
          acc.has_lg = true;
          acc.adc_lg = adc;
        }
        break;
      }

      case 1: {
        // GTS header: tag in [27:0]
        if (!gate.started) break;
        gts_gate_number = gate.gate_number;
        gts_wordid = 1;
        gts_inGate = gate.afterGateTrailer ? 0 : 1;
        gts_tag = (uint32_t)getbits(w, 0, 28);
        gts_time20 = 0;
        gts_data = 0;
        gts_time_us_assumed = 0.0;
        tGTS->Fill();
        break;
      }

      case 4: {
        // GTS trailer 1: tag in [27:0]
        if (!gate.started) break;

        uint32_t tag = (uint32_t)getbits(w, 0, 28);

        gts_gate_number = gate.gate_number;
        gts_wordid = 4;
        gts_inGate = gate.afterGateTrailer ? 0 : 1;
        gts_tag = tag;
        gts_time20 = 0;
        gts_data = 0;
        gts_time_us_assumed = 0.0;
        tGTS->Fill();

        // keep tag in state as "last seen"; time comes with ID 5
        gate.hasGTS_last = true;
        gate.gts_tag_last = tag;
        break;
      }

      case 5: {
        // GTS trailer 2: Data bit [27], time [19:0]
        if (!gate.started) break;

        int data = getbits(w, 27, 1);
        uint32_t t20 = (uint32_t)getbits(w, 0, 20);

        gts_gate_number = gate.gate_number;
        gts_wordid = 5;
        gts_inGate = gate.afterGateTrailer ? 0 : 1;
        // associate tag from most recent ID4 if available (otherwise leave as 0)
        gts_tag = gate.hasGTS_last ? gate.gts_tag_last : 0;
        gts_time20 = t20;
        gts_data = data;
        gts_time_us_assumed = (double)t20 * GTS_TICK_US_ASSUMED;
        tGTS->Fill();

        // update last + possibly first-inGate anchor
        gate.hasGTS_last = true;
        gate.gts_time_last = t20;
        gate.gts_data_last = data;

        if (!gate.afterGateTrailer && !gate.hasGTS_first_inGate) {
          gate.hasGTS_first_inGate = true;
          gate.gts_tag_first_inGate = gate.hasGTS_last ? gate.gts_tag_last : 0;
          gate.gts_time_first_inGate = t20;
          gate.gts_data_first_inGate = data;
        }
        break;
      }

      case 13: {
        // FEB trailer ends the gate (single-FEB assumption)
        if (!gate.started) break;
        gate.Eventdonetimeout = getbits(w, 18, 1);
        gate.D1FIFOfull       = getbits(w, 17, 1);
        gate.D0FIFOfull       = getbits(w, 16, 1);
        gate.Ndecodererror    = getbits(w, 0, 16);
        gate.hasFEBTrailer    = true;

        // Flush and reset for next gate
        flush_gate_to_trees();
        gate.reset();
        break;
      }

      default:
        // ignore other IDs for this single-FEB decoder
        break;
    }

    prevword = w;
  }

  // If file ends mid-gate, still flush what we have
  if (gate.started) {
    flush_gate_to_trees();
  }

  std::fclose(in);

  f->Write();
  f->Close();

  std::printf("Done.\n");
  std::printf("  ROOT file: %s\n", outRoot.c_str());
  std::printf("  Trees    : Gate, Hit, GTS\n");
  return 0;
}
