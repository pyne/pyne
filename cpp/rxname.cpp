#include "rxname.h"

std::string pyne::rxname::_names[NUM_RX_NAMES] = {
  "total",
  "elastic",
  "nonelastic",
  "inelastic",
  "misc",
  "continuum",
  "z_2nd",
  "z_2n",
  "z_3n",
  "fission",
  "fission_first",
  "fission_second",
  "fission_third",
  "na",
  "n3a",
  "z_2na",
  "z_3na",
  "absorption",
  "np",
  "n2a",
  "z_2n2a",
  "nd",
  "nt",
  "nHe3",
  "nd3a",
  "nt2a",
  "z_4n",
  "fission_fourth",
  "z_2np",
  "z_3np",
  "n2p",
  "npa",
  "n_0",
  "n_1",
  "n_2",
  "n_3",
  "n_4",
  "n_5",
  "n_6",
  "n_7",
  "n_8",
  "n_9",
  "n_10",
  "n_11",
  "n_12",
  "n_13",
  "n_14",
  "n_15",
  "n_16",
  "n_17",
  "n_18",
  "n_19",
  "n_20",
  "n_21",
  "n_22",
  "n_23",
  "n_24",
  "n_25",
  "n_26",
  "n_27",
  "n_28",
  "n_29",
  "n_30",
  "n_31",
  "n_32",
  "n_33",
  "n_34",
  "n_35",
  "n_36",
  "n_37",
  "n_38",
  "n_39",
  "n_40",
  "n_continuum",
  "disappeareance",
  "gamma",
  "p",
  "d",
  "t",
  "He3",
  "a",
  "z_2a",
  "z_3a",
  "z_2p",
  "pa",
  "t2a",
  "d2a",
  "pd",
  "pt",
  "da",
  "resonance_parameters",
  "n_total",
  "gamma_total",
  "p_total",
  "d_total",
  "t_total",
  "He3_total",
  "a_total",
  "pionp",
  "pion0",
  "pionm",
  "muonp",
  "muonm",
  "kaonp",
  "kaon0_long",
  "kaon0_short",
  "kaonm",
  "antip",
  "antin",
  "mubar",
  "epsilon",
  "y",
  "damage",
  "heading",
  "nubar",
  "fission_product_yield_independent",
  "nubar_delayed",
  "nubar_prompt",
  "decay",
  "energy_per_fission",
  "fission_product_yield_cumulative",
  "gamma_delayed",
  "stopping_power",
  "photon_total",
  "photon_coherent",
  "photon_incoherent",
  "scattering_factor_imag",
  "scattering_factor_real",
  "pair_prod_elec",
  "pair_prod",
  "pair_prod_nuc",
  "absorption_photoelectric",
  "photoexcitation",
  "scattering_electroatomic",
  "bremsstrahlung",
  "excitation_electroatomic",
  "atomic_relaxation",
  "k_photoelectric",
  "l1_photoelectric",
  "l2_photoelectric",
  "l3_photoelectric",
  "m1_photoelectric",
  "m2_photoelectric",
  "m3_photoelectric",
  "m4_photoelectric",
  "m5_photoelectric",
  "n1_photoelectric",
  "n2_photoelectric",
  "n3_photoelectric",
  "n4_photoelectric",
  "n5_photoelectric",
  "n6_photoelectric",
  "n7_photoelectric",
  "o1_photoelectric",
  "o2_photoelectric",
  "o3_photoelectric",
  "o4_photoelectric",
  "o5_photoelectric",
  "o6_photoelectric",
  "o7_photoelectric",
  "o8_photoelectric",
  "o9_photoelectric",
  "p1_photoelectric",
  "p2_photoelectric",
  "p3_photoelectric",
  "p4_photoelectric",
  "p5_photoelectric",
  "p6_photoelectric",
  "p7_photoelectric",
  "p8_photoelectric",
  "p9_photoelectric",
  "p10_photoelectric",
  "p11_photoelectric",
  "q1_photoelectric",
  "q2_photoelectric",
  "q3_photoelectric",
  "p_0",
  "p_1",
  "p_2",
  "p_3",
  "p_4",
  "p_5",
  "p_6",
  "p_7",
  "p_8",
  "p_9",
  "p_10",
  "p_11",
  "p_12",
  "p_13",
  "p_14",
  "p_15",
  "p_16",
  "p_17",
  "p_18",
  "p_19",
  "p_20",
  "p_21",
  "p_22",
  "p_23",
  "p_24",
  "p_25",
  "p_26",
  "p_27",
  "p_28",
  "p_29",
  "p_30",
  "p_31",
  "p_32",
  "p_33",
  "p_34",
  "p_35",
  "p_36",
  "p_37",
  "p_38",
  "p_39",
  "p_40",
  "p_41",
  "p_42",
  "p_43",
  "p_44",
  "p_45",
  "p_46",
  "p_47",
  "p_48",
  "p_continuum",
  "d_0",
  "d_1",
  "d_2",
  "d_3",
  "d_4",
  "d_5",
  "d_6",
  "d_7",
  "d_8",
  "d_9",
  "d_10",
  "d_11",
  "d_12",
  "d_13",
  "d_14",
  "d_15",
  "d_16",
  "d_17",
  "d_18",
  "d_19",
  "d_20",
  "d_21",
  "d_22",
  "d_23",
  "d_24",
  "d_25",
  "d_26",
  "d_27",
  "d_28",
  "d_29",
  "d_30",
  "d_31",
  "d_32",
  "d_33",
  "d_34",
  "d_35",
  "d_36",
  "d_37",
  "d_38",
  "d_39",
  "d_40",
  "d_41",
  "d_42",
  "d_43",
  "d_44",
  "d_45",
  "d_46",
  "d_47",
  "d_48",
  "d_continuum",
  "t_continuum",
  "He3_0",
  "He3_1",
  "He3_2",
  "He3_3",
  "He3_4",
  "He3_5",
  "He3_6",
  "He3_7",
  "He3_8",
  "He3_9",
  "He3_10",
  "He3_11",
  "He3_12",
  "He3_13",
  "He3_14",
  "He3_15",
  "He3_16",
  "He3_17",
  "He3_18",
  "He3_19",
  "He3_20",
  "He3_21",
  "He3_22",
  "He3_23",
  "He3_24",
  "He3_25",
  "He3_26",
  "He3_27",
  "He3_28",
  "He3_29",
  "He3_30",
  "He3_31",
  "He3_32",
  "He3_33",
  "He3_34",
  "He3_35",
  "He3_36",
  "He3_37",
  "He3_38",
  "He3_39",
  "He3_40",
  "He3_41",
  "He3_42",
  "He3_43",
  "He3_44",
  "He3_45",
  "He3_46",
  "He3_47",
  "He3_48",
  "He3_continuum",
  "a_continuum"
  };
std::set<std::string> pyne::rxname::names(pyne::rxname::_names, 
                                          pyne::rxname::_names+NUM_RX_NAMES);


std::map<std::string, unsigned int> pyne::rxname::altnames;
std::map<unsigned int, std::string> pyne::rxname::id_name;
std::map<std::string, unsigned int> pyne::rxname::name_id;
std::map<unsigned int, unsigned int> pyne::rxname::id_mt;
std::map<unsigned int, unsigned int> pyne::rxname::mt_id;
std::map<unsigned int, std::string> pyne::rxname::labels;
std::map<unsigned int, std::string> pyne::rxname::docs;
std::map<std::string, std::map<int, std::map<int, unsigned int> > > pyne::rxname::zadelta;

void * pyne::rxname::_fill_maps()
{
  std::string rx;
  unsigned int rxid;
  unsigned int _mts [NUM_RX_NAMES] = {
    1,
    2,
    3,
    4,
    5,
    10,
    11,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    27,
    28,
    29,
    30,
    32,
    33,
    34,
    35,
    36,
    37,
    38,
    41,
    42,
    44,
    45,
    50,
    51,
    52,
    53,
    54,
    55,
    56,
    57,
    58,
    59,
    60,
    61,
    62,
    63,
    64,
    65,
    66,
    67,
    68,
    69,
    70,
    71,
    72,
    73,
    74,
    75,
    76,
    77,
    78,
    79,
    80,
    81,
    82,
    83,
    84,
    85,
    86,
    87,
    88,
    89,
    90,
    91,
    101,
    102,
    103,
    104,
    105,
    106,
    107,
    108,
    109,
    111,
    112,
    113,
    114,
    115,
    116,
    117,
    151,
    201,
    202,
    203,
    204,
    205,
    206,
    207,
    208,
    209,
    210,
    211,
    212,
    213,
    214,
    215,
    216,
    217,
    218,
    251,
    252,
    253,
    444,
    451,
    452,
    454,
    455,
    456,
    457,
    458,
    459,
    460,
    500,
    501,
    502,
    504,
    505,
    506,
    515,
    516,
    517,
    522,
    523,
    526,
    527,
    528,
    533,
    534,
    535,
    536,
    537,
    538,
    539,
    540,
    541,
    542,
    543,
    544,
    545,
    546,
    547,
    548,
    549,
    550,
    551,
    552,
    553,
    554,
    555,
    556,
    557,
    558,
    559,
    560,
    561,
    562,
    563,
    564,
    565,
    566,
    567,
    568,
    569,
    570,
    571,
    572,
    600,
    601,
    602,
    603,
    604,
    605,
    606,
    607,
    608,
    609,
    610,
    611,
    612,
    613,
    614,
    615,
    616,
    617,
    618,
    619,
    620,
    621,
    622,
    623,
    624,
    625,
    626,
    627,
    628,
    629,
    630,
    631,
    632,
    633,
    634,
    635,
    636,
    637,
    638,
    639,
    640,
    641,
    642,
    643,
    644,
    645,
    646,
    647,
    648,
    649,
    650,
    651,
    652,
    653,
    654,
    655,
    656,
    657,
    658,
    659,
    660,
    661,
    662,
    663,
    664,
    665,
    666,
    667,
    668,
    669,
    670,
    671,
    672,
    673,
    674,
    675,
    676,
    677,
    678,
    679,
    680,
    681,
    682,
    683,
    684,
    685,
    686,
    687,
    688,
    689,
    690,
    691,
    692,
    693,
    694,
    695,
    696,
    697,
    698,
    699,
    749,
    750,
    751,
    752,
    753,
    754,
    755,
    756,
    757,
    758,
    759,
    760,
    761,
    762,
    763,
    764,
    765,
    766,
    767,
    768,
    769,
    770,
    771,
    772,
    773,
    774,
    775,
    776,
    777,
    778,
    779,
    780,
    781,
    782,
    783,
    784,
    785,
    786,
    787,
    788,
    789,
    790,
    791,
    792,
    793,
    794,
    795,
    796,
    797,
    798,
    799,
    849,
  };
  std::string _labels[NUM_RX_NAMES] = {
    "(z,total)",
    "(z,elastic)",
    "(z,nonelastic)",
    "(z,inelastic)",
    "(misc)",
    "(z,continuum)",
    "(z,2nd )",
    "(z,2n)",
    "(z,3n)",
    "(z,fission)",
    "(z,f)",
    "(z,nf)",
    "(z,2nf)",
    "(z,n+a)",
    "(z,n+3a)",
    "(z,2n+a)",
    "(z,3n+a)",
    "(z,abs) Absorption",
    "(z,n+p)",
    "(z,n+2a)",
    "(z,2n+2a)",
    "(z,nd)",
    "(z,nt)",
    "(z,n+He3)",
    "(z,n+d+3a)",
    "(z,n+t+2a)",
    "(z,4n)",
    "(z,3nf)",
    "(z,2n+p)",
    "(z,3n+p)",
    "(z,n+2p)",
    "(z,npa)",
    "(z,n0)",
    "(z,n1)",
    "(z,n2)",
    "(z,n3)",
    "(z,n4)",
    "(z,n5)",
    "(z,n6)",
    "(z,n7)",
    "(z,n8)",
    "(z,n9)",
    "(z,n10)",
    "(z,n11)",
    "(z,n12)",
    "(z,n13)",
    "(z,n14)",
    "(z,n15)",
    "(z,n16)",
    "(z,n17)",
    "(z,n18)",
    "(z,n19)",
    "(z,n20)",
    "(z,n21)",
    "(z,n22)",
    "(z,n23)",
    "(z,n24)",
    "(z,n25)",
    "(z,n26)",
    "(z,n27)",
    "(z,n28)",
    "(z,n29)",
    "(z,n30)",
    "(z,n31)",
    "(z,n32)",
    "(z,n33)",
    "(z,n34)",
    "(z,n35)",
    "(z,n36)",
    "(z,n37)",
    "(z,n38)",
    "(z,n39)",
    "(z,n40)",
    "(z,nc)",
    "(z,disap) Neutron disappeareance",
    "(z,gamma)",
    "(z,p)",
    "(z,d)",
    "(z,t)",
    "(z,3He)",
    "(z,a)",
    "(z,2a)",
    "(z,3a)",
    "(z,2p)",
    "(z,pa)",
    "(z,t2a)",
    "(z,d2a)",
    "(z,pd)",
    "(z,pt)",
    "(z,da)",
    "Resonance Parameters",
    "(z,Xn)",
    "(z,Xgamma)",
    "(z,Xp)",
    "(z,Xd)",
    "(z,Xt)",
    "(z,X3He)",
    "(z,Xa)",
    "(z,Xpi+) Total pi+ meson production",
    "(z,Xpi0) Total pi0 meson production",
    "(z,Xpi-) Total pi- meson production",
    "(z,Xmu+) Total anti-muon production",
    "(z,Xmu-) Total muon production",
    "(z,Xk+) Total positive kaon production",
    "(z,Xk0long) Total long-lived neutral kaon production",
    "(z,Xk0short) Total short-lived neutral kaon production",
    "(z,Xk-) Total negative kaon production",
    "(z,Xp-) Total anti-proton production",
    "(z,Xn-) Total anti-neutron production",
    "Average cosine of scattering angle",
    "Average logarithmic energy decrement",
    "Average xi^2/(2*xi)",
    "(damage)",
    "Desciptive Data",
    "Total Neutrons per Fission",
    "Independent fission product yield",
    "Delayed Neutron Data",
    "Prompt Neutrons per Fission",
    "Radioactive Decay Data",
    "Energy Release Due to Fission",
    "Cumulative Fission Product Yield",
    "Delayed Photon Data",
    "Total charged-particle stopping power",
    "Total photon interaction",
    "Photon coherent scattering",
    "Photon incoherent scattering",
    "Imaginary scattering factor",
    "Real scattering factor",
    "Pair production, electron field",
    "Total pair production",
    "Pair production, nuclear field",
    "Photoelectric absorption",
    "Photo-excitation cross section",
    "Electro-atomic scattering",
    "Electro-atomic bremsstrahlung",
    "Electro-atomic excitation cross section",
    "Atomic relaxation data",
    "K (1s1/2) subshell",
    "L1 (2s1/2) subshell",
    "L2 (2p1/2) subshell",
    "L3 (2p3/2) subshell",
    "M1 (3s1/2) subshell",
    "M2 (3p1/2) subshell",
    "M3 (3p3/2) subshell",
    "M4 (3d1/2) subshell",
    "M5 (3d1/2) subshell",
    "N1 (4s1/2) subshell",
    "N2 (4p1/2) subshell",
    "N3 (4p3/2) subshell",
    "N4 (4d3/2) subshell",
    "N5 (4d5/2) subshell",
    "N6 (4f5/2) subshell",
    "N7 (4f7/2) subshell",
    "O1 (5s1/2) subshell",
    "O2 (5p1/2) subshell",
    "O3 (5p3/2) subshell",
    "O4 (5d3/2) subshell",
    "O5 (5d5/2) subshell",
    "O6 (5f5/2) subshell",
    "O7 (5f7/2) subshell",
    "O8 (5g7/2) subshell",
    "O9 (5g9/2) subshell",
    "P1 (6s1/2) subshell",
    "P2 (6p1/2) subshell",
    "P3 (6p3/2) subshell",
    "P4 (6d3/2) subshell",
    "P5 (6d5/2) subshell",
    "P6 (6f5/2) subshell",
    "P7 (6f7/2) subshell",
    "P8 (6g7/2) subshell",
    "P9 (6g9/2) subshell",
    "P10 (6h9/2) subshell",
    "P11 (6h11/2) subshell",
    "Q1 (7s1/2) subshell",
    "Q2 (7p1/2) subshell",
    "Q3 (7p3/2) subshell",
    "(z,p0)",
    "(z,p1)",
    "(z,p2)",
    "(z,p3)",
    "(z,p4)",
    "(z,p5)",
    "(z,p6)",
    "(z,p7)",
    "(z,p8)",
    "(z,p9)",
    "(z,p10)",
    "(z,p11)",
    "(z,p12)",
    "(z,p13)",
    "(z,p14)",
    "(z,p15)",
    "(z,p16)",
    "(z,p17)",
    "(z,p18)",
    "(z,p19)",
    "(z,p20)",
    "(z,p21)",
    "(z,p22)",
    "(z,p23)",
    "(z,p24)",
    "(z,p25)",
    "(z,p26)",
    "(z,p27)",
    "(z,p28)",
    "(z,p29)",
    "(z,p30)",
    "(z,p31)",
    "(z,p32)",
    "(z,p33)",
    "(z,p34)",
    "(z,p35)",
    "(z,p36)",
    "(z,p37)",
    "(z,p38)",
    "(z,p39)",
    "(z,p40)",
    "(z,p41)",
    "(z,p42)",
    "(z,p43)",
    "(z,p44)",
    "(z,p45)",
    "(z,p46)",
    "(z,p47)",
    "(z,p48)",
    "(z,pc)",
    "(z,d0)",
    "(z,d1)",
    "(z,d2)",
    "(z,d3)",
    "(z,d4)",
    "(z,d5)",
    "(z,d6)",
    "(z,d7)",
    "(z,d8)",
    "(z,d9)",
    "(z,d10)",
    "(z,d11)",
    "(z,d12)",
    "(z,d13)",
    "(z,d14)",
    "(z,d15)",
    "(z,d16)",
    "(z,d17)",
    "(z,d18)",
    "(z,d19)",
    "(z,d20)",
    "(z,d21)",
    "(z,d22)",
    "(z,d23)",
    "(z,d24)",
    "(z,d25)",
    "(z,d26)",
    "(z,d27)",
    "(z,d28)",
    "(z,d29)",
    "(z,d30)",
    "(z,d31)",
    "(z,d32)",
    "(z,d33)",
    "(z,d34)",
    "(z,d35)",
    "(z,d36)",
    "(z,d37)",
    "(z,d38)",
    "(z,d39)",
    "(z,d40)",
    "(z,d41)",
    "(z,d42)",
    "(z,d43)",
    "(z,d44)",
    "(z,d45)",
    "(z,d46)",
    "(z,d47)",
    "(z,d48)",
    "(z,dc)",
    "(z,tc)",
    "(z,3He0)",
    "(z,3He1)",
    "(z,3He2)",
    "(z,3He3)",
    "(z,3He4)",
    "(z,3He5)",
    "(z,3He6)",
    "(z,3He7)",
    "(z,3He8)",
    "(z,3He9)",
    "(z,3He10)",
    "(z,3He11)",
    "(z,3He12)",
    "(z,3He13)",
    "(z,3He14)",
    "(z,3He15)",
    "(z,3He16)",
    "(z,3He17)",
    "(z,3He18)",
    "(z,3He19)",
    "(z,3He20)",
    "(z,3He21)",
    "(z,3He22)",
    "(z,3He23)",
    "(z,3He24)",
    "(z,3He25)",
    "(z,3He26)",
    "(z,3He27)",
    "(z,3He28)",
    "(z,3He29)",
    "(z,3He30)",
    "(z,3He31)",
    "(z,3He32)",
    "(z,3He33)",
    "(z,3He34)",
    "(z,3He35)",
    "(z,3He36)",
    "(z,3He37)",
    "(z,3He38)",
    "(z,3He39)",
    "(z,3He40)",
    "(z,3He41)",
    "(z,3He42)",
    "(z,3He43)",
    "(z,3He44)",
    "(z,3He45)",
    "(z,3He46)",
    "(z,3He47)",
    "(z,3He48)",
    "(z,3Hec)",
    "(z,ac)",
  };
  std::string _docs[NUM_RX_NAMES] = {
    "(n,total) Neutron total",
    "(z,z0) Elastic scattering",
    "(z,nonelas) Nonelastic neutron",
    "(z,n) One neutron in exit channel",
    "(z,anything) Miscellaneous",
    "(z,contin) Total continuum reaction",
    "(z,2nd) Production of 2n and d",
    "(z,2n) Production of 2n",
    "(z,3n) Production of 3n",
    "(z,fiss) Particle-induced fission",
    "(z,f) First-chance fission",
    "(z,nf) Second chance fission",
    "(z,2nf) Third-chance fission",
    "(z,na) Production of n and alpha",
    "(z,n3a) Production of n and 3 alphas",
    "(z,2na) Production of 2n and alpha",
    "(z,3na) Production of 3n and alpha",
    "(n,abs) Absorption",
    "(z,np) Production of n and p",
    "(z,n2a) Production of n and 2 alphas",
    "(z,2n2a) Production of 2n and 2 alphas",
    "(z,nd) Production of n and d",
    "(z,nt) Production of n and t",
    "(z,n3He) Production of n and He-3",
    "(z,nd2a) Production of n, d, and alpha",
    "(z,nt2a) Production of n, t, and 2 alphas",
    "(z,4n) Production of 4n",
    "(z,3nf) Fourth-chance fission",
    "(z,2np) Production of 2n and p",
    "(z,3np) Production of 3n and p",
    "(z,n2p) Production of n and 2p",
    "(z,npa) Production of n, p, and alpha",
    "(z,n0) Production of n, ground state",
    "(z,n1) Production of n, 1st excited state",
    "(z,n2) Production of n, 2nd excited state",
    "(z,n3) Production of n, 3rd excited state",
    "(z,n4) Production of n, 4th excited state",
    "(z,n5) Production of n, 5th excited state",
    "(z,n6) Production of n, 6th excited state",
    "(z,n7) Production of n, 7th excited state",
    "(z,n8) Production of n, 8th excited state",
    "(z,n9) Production of n, 9th excited state",
    "(z,n10) Production of n, 10th excited state",
    "(z,n11) Production of n, 11th excited state",
    "(z,n12) Production of n, 12th excited state",
    "(z,n13) Production of n, 13th excited state",
    "(z,n14) Production of n, 14th excited state",
    "(z,n15) Production of n, 15th excited state",
    "(z,n16) Production of n, 16th excited state",
    "(z,n17) Production of n, 17th excited state",
    "(z,n18) Production of n, 18th excited state",
    "(z,n19) Production of n, 19th excited state",
    "(z,n20) Production of n, 20th excited state",
    "(z,n21) Production of n, 21st excited state",
    "(z,n22) Production of n, 22nd excited state",
    "(z,n23) Production of n, 23rd excited state",
    "(z,n24) Production of n, 24th excited state",
    "(z,n25) Production of n, 25th excited state",
    "(z,n26) Production of n, 26th excited state",
    "(z,n27) Production of n, 27th excited state",
    "(z,n28) Production of n, 28th excited state",
    "(z,n29) Production of n, 29th excited state",
    "(z,n30) Production of n, 30th excited state",
    "(z,n31) Production of n, 31st excited state",
    "(z,n32) Production of n, 32nd excited state",
    "(z,n33) Production of n, 33rd excited state",
    "(z,n34) Production of n, 34th excited state",
    "(z,n35) Production of n, 35th excited state",
    "(z,n36) Production of n, 36th excited state",
    "(z,n37) Production of n, 37th excited state",
    "(z,n38) Production of n, 38th excited state",
    "(z,n39) Production of n, 39th excited state",
    "(z,n40) Production of n, 40th excited state",
    "(z,nc) Production of n in continuum",
    "(n,disap) Neutron disappeareance",
    "(z,gamma) Radiative capture",
    "(z,p) Production of p",
    "(z,d) Production of d",
    "(z,t) Production of t",
    "(z,3He) Production of He-3",
    "(z,a) Production of alpha",
    "(z,2a) Production of 2 alphas",
    "(z,3a) Production of 3 alphas",
    "(z,2p) Production of 2p",
    "(z,pa) Production of p and alpha",
    "(z,t2a) Production of t and 2 alphas",
    "(z,d2a) Production of d and 2 alphas",
    "(z,pd) Production of p and d",
    "(z,pt) Production of p and t",
    "(z,da) Production of d and a",
    "Resonance Parameters",
    "(z,Xn) Total neutron production",
    "(z,Xgamma) Total gamma production",
    "(z,Xp) Total proton production",
    "(z,Xd) Total deuteron production",
    "(z,Xt) Total triton production",
    "(z,X3He) Total He-3 production",
    "(z,Xa) Total alpha production",
    "(z,Xpi+) Total pi+ meson production",
    "(z,Xpi0) Total pi0 meson production",
    "(z,Xpi-) Total pi- meson production",
    "(z,Xmu+) Total anti-muon production",
    "(z,Xmu-) Total muon production",
    "(z,Xk+) Total positive kaon production",
    "(z,Xk0long) Total long-lived neutral kaon production",
    "(z,Xk0short) Total short-lived neutral kaon production",
    "(z,Xk-) Total negative kaon production",
    "(z,Xp-) Total anti-proton production",
    "(z,Xn-) Total anti-neutron production",
    "Average cosine of scattering angle",
    "Average logarithmic energy decrement",
    "Average xi^2/(2*xi)",
    "(damage)",
    "Desciptive Data",
    "Total Neutrons per Fission",
    "Independent fission product yield",
    "Delayed Neutron Data",
    "Prompt Neutrons per Fission",
    "Radioactive Decay Data",
    "Energy Release Due to Fission",
    "Cumulative Fission Product Yield",
    "Delayed Photon Data",
    "Total charged-particle stopping power",
    "Total photon interaction",
    "Photon coherent scattering",
    "Photon incoherent scattering",
    "Imaginary scattering factor",
    "Real scattering factor",
    "Pair production, electron field",
    "Total pair production",
    "Pair production, nuclear field",
    "Photoelectric absorption",
    "Photo-excitation cross section",
    "Electro-atomic scattering",
    "Electro-atomic bremsstrahlung",
    "Electro-atomic excitation cross section",
    "Atomic relaxation data",
    "K (1s1/2) subshell",
    "L1 (2s1/2) subshell",
    "L2 (2p1/2) subshell",
    "L3 (2p3/2) subshell",
    "M1 (3s1/2) subshell",
    "M2 (3p1/2) subshell",
    "M3 (3p3/2) subshell",
    "M4 (3d1/2) subshell",
    "M5 (3d1/2) subshell",
    "N1 (4s1/2) subshell",
    "N2 (4p1/2) subshell",
    "N3 (4p3/2) subshell",
    "N4 (4d3/2) subshell",
    "N5 (4d5/2) subshell",
    "N6 (4f5/2) subshell",
    "N7 (4f7/2) subshell",
    "O1 (5s1/2) subshell",
    "O2 (5p1/2) subshell",
    "O3 (5p3/2) subshell",
    "O4 (5d3/2) subshell",
    "O5 (5d5/2) subshell",
    "O6 (5f5/2) subshell",
    "O7 (5f7/2) subshell",
    "O8 (5g7/2) subshell",
    "O9 (5g9/2) subshell",
    "P1 (6s1/2) subshell",
    "P2 (6p1/2) subshell",
    "P3 (6p3/2) subshell",
    "P4 (6d3/2) subshell",
    "P5 (6d5/2) subshell",
    "P6 (6f5/2) subshell",
    "P7 (6f7/2) subshell",
    "P8 (6g7/2) subshell",
    "P9 (6g9/2) subshell",
    "P10 (6h9/2) subshell",
    "P11 (6h11/2) subshell",
    "Q1 (7s1/2) subshell",
    "Q2 (7p1/2) subshell",
    "Q3 (7p3/2) subshell",
    "(n,p0)",
    "(n,p1)",
    "(n,p2)",
    "(n,p3)",
    "(n,p4)",
    "(n,p5)",
    "(n,p6)",
    "(n,p7)",
    "(n,p8)",
    "(n,p9)",
    "(n,p10)",
    "(n,p11)",
    "(n,p12)",
    "(n,p13)",
    "(n,p14)",
    "(n,p15)",
    "(n,p16)",
    "(n,p17)",
    "(n,p18)",
    "(n,p19)",
    "(n,p20)",
    "(n,p21)",
    "(n,p22)",
    "(n,p23)",
    "(n,p24)",
    "(n,p25)",
    "(n,p26)",
    "(n,p27)",
    "(n,p28)",
    "(n,p29)",
    "(n,p30)",
    "(n,p31)",
    "(n,p32)",
    "(n,p33)",
    "(n,p34)",
    "(n,p35)",
    "(n,p36)",
    "(n,p37)",
    "(n,p38)",
    "(n,p39)",
    "(n,p40)",
    "(n,p41)",
    "(n,p42)",
    "(n,p43)",
    "(n,p44)",
    "(n,p45)",
    "(n,p46)",
    "(n,p47)",
    "(n,p48)",
    "(n,pc)",
    "(n,d0)",
    "(n,d1)",
    "(n,d2)",
    "(n,d3)",
    "(n,d4)",
    "(n,d5)",
    "(n,d6)",
    "(n,d7)",
    "(n,d8)",
    "(n,d9)",
    "(n,d10)",
    "(n,d11)",
    "(n,d12)",
    "(n,d13)",
    "(n,d14)",
    "(n,d15)",
    "(n,d16)",
    "(n,d17)",
    "(n,d18)",
    "(n,d19)",
    "(n,d20)",
    "(n,d21)",
    "(n,d22)",
    "(n,d23)",
    "(n,d24)",
    "(n,d25)",
    "(n,d26)",
    "(n,d27)",
    "(n,d28)",
    "(n,d29)",
    "(n,d30)",
    "(n,d31)",
    "(n,d32)",
    "(n,d33)",
    "(n,d34)",
    "(n,d35)",
    "(n,d36)",
    "(n,d37)",
    "(n,d38)",
    "(n,d39)",
    "(n,d40)",
    "(n,d41)",
    "(n,d42)",
    "(n,d43)",
    "(n,d44)",
    "(n,d45)",
    "(n,d46)",
    "(n,d47)",
    "(n,d48)",
    "(n,dc)",
    "(n,tc)",
    "(n,3He0)",
    "(n,3He1)",
    "(n,3He2)",
    "(n,3He3)",
    "(n,3He4)",
    "(n,3He5)",
    "(n,3He6)",
    "(n,3He7)",
    "(n,3He8)",
    "(n,3He9)",
    "(n,3He10)",
    "(n,3He11)",
    "(n,3He12)",
    "(n,3He13)",
    "(n,3He14)",
    "(n,3He15)",
    "(n,3He16)",
    "(n,3He17)",
    "(n,3He18)",
    "(n,3He19)",
    "(n,3He20)",
    "(n,3He21)",
    "(n,3He22)",
    "(n,3He23)",
    "(n,3He24)",
    "(n,3He25)",
    "(n,3He26)",
    "(n,3He27)",
    "(n,3He28)",
    "(n,3He29)",
    "(n,3He30)",
    "(n,3He31)",
    "(n,3He32)",
    "(n,3He33)",
    "(n,3He34)",
    "(n,3He35)",
    "(n,3He36)",
    "(n,3He37)",
    "(n,3He38)",
    "(n,3He39)",
    "(n,3He40)",
    "(n,3He41)",
    "(n,3He42)",
    "(n,3He43)",
    "(n,3He44)",
    "(n,3He45)",
    "(n,3He46)",
    "(n,3He47)",
    "(n,3He48)",
    "(n,3Hec)",
    "(n,ac)",
  };

  // fill the maps
  for (int i = 0; i < NUM_RX_NAMES; i++)
  {
    rx = _names[i];
    rxid = pyne::rxname::hash(rx);
    id_name[rxid] = rx;
    name_id[rx] = rxid;
    if (0 < _mts[i]) {
      id_mt[rxid] = _mts[i];
      mt_id[_mts[i]] = rxid;
    };
    labels[rxid] = _labels[i];
    docs[rxid] = _docs[i];
  };

  // set alternative names
  altnames["tot"] = name_id["total"];
  altnames["abs"] = name_id["absorption"];
  altnames["fis"] = name_id["fission"];
  altnames["fiss"] = name_id["fission"];
  altnames["alpha"] = name_id["a"];
  altnames["deut"] = name_id["d"];
  altnames["deuteron"] = name_id["d"];
  altnames["deuterium"] = name_id["d"];
  altnames["trit"] = name_id["t"];
  altnames["triton"] = name_id["t"];
  altnames["tritium"] = name_id["t"];
  altnames["proton"] = name_id["p"];
  altnames["he3"] = name_id["He3"];
  altnames["HE3"] = name_id["He3"];
  altnames["3HE"] = name_id["He3"];
  altnames["3He"] = name_id["He3"];
  altnames["3he"] = name_id["He3"];
  altnames["he-3"] = name_id["He3"];
  altnames["HE-3"] = name_id["He3"];

  // set the nulcide difference mappings, zadelta
  // zadelta[incident particle type "n", "p", ...][delta Z num][delta A num][rxid]
  zadelta["n"][0][1] = name_id["absorption"];
  zadelta["n"][-1][0] = name_id["p"];
  zadelta["n"][-1][-1] = name_id["d"];
  zadelta["p"][1][1] = name_id["absorption"];

};
void * pyne::rxname::_ = pyne::rxname::_fill_maps();



unsigned int pyne::rxname::hash(std::string s){return pyne::rxname::hash(s.c_str());};
unsigned int pyne::rxname::hash(const char * s)
{
  // Modified from http://cboard.cprogramming.com/tech-board/114650-string-hashing-algorithm.html#post853145
  // starting from h = 32*2^5 > 1000, rather than 0, to reserve space for MT numbers
  int c;
  unsigned int h = 32; 
  while((c = *s++))
  {
    h = ((h << 5) + h) ^ c;
  }
  return h;
};


// ************************
// *** name functions *****
// ************************

std::string pyne::rxname::name(char * s){return pyne::rxname::name(std::string(s));};
std::string pyne::rxname::name(std::string s)
{
  if (0 < names.count(s))
    return s;
  if (0 < altnames.count(s))
    return id_name[altnames[s]];
  // see if id in string form
  int i = 0;
  int I = s.length();
  int found = 0;
  while(0 <= found && i < I)
  {
    found = pyne::digits.find(s[i]);
    i++;
  }
  if (0<=found)
    return pyne::rxname::name(atoi(s.c_str()));
  // dead...
  throw NotAReaction(s, "???");
};


std::string pyne::rxname::name(int n){return pyne::rxname::name((unsigned int) n);};
std::string pyne::rxname::name(unsigned int n)
{
  if (0 < id_name.count(n))
    return id_name[n];
  if (0 < mt_id.count(n))
    return id_name[mt_id[n]];
  throw NotAReaction(n, "???");
};


std::string pyne::rxname::name(int from_nuc, int to_nuc, std::string z)
{
  // This assumes nuclides are in zzaaam form
  int dz = (to_nuc/10000) - (from_nuc/10000);
  int da = ((to_nuc/10)%1000) - ((from_nuc/10)%1000);
  if (0 == zadelta.count(z))
    throw IndeterminateReactionForm("z=" + z, "???");
  if (0 == zadelta[z].count(dz))
    throw IndeterminateReactionForm(pyne::to_str(from_nuc) + ", " + \
                                    pyne::to_str(to_nuc), "delta Z = " + \
                                    pyne::to_str(dz));
  if (0 == zadelta[z][dz].count(da))
    throw IndeterminateReactionForm(pyne::to_str(from_nuc) + ", " + \
                                    pyne::to_str(to_nuc), "delta A = " + \
                                    pyne::to_str(da));
  return id_name[zadelta[z][dz][da]];
};

std::string pyne::rxname::name(std::string from_nuc, int to_nuc, std::string z)
{
  return pyne::rxname::name(pyne::nucname::zzaaam(from_nuc), 
                            pyne::nucname::zzaaam(to_nuc), z);
};

std::string pyne::rxname::name(int from_nuc, std::string to_nuc, std::string z)
{
  return pyne::rxname::name(pyne::nucname::zzaaam(from_nuc), 
                            pyne::nucname::zzaaam(to_nuc), z);
};

std::string pyne::rxname::name(std::string from_nuc, std::string to_nuc, std::string z)
{
  return pyne::rxname::name(pyne::nucname::zzaaam(from_nuc), 
                            pyne::nucname::zzaaam(to_nuc), z);
};



// **********************
// *** id functions *****
// **********************
unsigned int pyne::rxname::id(int x)
{
  return name_id[pyne::rxname::name(x)];
};
  
unsigned int pyne::rxname::id(unsigned int x)
{
  if (0 < id_name.count(x))
    return x;
  if (0 < mt_id.count(x))
    return mt_id[x];
  return name_id[pyne::rxname::name(x)];
};
  
unsigned int pyne::rxname::id(char * x)
{
  return name_id[pyne::rxname::name(x)];
};
  
unsigned int pyne::rxname::id(std::string x)
{
  if (0 < names.count(x))
    return name_id[x];
  if (0 < altnames.count(x))
    return altnames[x];
  return name_id[pyne::rxname::name(x)];  
};
  
unsigned int pyne::rxname::id(int from_nuc, int to_nuc, std::string z)
{
  // This assumes nuclides are in zzaaam form
  int dz = (to_nuc/10000) - (from_nuc/10000);
  int da = ((to_nuc/10)%1000) - ((from_nuc/10)%1000);
  if (0 == zadelta.count(z))
    throw IndeterminateReactionForm("z=" + z, "???");
  if (0 == zadelta[z].count(dz))
    throw IndeterminateReactionForm(pyne::to_str(from_nuc) + ", " + \
                                    pyne::to_str(to_nuc), "delta Z = " + \
                                    pyne::to_str(dz));
  if (0 == zadelta[z][dz].count(da))
    throw IndeterminateReactionForm(pyne::to_str(from_nuc) + ", " + \
                                    pyne::to_str(to_nuc), "delta A = " + \
                                    pyne::to_str(da));
  return zadelta[z][dz][da];
};
  
unsigned int pyne::rxname::id(int from_nuc, std::string to_nuc, std::string z)
{
  return pyne::rxname::id(pyne::nucname::zzaaam(from_nuc), 
                          pyne::nucname::zzaaam(to_nuc), z);
};
  
unsigned int pyne::rxname::id(std::string from_nuc, int to_nuc, std::string z)
{
  return pyne::rxname::id(pyne::nucname::zzaaam(from_nuc), 
                          pyne::nucname::zzaaam(to_nuc), z);
};
  
unsigned int pyne::rxname::id(std::string from_nuc, std::string to_nuc, std::string z)
{
  return pyne::rxname::id(pyne::nucname::zzaaam(from_nuc), 
                          pyne::nucname::zzaaam(to_nuc), z);
};


// **********************
// *** MT functions *****
// **********************
unsigned int pyne::rxname::mt(int x)
{
  unsigned int rxid = pyne::rxname::id(x);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
};
  
unsigned int pyne::rxname::mt(unsigned int x)
{
  unsigned int rxid = pyne::rxname::id(x);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
};
  
unsigned int pyne::rxname::mt(char * x)
{
  unsigned int rxid = pyne::rxname::id(x);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
};
  
unsigned int pyne::rxname::mt(std::string x)
{
  unsigned int rxid = pyne::rxname::id(x);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
};
  
unsigned int pyne::rxname::mt(int from_nuc, int to_nuc, std::string z)
{
  unsigned int rxid = pyne::rxname::id(from_nuc, to_nuc, z);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
};
  
unsigned int pyne::rxname::mt(int from_nuc, std::string to_nuc, std::string z)
{
  unsigned int rxid = pyne::rxname::id(from_nuc, to_nuc, z);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
};
  
unsigned int pyne::rxname::mt(std::string from_nuc, int to_nuc, std::string z)
{
  unsigned int rxid = pyne::rxname::id(from_nuc, to_nuc, z);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
};
  
unsigned int pyne::rxname::mt(std::string from_nuc, std::string to_nuc, std::string z)
{
  unsigned int rxid = pyne::rxname::id(from_nuc, to_nuc, z);
  if (0 == id_mt.count(rxid))
    throw NotAReaction();
  return id_mt[rxid];
};


// *************************
// *** label functions *****
// *************************
std::string pyne::rxname::label(int x)
{
  return labels[pyne::rxname::id(x)];
};
  
std::string pyne::rxname::label(unsigned int x)
{
  return labels[pyne::rxname::id(x)];
};
  
std::string pyne::rxname::label(char * x)
{
  return labels[pyne::rxname::id(x)];
};
  
std::string pyne::rxname::label(std::string x)
{
  return labels[pyne::rxname::id(x)];
};
  
std::string pyne::rxname::label(int from_nuc, int to_nuc, std::string z)
{
  return labels[pyne::rxname::id(from_nuc, to_nuc, z)];
};
  
std::string pyne::rxname::label(int from_nuc, std::string to_nuc, std::string z)
{
  return labels[pyne::rxname::id(from_nuc, to_nuc, z)];
};
  
std::string pyne::rxname::label(std::string from_nuc, int to_nuc, std::string z)
{
  return labels[pyne::rxname::id(from_nuc, to_nuc, z)];
};
  
std::string pyne::rxname::label(std::string from_nuc, std::string to_nuc, std::string z)
{
  return labels[pyne::rxname::id(from_nuc, to_nuc, z)];
};
