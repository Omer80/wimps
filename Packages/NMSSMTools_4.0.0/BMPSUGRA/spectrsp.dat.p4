# NMSSMTools OUTPUT IN SLHA FORMAT
# Info about spectrum calculator
BLOCK SPINFO   # Program information
     1   NMSSMTools # Spectrum calculator
     2   4.0.0      # Version number
     8   0          # Higgs mass precision
     3   # Relic density too small (WMAP)
     3   # Muon magn. mom. more than 2 sigma away
# Input parameters
BLOCK MODSEL
  3  1   # NMSSM particle content
  1  1   # IMOD
 10  0   # ISCAN
  9  1   # Call micrOmegas
  8  0   # Precision for Higgs masses
 13  0   # Sparticle decays via NMSDECAY
BLOCK SMINPUTS
     1     1.27920000E+02   # ALPHA_EM^-1(MZ)
     2     1.16639000E-05   # GF
     3     1.17200000E-01   # ALPHA_S(MZ)
     4     9.11870000E+01   # MZ
     5     4.21400000E+00   # MB(MB)
     6     1.71400000E+02   # MTOP (POLE MASS)
     7     1.77700000E+00   # MTAU
# SMINPUTS Beyond SLHA:
# MW:     0.80420000E+02
# MS:     0.19000000E+00
# MC:     0.14000000E+01
# VUS:     0.22000000E+00
# VCB:     0.40000000E-01
# VUB:     0.40000000E-02
BLOCK MINPAR
     1     7.85100000E+02   # M0(MGUT)
     2     7.73000000E+02   # M12(MGUT)
     3     3.05000000E+00   # TANBETA(MZ)
     4    -1.00000000E+00   # SIGMU
     5    -2.25000000E+03   # A0(MGUT)
BLOCK EXTPAR
    21     7.25000000E+05   # MHD^2 AT THE GUT SCALE
    22     4.65000000E+06   # MHU^2 AT THE GUT SCALE
    61     4.90000000E-01   # LAMBDA AT THE SUSY SCALE
#   62     5.55501507E-02   # KAPPA AT THE SUSY SCALE
    64    -9.65000000E+02   # AKAPPA AT THE GUT SCALE
#   65    -2.07072729E+02   # MUEFF AT THE SUSY SCALE
#   70     1.41723919E+06   # MS^2 AT THE GUT SCALE
# 
BLOCK MASS   # Mass spectrum 
#  PDG Code     mass             particle 
        25     4.05478284E+01   # lightest neutral scalar
        35     1.20687645E+02   # second neutral scalar
        45     6.82335413E+02   # third neutral scalar
        36     1.14761398E+02   # lightest pseudoscalar
        46     6.82062899E+02   # second pseudoscalar
        37     6.77483312E+02   # charged Higgs
   1000001     1.74112125E+03   #  ~d_L
   2000001     1.67301586E+03   #  ~d_R
   1000002     1.73971100E+03   #  ~u_L
   2000002     1.73584790E+03   #  ~u_R
   1000003     1.74112125E+03   #  ~s_L
   2000003     1.67301586E+03   #  ~s_R
   1000004     1.73971100E+03   #  ~c_L
   2000004     1.73584790E+03   #  ~c_R
   1000005     1.26244506E+03   #  ~b_1
   2000005     1.62431184E+03   #  ~b_2
   1000006     4.10807083E+02   #  ~t_1
   2000006     1.28604831E+03   #  ~t_2
   1000011     9.84157355E+02   #  ~e_L
   2000011     7.08531150E+02   #  ~e_R
   1000012     9.81654479E+02   #  ~nue_L
   1000013     9.84157355E+02   #  ~mu_L
   2000013     7.08531150E+02   #  ~mu_R
   1000014     9.81654479E+02   #  ~numu_L
   1000015     7.05330000E+02   #  ~tau_1
   2000015     9.83035912E+02   #  ~tau_2
   1000016     9.80516187E+02   #  ~nutau_L
   1000021     1.75511624E+03   #  ~g
   1000022    -6.07024641E+01   # neutralino(1)
   1000023     2.25299854E+02   # neutralino(2)
   1000025    -2.26963493E+02   # neutralino(3)
   1000035     3.33994199E+02   # neutralino(4)
   1000045     6.38987227E+02   # neutralino(5)
   1000024     2.12648951E+02   # chargino(1)
   1000037    -6.38961355E+02   # chargino(2)
# 
# Low energy observables
BLOCK LOWEN
# Exp. 2 Sigma: 3.04E-4 < BR(b -> s gamma) < 4.06E-4:
     1     3.90035573E-04   # BR(b -> s gamma)
    11     4.35921263E-04   # (BR(b -> s gamma)+Theor.Err.)
    12     3.24827761E-04   # (BR(b -> s gamma)-Theor.Err.)
# Exp. 2 Sigma: 4.99E-1 < Delta M_d < 5.15E-1:
     2     6.29647995E-01   # Delta M_d in ps^-1
    21     1.10787566E+00   # Delta M_d +Theor.Err.
    22     1.67614180E-01   # Delta M_d -Theor.Err.
# Exp. 2 Sigma: 1.7633E+1 < Delta Ms < 1.7805E+1:
     3     2.18144808E+01   # Delta M_s in ps^-1
    31     2.90038912E+01   # Delta M_s +Theor.Err.
    32     1.48609705E+01   # Delta M_s -Theor.Err.
# Exp. 2 Sigma: 2.0E-9 < BR(Bs->mu+mu-) < 4.7E-9:
     4     3.54041393E-09   # BR(Bs -> mu+mu-)
    41     6.01332642E-09   # BR(Bs -> mu+mu-)+Theor.Err.
    42     1.71863451E-09   # BR(Bs -> mu+mu-)-Theor.Err.
# Exp. 2 Sigma: 1.07E-4 < BR(B+ > tau+ + nu_tau) < 2.27E-4:
     5     1.31651752E-04   # BR(B+ -> tau+ + nu_tau)
    51     2.63355000E-04   # BR(B+ -> tau+ + nu_tau) + Theor.Err.
    52     5.68177368E-05   # BR(B+ -> tau+ + nu_tau) - Theor.Err.
# 
# BSM contr. to the muon anomalous magn. moment:
     6    -8.46772014E-11   # Del_a_mu
    61     1.97137600E-10   # Del_a_mu + Theor.Err.
    62    -3.66492002E-10   # Del_a_mu - Theor.Err.
# Minimal Exp.-SM (2 sigma):  8.77306222E-10
# Maximal Exp.-SM (2 sigma):  4.61144414E-09
# 
# Omega h^2 (allowed: 0.110 < Omega h^2 < 0.127):
    10     8.46354354E-02   # Omega h^2
# Channels which contribute to 1/(omega) more than 0.1%.
# Relative contributions in % are displayed
  50.7% ~o1 ~o1 ->h1 h1
  16.2% ~o1 ~o1 ->b B
   5.4% ~o1 ~o1 ->d D
   5.4% ~o1 ~o1 ->s S
   4.7% ~o1 ~o1 ->c C
   4.2% ~o1 ~o1 ->u U
   2.4% ~o1 ~o1 ->ne Ne
   2.4% ~o1 ~o1 ->nm Nm
   2.4% ~o1 ~o1 ->nl Nl
   2.4% ~o1 ~o1 ->l L
   1.3% ~o1 ~o1 ->G G
   1.2% ~o1 ~o1 ->e E
   1.2% ~o1 ~o1 ->m M
# 
BLOCK HMIX Q=  8.12688087E+02 # (STOP/SBOTTOM MASSES)
     1    -2.05609588E+02   # MUEFF
     2     3.04133642E+00   # TAN(BETA)
     3     2.42978442E+02   # V(Q)
     4     4.64030810E+05   # MA^2
     5     1.86917053E+04   # MP^2
# 
# 3*3 Higgs mixing
BLOCK NMHMIX
  1  1    -1.02148240E-01   # S_(1,1)
  1  2     1.13161264E-02   # S_(1,2)
  1  3     9.94704822E-01   # S_(1,3)
  2  1     3.15673364E-01   # S_(2,1)
  2  2     9.48621460E-01   # S_(2,2)
  2  3     2.16252677E-02   # S_(2,3)
  3  1     9.43353627E-01   # S_(3,1)
  3  2    -3.16210800E-01   # S_(3,2)
  3  3     1.00472212E-01   # S_(3,3)
# 
# 3*3 Pseudoscalar Higgs mixing
BLOCK NMAMIX
  1  1     1.04435869E-01   # P_(1,1)
  1  2     3.42412687E-02   # P_(1,2)
  1  3     9.93941993E-01   # P_(1,3)
  2  1     9.44473043E-01   # P_(2,1)
  2  2     3.09663293E-01   # P_(2,2)
  2  3    -1.09905938E-01   # P_(2,3)
# 
# 3rd generation sfermion mixing
BLOCK STOPMIX  # Stop mixing matrix
  1  1     1.71681820E-01   # Rst_(1,1)
  1  2     9.85152452E-01   # Rst_(1,2)
  2  1    -9.85152452E-01   # Rst_(2,1)
  2  2     1.71681820E-01   # Rst_(2,2)
BLOCK SBOTMIX  # Sbottom mixing matrix
  1  1     9.99972399E-01   # Rsb_(1,1)
  1  2     7.42978843E-03   # Rsb_(1,2)
  2  1    -7.42978843E-03   # Rsb_(2,1)
  2  2     9.99972399E-01   # Rsb_(2,2)
BLOCK STAUMIX  # Stau mixing matrix
  1  1     7.19460171E-03   # Rsl_(1,1)
  1  2     9.99974119E-01   # Rsl_(1,2)
  2  1    -9.99974119E-01   # Rsl_(2,1)
  2  2     7.19460171E-03   # Rsl_(2,2)
# 
# Gaugino-Higgsino mixing
BLOCK NMNMIX  # 5*5 Neutralino Mixing Matrix
  1  1    -3.47802641E-02   # N_(1,1)
  1  2     3.66699123E-02   # N_(1,2)
  1  3     3.92331920E-02   # N_(1,3)
  1  4     3.44188594E-01   # N_(1,4)
  1  5     9.36718004E-01   # N_(1,5)
  2  1    -1.59135686E-01   # N_(2,1)
  2  2     7.58942998E-02   # N_(2,2)
  2  3     7.03956361E-01   # N_(2,3)
  2  4     6.32492105E-01   # N_(2,4)
  2  5    -2.70767553E-01   # N_(2,5)
  3  1    -6.59394386E-02   # N_(3,1)
  3  2     7.99712937E-02   # N_(3,2)
  3  3    -7.05900717E-01   # N_(3,3)
  3  4     6.65118331E-01   # N_(3,4)
  3  5    -2.20405005E-01   # N_(3,5)
  4  1     9.84311670E-01   # N_(4,1)
  4  2     3.48260565E-02   # N_(4,2)
  4  3     6.78952300E-02   # N_(4,3)
  4  4     1.57049310E-01   # N_(4,4)
  4  5    -2.53659868E-02   # N_(4,5)
  5  1    -1.57699919E-02   # N_(5,1)
  5  2     9.92616296E-01   # N_(5,2)
  5  3    -7.83467323E-04   # N_(5,3)
  5  4    -1.20171000E-01   # N_(5,4)
  5  5     4.74485197E-03   # N_(5,5)
# 
BLOCK UMIX  # Chargino U Mixing Matrix
  1  1     9.52461143E-04   # U_(1,1)
  1  2     9.99999546E-01   # U_(1,2)
  2  1    -9.99999546E-01   # U_(2,1)
  2  2     9.52461143E-04   # U_(2,2)
# 
BLOCK VMIX  # Chargino V Mixing Matrix
  1  1     1.69447287E-01   # V_(1,1)
  1  2    -9.85539252E-01   # V_(1,2)
  2  1     9.85539252E-01   # V_(2,1)
  2  2     1.69447287E-01   # V_(2,2)
# 
# Higgs reduced couplings
# (as compared to a SM Higgs with same mass)
BLOCK REDCOUP
# H1
  1  1     1.19088345E-02   # U-type fermions
  1  2    -3.27870390E-01   # D-type fermions
  1  3    -2.10714348E-02   # W,Z bosons
  1  4     1.70596504E-01   # Gluons
  1  5     6.70848639E-02   # Photons
# H2
  2  1     9.98307693E-01   # U-type fermions
  2  2     1.01323282E+00   # D-type fermions
  2  3     9.99756383E-01   # W,Z bosons
  2  4     9.73278343E-01   # Gluons
  2  5     9.95684275E-01   # Photons
# H3
  3  1    -3.32773068E-01   # U-type fermions
  3  2     3.02793003E+00   # D-type fermions
  3  3    -6.57038965E-03   # W,Z bosons
  3  4     3.18737366E-01   # Gluons
  3  5     1.37688282E+00   # Photons
# A1
  4  1     3.60347339E-02   # U-type fermions
  4  2     3.35213112E-01   # D-type fermions
  4  3     0.00000000E+00   # W,Z bosons
  4  4     3.20495337E-02   # Gluons
  4  5     2.77151439E-01   # Photons
# A2
  5  1     3.25882621E-01   # U-type fermions
  5  2     3.03152308E+00   # D-type fermions
  5  3     0.00000000E+00   # W,Z bosons
  5  4     3.32424156E-01   # Gluons
  5  5     2.86878692E-01   # Photons
# 
# GAUGE AND YUKAWA COUPLINGS AT THE SUSY SCALE
BLOCK GAUGE Q=  1.65984760E+03 # (SUSY SCALE)
     1     3.64125773E-01   # g1(Q,DR_bar)
     2     6.42237552E-01   # g2(Q,DR_bar)
     3     1.03603577E+00   # g3(Q,DR_bar)
BLOCK YU Q=  1.65984760E+03 # (SUSY SCALE)
  3  3     8.82775777E-01   # HTOP(Q,DR_bar)
BLOCK YD Q=  1.65984760E+03 # (SUSY SCALE)
  3  3     4.38617192E-02   # HBOT(Q,DR_bar)
BLOCK YE Q=  1.65984760E+03 # (SUSY SCALE)
  3  3     3.18978317E-02   # HTAU(Q,DR_bar)
# 
# SOFT TRILINEAR COUPLINGS AT THE SUSY SCALE
BLOCK AU Q=  1.65984760E+03 # (SUSY SCALE)
  3  3    -1.76253201E+03   # ATOP
BLOCK AD Q=  1.65984760E+03 # (SUSY SCALE)
  3  3    -3.72584147E+03   # ABOT
BLOCK AE Q=  1.65984760E+03 # (SUSY SCALE)
  2  2    -2.53008870E+03   # AMUON
  3  3    -2.52360581E+03   # ATAU
# 
# SOFT MASSES AT THE SUSY SCALE
BLOCK MSOFT Q=  1.65984760E+03 # (SUSY SCALE)
     1     3.32925617E+02   # M1
     2     6.12323801E+02   # M2
     3     1.66620866E+03   # M3
    21     3.74590750E+05   # M_HD^2
    22     8.10124535E+04   # M_HU^2
    31     9.83284114E+02   # M_eL
    32     9.83284114E+02   # M_muL
    33     9.82147711E+02   # M_tauL
    34     7.07483079E+02   # M_eR
    35     7.07483079E+02   # M_muR
    36     7.04292154E+02   # M_tauR
    41     1.67813712E+03   # M_q1L
    42     1.67813712E+03   # M_q2L
    43     1.25632496E+03   # M_q3L
    44     1.67375892E+03   # M_uR
    45     1.67375892E+03   # M_cR
    46     5.25709467E+02   # M_tR
    47     1.60829688E+03   # M_dR
    48     1.60829688E+03   # M_sR
    49     1.60588059E+03   # M_bR
# 
# NMSSM SPECIFIC PARAMETERS THE SUSY SCALE
BLOCK NMSSMRUN Q=  1.65984760E+03 # (SUSY SCALE)
     1     4.90000000E-01   # LAMBDA(Q,DR_bar)
     2     5.55501507E-02   # KAPPA(Q,DR_bar)
     3    -6.84453379E+02   # ALAMBDA
     4     1.51737093E+02   # AKAPPA
     5    -2.07072729E+02   # MUEFF
     6     0.00000000E+00   # XIF
     7     0.00000000E+00   # XIS
     8     0.00000000E+00   # MU'
     9     0.00000000E+00   # MS'^2
    10     6.09537505E+03   # MS^2
    12     0.00000000E+00   # M3H^2
# 
# GAUGE AND YUKAWA COUPLINGS AT THE GUT SCALE
BLOCK GAUGE Q=  1.74778596E+16 # (GUT SCALE)
     1     7.09796582E-01   # g1(MGUT,DR_bar), GUT normalization
     2     7.09796580E-01   # g2(MGUT,DR_bar)
     3     6.99961387E-01   # g3(MGUT,DR_bar)
BLOCK YU Q=  1.74778596E+16 # (GUT SCALE)
  3  3     6.46683857E-01   # HTOP(MGUT,DR_bar)
BLOCK YD Q=  1.74778596E+16 # (GUT SCALE)
  3  3     1.91120641E-02   # HBOT(MGUT,DR_bar)
BLOCK YE Q=  1.74778596E+16 # (GUT SCALE)
  3  3     2.35541502E-02   # HTAU(MGUT,DR_bar)
# 
# SOFT TRILINEAR COUPLINGS AT THE GUT SCALE
BLOCK AU Q=  1.74778596E+16 # (GUT SCALE)
  3  3    -2.25000041E+03   # ATOP
BLOCK AD Q=  1.74778596E+16 # (GUT SCALE)
  3  3    -2.24999958E+03   # ABOT
BLOCK AE Q=  1.74778596E+16 # (GUT SCALE)
  2  2    -2.24999952E+03   # AMUON
  3  3    -2.24999951E+03   # ATAU
# 
# SOFT MASSES SQUARED AT THE GUT SCALE
BLOCK MSOFT Q=  1.74778596E+16 # (GUT SCALE)
     1     7.73000128E+02   # M1
     2     7.73000051E+02   # M2
     3     7.72999931E+02   # M3
    21     7.24996947E+05   # M_HD^2
    22     4.64999906E+06   # M_HU^2
    31     7.85099926E+02   # M_eL
    32     7.85099926E+02   # M_muL
    33     7.85099926E+02   # M_tauL
    34     7.85100020E+02   # M_eR
    35     7.85100020E+02   # M_muR
    36     7.85100019E+02   # M_tauR
    41     7.85099895E+02   # M_q1L
    42     7.85099895E+02   # M_q2L
    43     7.85100318E+02   # M_q3L
    44     7.85099899E+02   # M_uR
    45     7.85099899E+02   # M_cR
    46     7.85100772E+02   # M_tR
    47     7.85099944E+02   # M_dR
    48     7.85099944E+02   # M_sR
    49     7.85099944E+02   # M_bR
# 
# NMSSM SPECIFIC PARAMETERS AT THE GUT SCALE
BLOCK NMSSMRUN Q=  1.74778596E+16 # (GUT SCALE)
     1     6.53581442E-01   # LAMBDA(MGUT,DR_bar)
     2     8.13766496E-02   # KAPPA(MGUT,DR_bar)
     3    -2.24999805E+03   # ALAMBDA
     4    -9.64995044E+02   # AKAPPA
     6     0.00000000E+00   # XIF
     7     0.00000000E+00   # XIS
     8     0.00000000E+00   # MU'
     9     0.00000000E+00   # MS'^2
    10     1.41723919E+06   # MS^2
    12     0.00000000E+00   # M3H^2
# 
# FINE-TUNING parameter d(ln Mz^2)/d(ln PG^2)
# BLOCK FINETUNING
     1     2.88607547E+02   # PG=MHU
     2    -1.22479674E+01   # PG=MHD
     3     1.27711967E+02   # PG=MS
     4    -5.75188335E+01   # PG=M0
     5     0.00000000E+00   # PG=M1
     6     0.00000000E+00   # PG=M2
     7     0.00000000E+00   # PG=M3
     8    -1.63603230E+02   # PG=M12
     9     0.00000000E+00   # PG=ALAMBDA
    10     1.45345331E+00   # PG=AKAPPA
    11    -1.72670773E+02   # PG=A0
    12     0.00000000E+00   # PG=XIF
    13     0.00000000E+00   # PG=XIS
    14     0.00000000E+00   # PG=MUP
    15     0.00000000E+00   # PG=MSP
    16     0.00000000E+00   # PG=M3H
    17    -1.10435412E+02   # PG=LAMBDA
    18    -6.93091869E-01   # PG=KAPPA
    19     1.79030983E+01   # PG=HTOP
    20    -6.21403673E+02   # PG=G0
    21    -1.16996565E+01   # PG=MGUT
    22     2.88607547E+02   # MAX
    23                  1   # IMAX
# 
# REDUCED CROSS SECTIONS AT LHC
BLOCK LHCCROSSSECTIONS
    11     4.70208200E-04   # VBF/VH -> H1 -> tautau
    12     3.08206813E-02   # ggF -> H1 -> tautau
    13     4.62670980E-04   # VBF/VH -> H1 -> bb
    14     1.47782348E-04   # ttH -> H1 -> bb
    15     0.00000000E+00   # VBF/VH -> H1 -> ZZ/WW
    16     0.00000000E+00   # ggF -> H1 -> ZZ/WW
    17     1.96861252E-05   # VBF/VH -> H1 -> gammagamma
    18     1.29036412E-03   # ggF -> H1 -> gammagamma
    21     1.93240235E-01   # VBF/VH -> H2 -> tautau
    22     1.83140040E-01   # ggF -> H2 -> tautau
    23     1.93068652E-01   # VBF/VH -> H2 -> bb
    24     1.92509528E-01   # ttH -> H2 -> bb
    25     1.88134060E-01   # VBF/VH -> H2 -> ZZ/WW
    26     1.78300753E-01   # ggF -> H2 -> ZZ/WW
    27     1.86615498E-01   # VBF/VH -> H2 -> gammagamma
    28     1.76861562E-01   # ggF -> H2 -> gammagamma
    31     8.73769964E-03   # VBF/VH -> H3 -> tautau
    32     2.05627322E+01   # ggF -> H3 -> tautau
    33     5.73277347E-03   # VBF/VH -> H3 -> bb
    34     1.47054687E+01   # ttH -> H3 -> bb
    35     4.11422172E-08   # VBF/VH -> H3 -> ZZ/WW
    36     9.68214097E-05   # ggF -> H3 -> ZZ/WW
    37     1.80686105E-03   # VBF/VH -> H3 -> gammagamma
    38     4.25214890E+00   # ggF -> H3 -> gammagamma
# 
# PARAMETERS OF THE EFFECTIVE LAGRANGIAN IN THE HIGGS SECTOR
BLOCK  EFFECTIVE_COUPLINGS
 X  0.00000000E+00
 A1 -0.24592241E-04
 A2 -0.24452730E-04
 A3  0.00000000E+00
 A4  0.29179759E+01
 A5 -0.30728543E+03
 A6 -0.13573530E+01
 B1  0.00000000E+00
 B2  0.00000000E+00
 L1  0.12405564E+00
 L2  0.26858455E+00
 L3  0.54766272E-01
 L4  0.21071149E-01
 L5 -0.35770521E-03
 L6  0.47992624E-03
 L7  0.12846808E-02
 K1  0.23179998E+00
 K2  0.23048499E+00
 K3  0.48976618E-02
 K4 -0.23706715E-04
 K5  0.00000000E+00
 K6  0.26356138E-01
 K7  0.00000000E+00
 K8  0.00000000E+00
 DELMB  0.33698773E-02
 XVEV -0.42536789E+03
