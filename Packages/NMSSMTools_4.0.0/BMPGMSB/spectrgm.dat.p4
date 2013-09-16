# NMSSMTools OUTPUT IN SLHA FORMAT
# Info about spectrum calculator
BLOCK SPINFO   # Program information
     1   NMSSMTools # Spectrum calculator
     2   4.0.0      # Version number
     8   0          # Higgs mass precision
     3   # Excluded by ee -> hZ, h -> AA -> 4bs
     3   # Landau Pole below MGUT
     3   # Muon magn. mom. more than 2 sigma away
# Input parameters
BLOCK MODSEL
  3  1   # NMSSM particle content
  1  2   # IMOD
 10  0   # ISCAN
  8  0   # Precision for Higgs masses
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
     1     5.00000000E+04   # MSUSYEFF
     2     7.50000000E+06   # MMESS
     3     1.90000000E+00   # TANBETA(MZ)
     4     1.00000000E+00   # SIGMU
     5     2.00000000E+00   # N5
BLOCK EXTPAR
    61     6.00000000E-01   # LAMBDA AT THE SUSY SCALE
#   62     3.96890557E-01   # KAPPA AT THE SUSY SCALE
    64     0.00000000E+00   # ALAMBDA AT THE MES SCALE
#   70    -3.45425386E+05   # MS^2 AT THE MES SCALE
# 
BLOCK MASS   # Mass spectrum 
#  PDG Code     mass             particle 
        25     9.72886037E+01   # lightest neutral scalar
        35     6.88234385E+02   # second neutral scalar
        45     7.55620380E+02   # third neutral scalar
        36     1.27642251E+01   # lightest pseudoscalar
        46     7.15397054E+02   # second pseudoscalar
        37     7.04810113E+02   # charged Higgs
   1000001     8.41119012E+02   #  ~d_L
   2000001     8.01875786E+02   #  ~d_R
   1000002     8.39053062E+02   #  ~u_L
   2000002     8.04898852E+02   #  ~u_R
   1000003     8.41119012E+02   #  ~s_L
   2000003     8.01875786E+02   #  ~s_R
   1000004     8.39053062E+02   #  ~c_L
   2000004     8.04898852E+02   #  ~c_R
   1000005     7.87524110E+02   #  ~b_1
   2000005     7.99411738E+02   #  ~b_2
   1000006     6.70185169E+02   #  ~t_1
   2000006     8.32043345E+02   #  ~t_2
   1000011     2.71718158E+02   #  ~e_L
   2000011     1.37226595E+02   #  ~e_R
   1000012     2.65186333E+02   #  ~nue_L
   1000013     2.71718158E+02   #  ~mu_L
   2000013     1.37226595E+02   #  ~mu_R
   1000014     2.65186333E+02   #  ~numu_L
   1000015     1.36925682E+02   #  ~tau_1
   2000015     2.71836756E+02   #  ~tau_2
   1000016     2.65171156E+02   #  ~nutau_L
   1000021     7.86598868E+02   #  ~g
   1000022     1.31236846E+02   # neutralino(1)
   1000023     2.49044306E+02   # neutralino(2)
   1000025    -5.63923193E+02   # neutralino(3)
   1000035     5.73112305E+02   # neutralino(4)
   1000045     7.47632439E+02   # neutralino(5)
   1000024     2.48491124E+02   # chargino(1)
   1000037     5.75214442E+02   # chargino(2)
# 
# Low energy observables
BLOCK LOWEN
# Exp. 2 Sigma: 3.04E-4 < BR(b -> s gamma) < 4.06E-4:
     1     3.47188106E-04   # BR(b -> s gamma)
    11     3.82049409E-04   # (BR(b -> s gamma)+Theor.Err.)
    12     2.92534523E-04   # (BR(b -> s gamma)-Theor.Err.)
# Exp. 2 Sigma: 4.99E-1 < Delta M_d < 5.15E-1:
     2     6.23575362E-01   # Delta M_d in ps^-1
    21     1.09433119E+00   # Delta M_d +Theor.Err.
    22     1.66445799E-01   # Delta M_d -Theor.Err.
# Exp. 2 Sigma: 1.7633E+1 < Delta Ms < 1.7805E+1:
     3     2.15991002E+01   # Delta M_s in ps^-1
    31     2.86436700E+01   # Delta M_s +Theor.Err.
    32     1.47534460E+01   # Delta M_s -Theor.Err.
# Exp. 2 Sigma: 2.0E-9 < BR(Bs->mu+mu-) < 4.7E-9:
     4     3.34674332E-09   # BR(Bs -> mu+mu-)
    41     5.68438052E-09   # BR(Bs -> mu+mu-)+Theor.Err.
    42     1.62462036E-09   # BR(Bs -> mu+mu-)-Theor.Err.
# Exp. 2 Sigma: 1.07E-4 < BR(B+ > tau+ + nu_tau) < 2.27E-4:
     5     1.31741528E-04   # BR(B+ -> tau+ + nu_tau)
    51     2.63534589E-04   # BR(B+ -> tau+ + nu_tau) + Theor.Err.
    52     5.68564822E-05   # BR(B+ -> tau+ + nu_tau) - Theor.Err.
# 
# BSM contr. to the muon anomalous magn. moment:
     6     2.92935556E-10   # Del_a_mu
    61     5.81957333E-10   # Del_a_mu + Theor.Err.
    62     3.91377886E-12   # Del_a_mu - Theor.Err.
# Minimal Exp.-SM (2 sigma):  8.77306222E-10
# Maximal Exp.-SM (2 sigma):  4.61144414E-09
# 
BLOCK HMIX Q=  7.09468602E+02 # (STOP/SBOTTOM MASSES)
     1     5.54365349E+02   # MUEFF
     2     1.90000374E+00   # TAN(BETA)
     3     2.43876894E+02   # V(Q)
     4     5.11938644E+05   # MA^2
     5     1.09681551E+04   # MP^2
# 
# 3*3 Higgs mixing
BLOCK NMHMIX
  1  1     4.58023707E-01   # S_(1,1)
  1  2     8.83732857E-01   # S_(1,2)
  1  3    -9.60756042E-02   # S_(1,3)
  2  1     7.82762156E-01   # S_(2,1)
  2  2    -3.49730452E-01   # S_(2,2)
  2  3     5.14754329E-01   # S_(2,3)
  3  1    -4.21304749E-01   # S_(3,1)
  3  2     3.10974033E-01   # S_(3,2)
  3  3     8.51937474E-01   # S_(3,3)
# 
# 3*3 Pseudoscalar Higgs mixing
BLOCK NMAMIX
  1  1     1.28450179E-01   # P_(1,1)
  1  2     6.76053573E-02   # P_(1,2)
  1  3     9.89408948E-01   # P_(1,3)
  2  1     8.75546008E-01   # P_(2,1)
  2  2     4.60813688E-01   # P_(2,2)
  2  3    -1.45154858E-01   # P_(2,3)
# 
# 3rd generation sfermion mixing
BLOCK STOPMIX  # Stop mixing matrix
  1  1     4.35956663E-01   # Rst_(1,1)
  1  2     8.99967659E-01   # Rst_(1,2)
  2  1    -8.99967659E-01   # Rst_(2,1)
  2  2     4.35956663E-01   # Rst_(2,2)
BLOCK SBOTMIX  # Sbottom mixing matrix
  1  1     9.80904380E-01   # Rsb_(1,1)
  1  2     1.94490611E-01   # Rsb_(1,2)
  2  1    -1.94490611E-01   # Rsb_(2,1)
  2  2     9.80904380E-01   # Rsb_(2,2)
BLOCK STAUMIX  # Stau mixing matrix
  1  1     3.54636387E-02   # Rsl_(1,1)
  1  2     9.99370967E-01   # Rsl_(1,2)
  2  1    -9.99370967E-01   # Rsl_(2,1)
  2  2     3.54636387E-02   # Rsl_(2,2)
# 
# Gaugino-Higgsino mixing
BLOCK NMNMIX  # 5*5 Neutralino Mixing Matrix
  1  1     9.92383262E-01   # N_(1,1)
  1  2    -5.72114944E-02   # N_(1,2)
  1  3     8.92362011E-02   # N_(1,3)
  1  4    -6.21630988E-02   # N_(1,4)
  1  5     8.65768392E-03   # N_(1,5)
  2  1     8.17087045E-02   # N_(2,1)
  2  2     9.68298345E-01   # N_(2,2)
  2  3    -1.82996209E-01   # N_(2,3)
  2  4     1.47760713E-01   # N_(2,4)
  2  5    -2.00290379E-02   # N_(2,5)
  3  1    -1.81842781E-02   # N_(3,1)
  3  2     2.86341120E-02   # N_(3,2)
  3  3     7.04260307E-01   # N_(3,3)
  3  4     7.04987882E-01   # N_(3,4)
  3  5     7.65436240E-02   # N_(3,5)
  4  1    -8.98446014E-02   # N_(4,1)
  4  2     2.40328246E-01   # N_(4,2)
  4  3     6.57816284E-01   # N_(4,3)
  4  4    -6.87596604E-01   # N_(4,4)
  4  5     1.69289481E-01   # N_(4,5)
  5  1     9.81978935E-03   # N_(5,1)
  5  2    -2.34004833E-02   # N_(5,2)
  5  3    -1.72755447E-01   # N_(5,3)
  5  4     6.71231496E-02   # N_(5,4)
  5  5     9.82347203E-01   # N_(5,5)
# 
BLOCK UMIX  # Chargino U Mixing Matrix
  1  1     9.62779900E-01   # U_(1,1)
  1  2    -2.70286635E-01   # U_(1,2)
  2  1     2.70286635E-01   # U_(2,1)
  2  2     9.62779900E-01   # U_(2,2)
# 
BLOCK VMIX  # Chargino V Mixing Matrix
  1  1     9.76809195E-01   # V_(1,1)
  1  2    -2.14111646E-01   # V_(1,2)
  2  1     2.14111646E-01   # V_(2,1)
  2  2     9.76809195E-01   # V_(2,2)
# 
# Higgs reduced couplings
# (as compared to a SM Higgs with same mass)
BLOCK REDCOUP
# H1
  1  1     9.98660480E-01   # U-type fermions
  1  2     9.83418605E-01   # D-type fermions
  1  3     9.95354216E-01   # W,Z bosons
  1  4     1.01431219E+00   # Gluons
  1  5     9.99660953E-01   # Photons
# H2
  2  1    -3.95212171E-01   # U-type fermions
  2  2     1.68066162E+00   # D-type fermions
  2  3     5.50858321E-02   # W,Z bosons
  2  4     3.98621891E-01   # Gluons
  2  5     2.33714292E+00   # Photons
# H3
  3  1     3.51415560E-01   # U-type fermions
  3  2    -9.04579658E-01   # D-type fermions
  3  3     7.89654044E-02   # W,Z bosons
  3  4     3.52754236E-01   # Gluons
  3  5     1.08034077E+00   # Photons
# A1
  4  1     7.63972937E-02   # U-type fermions
  4  2     2.75794230E-01   # D-type fermions
  4  3     0.00000000E+00   # W,Z bosons
  4  4     2.11215898E-01   # Gluons
  4  5     2.36099079E-01   # Photons
# A2
  5  1     5.20741552E-01   # U-type fermions
  5  2     1.87987700E+00   # D-type fermions
  5  3     0.00000000E+00   # W,Z bosons
  5  4     5.24105714E-01   # Gluons
  5  5     6.73187266E-01   # Photons
# 
# GAUGE AND YUKAWA COUPLINGS AT THE SUSY SCALE
BLOCK GAUGE Q=  7.94837786E+02 # (SUSY SCALE)
     1     3.62450067E-01   # g1(Q,DR_bar)
     2     6.44550633E-01   # g2(Q,DR_bar)
     3     1.07402573E+00   # g3(Q,DR_bar)
BLOCK YU Q=  7.94837786E+02 # (SUSY SCALE)
  3  3     9.69438571E-01   # HTOP(Q,DR_bar)
BLOCK YD Q=  7.94837786E+02 # (SUSY SCALE)
  3  3     3.02195368E-02   # HBOT(Q,DR_bar)
BLOCK YE Q=  7.94837786E+02 # (SUSY SCALE)
  3  3     2.14814073E-02   # HTAU(Q,DR_bar)
# 
# SOFT TRILINEAR COUPLINGS AT THE SUSY SCALE
BLOCK AU Q=  7.94837786E+02 # (SUSY SCALE)
  3  3    -3.14234444E+02   # ATOP
BLOCK AD Q=  7.94837786E+02 # (SUSY SCALE)
  3  3    -3.98573961E+02   # ABOT
BLOCK AE Q=  7.94837786E+02 # (SUSY SCALE)
  2  2    -4.67041938E+01   # AMUON
  3  3    -4.65733037E+01   # ATAU
# 
# SOFT MASSES AT THE SUSY SCALE
BLOCK MSOFT Q=  7.94837786E+02 # (SUSY SCALE)
     1     1.36914753E+02   # M1
     2     2.58025535E+02   # M2
     3     7.22856549E+02   # M3
    21     8.68723403E+04   # M_HD^2
    22    -1.84726048E+05   # M_HU^2
    31     2.69449865E+02   # M_eL
    32     2.69449865E+02   # M_muL
    33     2.69434928E+02   # M_tauL
    34     1.33338770E+02   # M_eR
    35     1.33338770E+02   # M_muR
    36     1.33277643E+02   # M_tauR
    41     8.13663868E+02   # M_q1L
    42     8.13663868E+02   # M_q2L
    43     7.61530481E+02   # M_q3L
    44     7.77471903E+02   # M_uR
    45     7.77471903E+02   # M_cR
    46     6.60965922E+02   # M_tR
    47     7.73633044E+02   # M_dR
    48     7.73633044E+02   # M_sR
    49     7.73522257E+02   # M_bR
# 
# NMSSM SPECIFIC PARAMETERS AT THE SUSY SCALE
BLOCK NMSSMRUN Q=  7.94837786E+02 # (SUSY SCALE)
     1     6.00000000E-01   # LAMBDA(Q,DR_bar)
     2     3.96890557E-01   # KAPPA(Q,DR_bar)
     3     1.30831520E+01   # ALAMBDA
     4     7.35415913E-01   # AKAPPA
     5     5.55242997E+02   # MUEFF
     6     0.00000000E+00   # XIF
     7     0.00000000E+00   # XIS
     8     0.00000000E+00   # MU'
     9     0.00000000E+00   # MS'^2
    10    -2.74752791E+05   # MS^2
    12     0.00000000E+00   # M3H^2
# 
# GAUGE AND YUKAWA COUPLINGS AT THE MES SCALE
BLOCK GAUGE Q=  7.50000000E+06 # (MES SCALE)
     1     3.97869496E-01   # g1(MMES,DR_bar)
     2     6.63676171E-01   # g2(MMES,DR_bar)
     3     9.11352480E-01   # g3(MMES,DR_bar)
BLOCK YU Q=  7.50000000E+06 # (MES SCALE)
  3  3     9.16727732E-01   # HTOP(MMES,DR_bar)
BLOCK YD Q=  7.50000000E+06 # (MES SCALE)
  3  3     2.22443932E-02   # HBOT(MMES,DR_bar)
BLOCK YE Q=  7.50000000E+06 # (MES SCALE)
  3  3     1.99570411E-02   # HTAU(MMES,DR_bar)
# 
# SOFT TRILINEAR COUPLINGS AT THE MES SCALE
BLOCK AU Q=  7.50000000E+06 # (MES SCALE)
  3  3     2.02036463E-06   # ATOP
BLOCK AD Q=  7.50000000E+06 # (MES SCALE)
  3  3     2.10350099E-06   # ABOT
BLOCK AE Q=  7.50000000E+06 # (MES SCALE)
  2  2     1.89740053E-07   # AMUON
  3  3     1.89568451E-07   # ATAU
# 
# SOFT MASSES SQUARED AT THE MES SCALE
BLOCK MSOFTMES Q=  7.50000000E+06 # (MES SCALE)
     1     1.67075791E+02   # M1
     2     2.78930457E+02   # M2
     3     5.25964276E+02   # M3
    21     6.25379491E+04   # M_HD^2
    22     6.25379504E+04   # M_HU^2
    31     2.50075886E+02   # M_eL
    32     2.50075886E+02   # M_muL
    33     2.50075886E+02   # M_tauL
    34     1.29415472E+02   # M_eR
    35     1.29415472E+02   # M_muR
    36     1.29415472E+02   # M_tauR
    41     6.53958949E+02   # M_q1L
    42     6.53958949E+02   # M_q2L
    43     6.53958949E+02   # M_q3L
    44     6.13424758E+02   # M_uR
    45     6.13424758E+02   # M_cR
    46     6.13424759E+02   # M_tR
    47     6.08857246E+02   # M_dR
    48     6.08857246E+02   # M_sR
    49     6.08857246E+02   # M_bR
# 
# NMSSM SPECIFIC PARAMETERS AT THE MES SCALE
BLOCK NMSSMRUN Q=  7.50000000E+06 # (MES SCALE)
     1     7.29761360E-01   # LAMBDA(MMES,DR_bar)
     2     4.93131150E-01   # KAPPA(MMES,DR_bar)
     3     7.74230535E-08   # ALAMBDA
     4    -3.76781551E-08   # AKAPPA
     6     0.00000000E+00   # XIF
     7     0.00000000E+00   # XIS
     8     0.00000000E+00   # MU'
     9     0.00000000E+00   # MS'^2
    10    -3.45425386E+05   # MS^2
    11     0.00000000E+00   # DH
    12     0.00000000E+00   # M3H^2
# 
# GAUGE COUPLINGS AT THE SCALE OF THE LANDAU SINGULARITY
BLOCK GAUGE MLANDAU=  4.36321281E+13 # (LANDAU SCALE)
     1     6.91105869E-01   # g1(MLANDAU), (GUT normalization)
     2     7.71863660E-01   # g2(MLANDAU)
     3     8.44768309E-01   # g3(MLANDAU)
# YUKAWA COUPLINGS AT THE LANDAU SCALE
BLOCK YU MLANDAU=  4.36321281E+13 # (MLANDAU)
  3  3     9.36117234E-01   # HTOP(MLANDAU)
BLOCK YD MLANDAU=  4.36321281E+13 # (MLANDAU)
  3  3     1.53304665E-02   # HBOT(MLANDAU)
BLOCK YE MLANDAU=  4.36321281E+13 # (MLANDAU)
  3  3     1.81472715E-02   # HTAU(MLANDAU)
# NMSSM SPECIFIC PARAMETERS AT THE LANDAU SCALE
BLOCK NMSSMRUN MLANDAU=  4.36321281E+13 # (LNDAU SCALE)
    61     3.24011682E+02   # LAMBDA(MLANDAU)
    62     5.93866712E+14   # KAPPA(MLANDAU)
# 
# FINE-TUNING parameter d(ln Mz^2)/d(ln PG^2)
# BLOCK FINETUNING
     1    -3.16191104E+01   # PM=MSUSYEFF
     2     3.40112118E+01   # PM=MS
     3     0.00000000E+00   # PM=DELTA_H
     4     0.00000000E+00   # PM=ALAMBDA
     5     0.00000000E+00   # PM=XIF
     6     0.00000000E+00   # PM=XIS
     7     0.00000000E+00   # PM=MUP
     8     0.00000000E+00   # PM=MSP
     9     5.42020111E+01   # PM=LAMBDA
    10    -4.88463056E+01   # PM=KAPPA
    11    -3.90108128E+01   # PM=HTOP
    12     1.36328571E+00   # PM=G1
    13     3.03206032E+01   # PM=G2
    14    -1.28938212E+02   # PM=G3
    15    -1.79630517E+00   # PM=MMES
    16     1.28938212E+02   # MAX
    17                 14   # IMAX
# 
# REDUCED CROSS SECTIONS AT LHC
BLOCK LHCCROSSSECTIONS
    11     2.57075718E-01   # VBF/VH -> H1 -> tautau
    12     2.66961743E-01   # ggF -> H1 -> tautau
    13     2.57289659E-01   # VBF/VH -> H1 -> bb
    14     2.59001774E-01   # ttH -> H1 -> bb
    15     2.63353769E-01   # VBF/VH -> H1 -> ZZ/WW
    16     2.73481220E-01   # ggF -> H1 -> ZZ/WW
    17     2.65653185E-01   # VBF/VH -> H1 -> gammagamma
    18     2.75869063E-01   # ggF -> H1 -> gammagamma
    21     2.02071716E-01   # VBF/VH -> H2 -> tautau
    22     1.05815184E+01   # ggF -> H2 -> tautau
    23     1.32534407E-01   # VBF/VH -> H2 -> bb
    24     6.82196412E+00   # ttH -> H2 -> bb
    25     2.17082488E-04   # VBF/VH -> H2 -> ZZ/WW
    26     1.13675599E-02   # ggF -> H2 -> ZZ/WW
    27     3.90787676E-01   # VBF/VH -> H2 -> gammagamma
    28     2.04636605E+01   # ggF -> H2 -> gammagamma
    31     1.52805915E-01   # VBF/VH -> H3 -> tautau
    32     3.04937556E+00   # ggF -> H3 -> tautau
    33     1.04767222E-01   # VBF/VH -> H3 -> bb
    34     2.07488330E+00   # ttH -> H3 -> bb
    35     1.16444843E-03   # VBF/VH -> H3 -> ZZ/WW
    36     2.32375859E-02   # ggF -> H3 -> ZZ/WW
    37     2.17968340E-01   # VBF/VH -> H3 -> gammagamma
    38     4.34974868E+00   # ggF -> H3 -> gammagamma
# 
# PARAMETERS OF THE EFFECTIVE LAGRANGIAN IN THE HIGGS SECTOR
BLOCK  EFFECTIVE_COUPLINGS
 X  0.00000000E+00
 A1  0.00000000E+00
 A2  0.00000000E+00
 A3  0.00000000E+00
 A4  0.91270979E-01
 A5  0.79941020E+01
 A6 -0.62239009E+00
 B1  0.00000000E+00
 B2  0.00000000E+00
 L1  0.12226389E+00
 L2  0.23589618E+00
 L3  0.77927635E-01
 L4  0.15814986E+00
 L5 -0.33767552E-03
 L6  0.17793538E-03
 L7 -0.51875698E-02
 K1  0.35197105E+00
 K2  0.35197244E+00
 K3  0.15684124E+00
 K4  0.00000000E+00
 K5  0.00000000E+00
 K6  0.23190570E+00
 K7  0.00000000E+00
 K8  0.00000000E+00
 DELMB  0.11982015E-01
 XVEV  0.92713474E+03
